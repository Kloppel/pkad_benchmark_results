#!/usr/bin/env python3
"""queue_all.py — Submit all PKAD benchmark jobs to the SLURM cluster.

Reads pkad_data/job_manifest.tsv (produced by prepare_benchmark.py) and
submits a SLURM job array covering all entries (crystal + AF3).

Usage
-----
    # Preview — print the sbatch command without submitting
    python3 queue_all.py --dry-run

    # Submit to cluster
    python3 queue_all.py

    # Resubmit only structures that have not yet produced a _results.txt
    python3 queue_all.py --resubmit-missing \
        --results-dir results/pkad_benchmark_20260619

    # After jobs finish: generate CHARMM build scripts (does NOT execute them)
    python3 queue_all.py --generate-charmm-scripts \
        --results-root results/pkad_benchmark_20260620 \
        --charmm-out-dir pkad_data/charmm_build_scripts

    # Combine: submit then generate scripts for already-finished runs
    python3 queue_all.py --generate-charmm-scripts \
        --results-root results/pkad_benchmark_20260620
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Optional

import csv

HERE = Path(__file__).parent
DEFAULT_MANIFEST = HERE / "pkad_data" / "job_manifest.tsv"
DEFAULT_JOB_SCRIPT = HERE / "pkad_benchmark.job"
BENCHMARKING_DIR = HERE.parent
REPO_ROOT = HERE.parents[1]

# Structures that require a 24 h TAPBS internal timeout (computed timeout too short).
TAPBS_TIMEOUT_STRUCTS: frozenset[str] = frozenset({
    "1bvi", "1i0e", "2pyo", "1q74", "1nlx", "1m56", "1tgu", "1bmf", "1qk1",
})

# Structures where the salt-bridge-opener PACs (-10 / +20) abort because no
# salt bridges are found.  Run only the pH-7 h_min PAC for these.
PH7_ONLY_STRUCTS: frozenset[str] = frozenset({
    "1a91", "1ans", "1beg", "1cvo", "1gb1", "1hic", "1jv8", "1m2c",
    "1poh", "1va3", "1vii", "1yrf", "2bus", "2gb1",
})


# ---------------------------------------------------------------------------
# Manifest loading
# ---------------------------------------------------------------------------

def load_manifest(path: Path) -> list[dict]:
    if not path.exists():
        print(f"[ERROR] Manifest not found: {path}", file=sys.stderr)
        print("        Run prepare_benchmark.py first.", file=sys.stderr)
        sys.exit(1)
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    n_crystal = sum(1 for r in rows if r.get("structure_type") == "crystal")
    n_af3 = sum(1 for r in rows if r.get("structure_type") == "af3")
    print(f"[manifest] {len(rows)} jobs  ({n_crystal} crystal, {n_af3} AF3)")
    return rows


# ---------------------------------------------------------------------------
# SLURM submission
# ---------------------------------------------------------------------------

def find_missing_indices(manifest: list[dict], results_dir: Path) -> list[int]:
    """Return 0-based manifest row indices for runs without a *_results.txt file."""
    missing = []
    for i, row in enumerate(manifest):
        run_name = (row["pdb_id"] if row.get("structure_type") == "crystal"
                    else f"{row['pdb_id']}_af3")
        if not (results_dir / f"{run_name}_results.txt").exists():
            missing.append(i)
    return missing


def partition_indices(
    manifest: list[dict],
    indices: list[int],
) -> tuple[list[int], list[int], list[int]]:
    """Split *indices* into (regular, tapbs_timeout, ph7_only) groups.

    Crystal entries whose pdb_id is in TAPBS_TIMEOUT_STRUCTS or PH7_ONLY_STRUCTS
    are routed to the matching special group; all others go to regular.
    AF3 entries are never routed to a special group (the special sets are
    crystal-only based on experimental evidence).
    """
    regular, tapbs_timeout, ph7_only = [], [], []
    for i in indices:
        row = manifest[i]
        pdb_id = row["pdb_id"]
        stype = row.get("structure_type", "crystal")
        if stype == "crystal" and pdb_id in TAPBS_TIMEOUT_STRUCTS:
            tapbs_timeout.append(i)
        elif stype == "crystal" and pdb_id in PH7_ONLY_STRUCTS:
            ph7_only.append(i)
        else:
            regular.append(i)
    return regular, tapbs_timeout, ph7_only


def submit_job_array(
    job_script: Path,
    n_tasks: int,
    max_concurrent: int = 50,
    dry_run: bool = False,
    extra_sbatch_args: Optional[list[str]] = None,
    array_indices: Optional[list[int]] = None,
) -> Optional[str]:
    """Submit a SLURM job array. Returns the job ID string, or None on dry-run.

    If array_indices is given, only those specific manifest rows are submitted
    (using SLURM's comma-separated index syntax).  Otherwise submits 0..n_tasks-1.
    """
    if array_indices is not None:
        array_spec = ",".join(str(i) for i in sorted(array_indices))
        array_spec += f"%{max_concurrent}"
    else:
        array_spec = f"0-{n_tasks - 1}%{max_concurrent}"
    cmd = [
        "sbatch",
        f"--array={array_spec}",
    ] + (extra_sbatch_args or []) + [str(job_script)]

    print(f"\nSLURM command: {' '.join(cmd)}")

    if dry_run:
        print("[dry-run] Not submitting.")
        return None

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] sbatch failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    output = result.stdout.strip()
    print(output)

    # Parse "Submitted batch job 12345"
    m = re.search(r"(\d+)$", output)
    job_id = m.group(1) if m else None
    if job_id:
        print(f"[submit] Job ID: {job_id}  ({n_tasks} tasks, max {max_concurrent} concurrent)")
        # Save for reference
        submitted_log = HERE / "pkad_data" / "submitted_jobs.txt"
        with open(submitted_log, "a") as fh:
            fh.write(f"{job_id}\t{n_tasks} tasks\t{job_script.name}\n")
        print(f"[write] {submitted_log}")
    return job_id


# ---------------------------------------------------------------------------
# CHARMM build script generation
# ---------------------------------------------------------------------------

def _find_protonated_pdb(results_root: Path, run_name: str) -> Optional[Path]:
    """Return the CHARMM-output protonated PDB from the KB3 result tree.

    KB3 writes: results_root/{run_name}/{run_name}_karlsberg/
                  initial_modelling/{run_name}_stripped_out.pdb
    """
    base = results_root / run_name / f"{run_name}_karlsberg" / "initial_modelling"
    candidate = base / f"{run_name}_stripped_out.pdb"
    if candidate.exists():
        return candidate
    # Also check without _stripped suffix
    candidate2 = base / f"{run_name}_out.pdb"
    if candidate2.exists():
        return candidate2
    return None


def _extract_disu_patches(inp_dir: Path, run_name: str) -> list[tuple[str, str]]:
    """Parse existing CHARMM .inp files for PATCH DISU lines.

    Returns list of (segid_res1, segid_res2) pairs, e.g. ("PROA 6", "PROA 38").
    """
    patches = []
    for inp_file in inp_dir.glob("*.inp"):
        text = inp_file.read_text(errors="replace")
        for line in text.splitlines():
            s = line.strip().upper()
            if s.startswith("PATCH") and "DISU" in s:
                # PATCH DISU PROA 6 PROA 38
                parts = s.split()
                # parts: ['PATCH', 'DISU', seg1, resid1, seg2, resid2]
                if len(parts) >= 6:
                    patches.append((f"{parts[2]} {parts[3]}", f"{parts[4]} {parts[5]}"))
    return patches


def generate_charmm_scripts(
    manifest: list[dict],
    results_root: Path,
    top_dir: Path,
    par_dir: Path,
    out_dir: Path,
) -> None:
    """Generate CHARMM .inp build scripts for completed KB3 runs.

    For each manifest row where the KB3-output protonated PDB exists, emits
    out_dir/{run_name}_build.inp that reads topology/parameters, builds PSF,
    and writes built PDB + PSF.  Does NOT execute the scripts.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # Locate topology and parameter files
    top_file = top_dir / "top_all36_prot.rtf"
    if not top_file.exists():
        # Try to find any RTF file
        rtf_files = list(top_dir.glob("top_all36*.rtf"))
        top_file = rtf_files[0] if rtf_files else top_dir / "top_all36_prot.rtf"

    par_file = par_dir / "par_all36m_prot.prm"
    if not par_file.exists():
        prm_files = list(par_dir.glob("par_all36*.prm"))
        par_file = prm_files[0] if prm_files else par_dir / "par_all36m_prot.prm"

    generated = 0
    skipped = 0

    for row in manifest:
        pdb_id = row["pdb_id"]
        stype = row.get("structure_type", "crystal")
        run_name = pdb_id if stype == "crystal" else f"{pdb_id}_af3"

        prot_pdb = _find_protonated_pdb(results_root, run_name)
        if prot_pdb is None:
            skipped += 1
            continue

        # Determine segment name from PDB (default PROA)
        segname = "PROA"

        # Extract DISU patches from existing CHARMM input
        inp_dir = prot_pdb.parent
        disu_patches = _extract_disu_patches(inp_dir, run_name)

        # Build CHARMM script content
        lines = [
            f"* CHARMM build script for {run_name}",
            f"* Generated by queue_all.py --generate-charmm-scripts",
            f"* Input: {prot_pdb}",
            "*",
            "",
            "bomlev -2",
            "",
            f"read rtf  card name {top_file}",
            f"read param card name {par_file}",
            "",
            f"read sequence pdb name {prot_pdb}",
            f"generate {segname} first NTER last CTER setup",
            "",
        ]

        for seg1, seg2 in disu_patches:
            lines.append(f"patch DISU {seg1} {seg2}")
        if disu_patches:
            lines.append("")

        built_pdb = out_dir / f"{run_name}_built.pdb"
        built_psf = out_dir / f"{run_name}_built.psf"
        lines += [
            "auto angles dihe",
            "",
            f"write coor pdb  name {built_pdb}",
            f"write psf  card name {built_psf}",
            "",
            "stop",
        ]

        out_file = out_dir / f"{run_name}_build.inp"
        out_file.write_text("\n".join(lines) + "\n")
        generated += 1

    print(f"[charmm scripts] generated={generated}  skipped(no result)={skipped}")
    print(f"[write] {out_dir}/")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                        help="Job manifest TSV (default: pkad_data/job_manifest.tsv)")
    parser.add_argument("--job-script", type=Path, default=DEFAULT_JOB_SCRIPT,
                        help="SLURM job script (default: pkad_benchmark.job)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print sbatch command without submitting.")
    parser.add_argument("--max-concurrent", type=int, default=50,
                        help="Max simultaneous SLURM array tasks (default: 50).")
    parser.add_argument("--generate-charmm-scripts", action="store_true",
                        help="Generate CHARMM .inp build scripts for finished runs.")
    parser.add_argument("--results-root", type=Path, default=None,
                        help="Root directory of completed KB3 results "
                             "(required for --generate-charmm-scripts).")
    parser.add_argument("--top-dir", type=Path,
                        default=REPO_ROOT / "benchmarking" / "test" / "top_c36",
                        help="CHARMM topology directory.")
    parser.add_argument("--par-dir", type=Path,
                        default=REPO_ROOT / "benchmarking" / "test" / "par_c36",
                        help="CHARMM parameter directory.")
    parser.add_argument("--charmm-out-dir", type=Path,
                        default=HERE / "pkad_data" / "charmm_build_scripts",
                        help="Output directory for CHARMM .inp files.")
    parser.add_argument("--resubmit-missing", action="store_true",
                        help="Only submit jobs for structures that have no *_results.txt "
                             "in --results-dir yet (skips already-converged runs).")
    parser.add_argument("--results-dir", type=Path, default=None,
                        help="Results directory to check for completed runs "
                             "(required with --resubmit-missing).")
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    if len(manifest) == 0:
        print("[WARN] Manifest is empty — no jobs to submit.")
        print("       Run download_pkad_data.py and prepare_benchmark.py first.")
        if not args.generate_charmm_scripts:
            sys.exit(0)

    # Determine which array indices to submit
    if args.resubmit_missing:
        if args.results_dir is None:
            print("[ERROR] --results-dir is required with --resubmit-missing", file=sys.stderr)
            sys.exit(1)
        if not args.results_dir.exists():
            print(f"[ERROR] Results directory not found: {args.results_dir}", file=sys.stderr)
            sys.exit(1)
        candidate_indices = find_missing_indices(manifest, args.results_dir)
        n_done = len(manifest) - len(candidate_indices)
        print(f"[resubmit] {n_done}/{len(manifest)} already converged — "
              f"{len(candidate_indices)} missing jobs to submit")
        if not candidate_indices:
            print("[resubmit] All jobs converged — nothing to submit.")
            sys.exit(0)
    else:
        candidate_indices = list(range(len(manifest)))

    # Split into regular / special groups so each group gets the right env vars.
    regular, tapbs_timeout, ph7_only = partition_indices(manifest, candidate_indices)
    print(f"[partition] regular={len(regular)}  "
          f"tapbs-timeout={len(tapbs_timeout)}  ph7-only={len(ph7_only)}")

    # Submit job array(s) (unless only generating CHARMM scripts)
    if not (args.generate_charmm_scripts and args.results_root and not args.dry_run):
        if not args.generate_charmm_scripts or args.dry_run:
            if regular:
                print(f"\n[submit] Regular structures ({len(regular)} jobs)")
                submit_job_array(
                    job_script=args.job_script,
                    n_tasks=len(manifest),
                    max_concurrent=args.max_concurrent,
                    dry_run=args.dry_run,
                    array_indices=regular,
                )
            if tapbs_timeout:
                print(f"\n[submit] TAPBS-timeout structures ({len(tapbs_timeout)} jobs,"
                      f" KB3_TAPBS_TIMEOUT=86400)")
                submit_job_array(
                    job_script=args.job_script,
                    n_tasks=len(manifest),
                    max_concurrent=min(args.max_concurrent, len(tapbs_timeout)),
                    dry_run=args.dry_run,
                    array_indices=tapbs_timeout,
                    extra_sbatch_args=["--export=ALL,KB3_TAPBS_TIMEOUT=86400"],
                )
            if ph7_only:
                print(f"\n[submit] No-salt-bridge structures ({len(ph7_only)} jobs,"
                      f" KB3_PH7_ONLY=1)")
                submit_job_array(
                    job_script=args.job_script,
                    n_tasks=len(manifest),
                    max_concurrent=min(args.max_concurrent, len(ph7_only)),
                    dry_run=args.dry_run,
                    array_indices=ph7_only,
                    extra_sbatch_args=["--export=ALL,KB3_PH7_ONLY=1"],
                )

    # Generate CHARMM build scripts for completed runs
    if args.generate_charmm_scripts:
        if args.results_root is None:
            print("[ERROR] --results-root is required for --generate-charmm-scripts",
                  file=sys.stderr)
            sys.exit(1)
        print(f"\nGenerating CHARMM build scripts from {args.results_root} …")
        generate_charmm_scripts(
            manifest=manifest,
            results_root=args.results_root,
            top_dir=args.top_dir,
            par_dir=args.par_dir,
            out_dir=args.charmm_out_dir,
        )


if __name__ == "__main__":
    main()
