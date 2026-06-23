"""
Generate a comprehensive overview table of all PKAD benchmark pKa results.

Reads:
  - pkad_data/experimental_pkas.tsv     — all experimental pKa values + metadata
  - pkad_data/super_table.csv           — full multi-dataset metadata, conflict flags, methods
  - pkad_data/job_manifest.tsv          — which structures were queued
  - pkad_data/benchmark_summary_*.tsv   — pre-computed benchmark summaries
  - results/*/[pdbid]_results.txt       — KB3 calculated pKa values per run
  - logs/pkad_*.out                     — SLURM environment and timing info
  - results/*/[pdbid]_karlsberg/        — convergence traces and setup log

Output:
  - overview/pkad_overview_table.csv    — one row per (pdb_id, chain, resname, resid, run_name)
  - overview/pkad_run_metadata.csv      — one row per (pdb_id, run_name) with job-level metadata
  - overview/pkad_statistics.txt        — per-dataset and overall RMSE/MAE/N statistics

Usage:
    python generate_overview_table.py [--results-root RESULTS_ROOT] [--output-dir OUTPUT_DIR]
"""

import argparse
import csv
import json
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths (relative to repo root)
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).parent.parent


def _repo(rel: str) -> Path:
    return REPO_ROOT / rel


# ---------------------------------------------------------------------------
# Residue name normalisation (KB3 internal names → standard PDB names)
# ---------------------------------------------------------------------------

KB3_TO_STANDARD = {
    "EPP": "GLU",
    "DPP": "ASP",
    "HSP": "HIS",
    "HSE": "HIS",
    "HSD": "HIS",
    "NTE": "NTE",   # N-terminus
    "CTE": "CTE",   # C-terminus
}

STANDARD_TO_KB3 = {v: k for k, v in KB3_TO_STANDARD.items() if k != v}


def norm_resname(name: str) -> str:
    return KB3_TO_STANDARD.get(name.upper(), name.upper())


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_results_txt(path: Path) -> dict:
    """
    Parse a KB3 *_results.txt file.

    Returns a dict keyed by (chain, resname, resid) -> pKa_value (float).
    Lines look like:   '       LYS-1_A: 11.72 '
    """
    pkas = {}
    pattern = re.compile(
        r"^\s*(\w+)-(\d+)_([A-Z]):\s*([\d\.\-]+)\s*$"
    )
    try:
        with open(path) as fh:
            for line in fh:
                m = pattern.match(line)
                if m:
                    resname, resid, chain, value = m.groups()
                    try:
                        pkas[(chain, resname, int(resid))] = float(value)
                    except ValueError:
                        pass
    except OSError:
        pass
    return pkas


def parse_experimental_pkas(path: Path) -> dict:
    """
    Returns dict keyed by (pdb_id, chain, resname, resid) -> expt_pka.
    """
    data = {}
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                raw_resid = row["resid"].split()[0].strip()
                key = (
                    row["pdb_id"].lower(),
                    row.get("chain", ""),
                    row["resname"].upper(),
                    int(raw_resid),
                )
                try:
                    data[key] = float(row["expt_pka"])
                except (ValueError, KeyError):
                    pass
    except OSError:
        pass
    return data


def parse_super_table(path: Path) -> dict:
    """
    Returns dict keyed by (pdb_id, chain, resname, resid) -> row dict.
    """
    data = {}
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                raw_resid = row["resid"].split()[0].strip() if row.get("resid") else "0"
                key = (
                    row["pdb_id"].lower(),
                    row.get("chain", ""),
                    row["resname"].upper(),
                    int(raw_resid) if raw_resid else 0,
                )
                data[key] = row
    except OSError:
        pass
    return data


def parse_job_manifest(path: Path) -> dict:
    """
    Returns dict keyed by pdb_id (lower) -> row dict.
    """
    data = {}
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                pid = row["pdb_id"].lower()
                data[pid] = row
    except OSError:
        pass
    return data


def parse_slurm_logs(logs_dir: Path) -> dict:
    """
    Parse SLURM .out files for per-structure metadata.

    Returns dict keyed by pdb_id -> {
        'slurm_job_id', 'task_id', 'datasets', 'structure_type',
        'pdb_path', 'conda_env', 'charmm_path', 'pac_workers',
        'tapbs_timeout', 'ph7_only', 'date_utc',
        'topology_files', 'parameter_files',
        'n_atoms_retained', 'pac_timeout_budget_s',
        'wall_seconds', 'protocol'
    }
    """
    results = {}
    if not logs_dir.is_dir():
        return results

    header_re = {
        "pdb": re.compile(r"^PDB\s+:\s+(\S+)\s+\((\w+)\).*datasets:\s*(.*)$"),
        "date": re.compile(r"^Date UTC\s+:\s+(.+)$"),
        "conda": re.compile(r"^Conda env\s+:\s+(\S+)$"),
        "charmm": re.compile(r"^CHARMM\s+:\s+(.+)$"),
        "pac_workers": re.compile(r"^PAC workers:\s+(\d+)$"),
        "results_root": re.compile(r"^Results root:\s+(.+)$"),
        "task_id": re.compile(r"^PKAD benchmark.*task\s+(\d+)$"),
    }
    topology_re = re.compile(r"Topology\s+:\s+\[(.+)\]")
    params_re = re.compile(r"Parameters\s+:\s+\[(.+)\]")
    atoms_re = re.compile(r"Retained (\d+) atoms")
    timeout_re = re.compile(r"PAC timeout budget:\s+([\d]+)s")
    protocol_re = re.compile(r"Protocol:\s+\[(.+)\]")
    tapbs_timeout_re = re.compile(r"KB3_TAPBS_TIMEOUT[=:]\s*(\d+)")
    ph7_only_re = re.compile(r"KB3_PH7_ONLY[=:]\s*(\d)")

    for log_file in sorted(logs_dir.glob("pkad_*.out")):
        m = re.match(r"pkad_(\d+)_(\d+)\.out", log_file.name)
        if not m:
            continue
        job_id, task_id = m.group(1), m.group(2)

        info = {
            "slurm_job_id": job_id,
            "task_id": task_id,
            "pdb_id": None,
            "structure_type": None,
            "datasets": None,
            "date_utc": None,
            "conda_env": None,
            "charmm_path": None,
            "pac_workers": None,
            "results_root": None,
            "topology_files": None,
            "parameter_files": None,
            "n_atoms_retained": None,
            "pac_timeout_budget_s": None,
            "protocol": None,
            "tapbs_timeout": None,
            "ph7_only": None,
        }

        try:
            with open(log_file) as fh:
                text = fh.read()
        except OSError:
            continue

        for line in text.splitlines():
            line = line.strip()
            mm = header_re["pdb"].match(line)
            if mm:
                info["pdb_id"] = mm.group(1).lower()
                info["structure_type"] = mm.group(2)
                info["datasets"] = mm.group(3).strip()
                continue
            mm = header_re["date"].match(line)
            if mm:
                info["date_utc"] = mm.group(1).strip()
                continue
            mm = header_re["conda"].match(line)
            if mm:
                info["conda_env"] = mm.group(1)
                continue
            mm = header_re["charmm"].match(line)
            if mm:
                info["charmm_path"] = mm.group(1).strip()
                continue
            mm = header_re["pac_workers"].match(line)
            if mm:
                info["pac_workers"] = int(mm.group(1))
                continue
            mm = header_re["results_root"].match(line)
            if mm:
                info["results_root"] = mm.group(1).strip()
                continue

        # Multi-line patterns — search whole text
        mm = topology_re.search(text)
        if mm:
            info["topology_files"] = mm.group(1).strip()
        mm = params_re.search(text)
        if mm:
            info["parameter_files"] = mm.group(1).strip()
        mm = atoms_re.search(text)
        if mm:
            info["n_atoms_retained"] = int(mm.group(1))
        mm = timeout_re.search(text)
        if mm:
            info["pac_timeout_budget_s"] = int(mm.group(1))
        mm = protocol_re.search(text)
        if mm:
            info["protocol"] = mm.group(1).strip()
        mm = tapbs_timeout_re.search(text)
        if mm:
            info["tapbs_timeout"] = int(mm.group(1))
        mm = ph7_only_re.search(text)
        if mm:
            info["ph7_only"] = bool(int(mm.group(1)))

        if info["pdb_id"]:
            # Multiple log files might exist per pdb (re-runs); keep most recent
            existing = results.get(info["pdb_id"])
            if existing is None or job_id >= existing.get("slurm_job_id", ""):
                results[info["pdb_id"]] = info

    return results


def parse_convergence(karlsberg_dir: Path) -> dict:
    """
    Parse convergence traces from a *_karlsberg/ directory.

    Returns dict keyed by pac_name -> {
        'iterations', 'final_hamming', 'stop_reason', 'convergence_trace'
    }
    """
    result = {}
    if not karlsberg_dir.is_dir():
        return result

    for pac_dir in sorted(karlsberg_dir.glob("pac_ph_*")):
        trace_file = pac_dir / "convergence_trace.tsv"
        if not trace_file.exists():
            continue
        rows = []
        try:
            with open(trace_file, newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    rows.append(row)
        except OSError:
            continue

        if rows:
            final = rows[-1]
            result[pac_dir.name] = {
                "iterations": len(rows),
                "final_hamming": int(final.get("hamming_distance", -1)),
                "stop_reason": final.get("stop_reason", ""),
                "convergence_trace": rows,
            }

    return result


def parse_log_dat(log_dat_path: Path) -> dict:
    """
    Extract structured information from a KB3 log.dat (ASCII) file.

    Returns a dict with setup metadata.
    """
    info = {}
    if not log_dat_path.exists():
        return info

    try:
        text = log_dat_path.read_text(errors="replace")
    except OSError:
        return info

    # Extract topology files
    m = re.search(r"-topology:\n((?:\s+\S+\n)+)", text)
    if m:
        info["topology"] = [l.strip() for l in m.group(1).strip().splitlines()]

    # Extract parameter files
    m = re.search(r"-parameters:\n((?:\s+\S+\n)+)", text)
    if m:
        info["parameters"] = [l.strip() for l in m.group(1).strip().splitlines()]

    # Extract protocol
    m = re.search(r"-protocol:\n((?:\s+.+\n)+)", text)
    if m:
        info["protocol"] = [l.strip() for l in m.group(1).strip().splitlines()]

    # Extract PAC timeout budget
    m = re.search(r"PAC timeout budget:\s+([\d]+)s\s+\((.+)\)", text)
    if m:
        info["pac_timeout_budget_s"] = int(m.group(1))
        info["pac_timeout_detail"] = m.group(2).strip()

    # Extract salt bridge cutoff
    m = re.search(r"-salt_bridge_cutoff:(\d+)", text)
    if m:
        info["salt_bridge_cutoff"] = int(m.group(1))

    return info


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def compute_stats(deltas: list) -> dict:
    """Compute RMSE, MAE, mean_delta, N from a list of (calc-expt) values."""
    if not deltas:
        return {"n": 0, "rmse": None, "mae": None, "mean_delta": None}
    n = len(deltas)
    mae = sum(abs(d) for d in deltas) / n
    rmse = math.sqrt(sum(d ** 2 for d in deltas) / n)
    mean_delta = sum(deltas) / n
    return {"n": n, "rmse": round(rmse, 3), "mae": round(mae, 3), "mean_delta": round(mean_delta, 3)}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_overview_table(results_root: Path, output_dir: Path, verbose: bool = True) -> None:
    pkad_data = _repo("pkad_data")

    # Load reference data
    exp_pkas = parse_experimental_pkas(pkad_data / "experimental_pkas.tsv")
    super_table = parse_super_table(pkad_data / "super_table.csv")
    job_manifest = parse_job_manifest(pkad_data / "job_manifest.tsv")
    slurm_info = parse_slurm_logs(_repo("logs"))

    if verbose:
        print(f"Loaded {len(exp_pkas)} experimental pKa entries")
        print(f"Loaded {len(super_table)} super_table rows")
        print(f"Loaded {len(job_manifest)} manifest entries")
        print(f"Parsed SLURM logs for {len(slurm_info)} structures")

    # Discover run directories
    run_dirs = sorted(
        [d for d in results_root.iterdir() if d.is_dir()],
        key=lambda d: d.name,
    ) if results_root.is_dir() else []

    if not run_dirs:
        print(f"WARNING: no run directories found in {results_root}", file=sys.stderr)

    # ---------------------------------------------------------------------------
    # Build per-residue rows
    # ---------------------------------------------------------------------------

    overview_rows = []
    run_meta_rows = []

    for run_dir in run_dirs:
        run_name = run_dir.name
        run_date = run_name.split("_")[-1] if "_" in run_name else ""

        # Find all results files
        results_files = sorted(run_dir.glob("*_results.txt"))

        for rf in results_files:
            pdb_id = rf.name.replace("_results.txt", "").lower()
            calc_pkas = parse_results_txt(rf)

            # Karlsberg directory (may not exist if results were cleaned up)
            karlsberg_dir = run_dir / f"{pdb_id}_karlsberg"
            conv = parse_convergence(karlsberg_dir)
            log_setup = parse_log_dat(karlsberg_dir / "log.dat") if karlsberg_dir.is_dir() else {}

            # Convergence summary across PAC runs
            conv_summary = {}
            if conv:
                conv_summary["pac_runs"] = list(conv.keys())
                conv_summary["total_iterations"] = sum(v["iterations"] for v in conv.values())
                conv_summary["all_converged"] = all(
                    v["stop_reason"] == "converged" for v in conv.values()
                )
                max_hamming = max(
                    (v["final_hamming"] for v in conv.values() if v["final_hamming"] >= 0),
                    default=None,
                )
                conv_summary["max_final_hamming"] = max_hamming

            # SLURM info for this structure
            slurm = slurm_info.get(pdb_id, {})
            manifest = job_manifest.get(pdb_id, {})

            # Structure type: prefer SLURM log, fallback to manifest
            structure_type = slurm.get("structure_type") or manifest.get("structure_type", "crystal")
            datasets = slurm.get("datasets") or manifest.get("datasets", "")

            # Build run-level metadata row
            topo = log_setup.get("topology") or slurm.get("topology_files") or ""
            params = log_setup.get("parameters") or slurm.get("parameter_files") or ""
            protocol = log_setup.get("protocol") or slurm.get("protocol") or ""

            run_meta_rows.append({
                "pdb_id": pdb_id,
                "run_name": run_name,
                "run_date": run_date,
                "structure_type": structure_type,
                "datasets": datasets,
                "source_results_file": str(rf.relative_to(REPO_ROOT)),
                "slurm_job_id": slurm.get("slurm_job_id", ""),
                "slurm_task_id": slurm.get("task_id", ""),
                "date_utc": slurm.get("date_utc", ""),
                "conda_env": slurm.get("conda_env", ""),
                "charmm_path": slurm.get("charmm_path", ""),
                "pac_workers": slurm.get("pac_workers", ""),
                "pac_timeout_budget_s": slurm.get("pac_timeout_budget_s", "") or log_setup.get("pac_timeout_budget_s", ""),
                "pac_timeout_detail": log_setup.get("pac_timeout_detail", ""),
                "tapbs_timeout_override": slurm.get("tapbs_timeout", ""),
                "ph7_only": slurm.get("ph7_only", ""),
                "protocol": protocol if isinstance(protocol, str) else "; ".join(protocol),
                "salt_bridge_cutoff": log_setup.get("salt_bridge_cutoff", ""),
                "topology_files": topo if isinstance(topo, str) else "; ".join(topo),
                "parameter_files": params if isinstance(params, str) else "; ".join(params),
                "n_atoms_retained": slurm.get("n_atoms_retained", ""),
                "n_titratable_calc": len(calc_pkas),
                "karlsberg_dir_exists": karlsberg_dir.is_dir(),
                "pac_runs_converged": conv_summary.get("all_converged", ""),
                "pac_total_iterations": conv_summary.get("total_iterations", ""),
                "pac_max_final_hamming": conv_summary.get("max_final_hamming", ""),
                "uniprot_id": manifest.get("uniprot_id", ""),
            })

            # Build per-residue rows: match calculated to experimental
            # For each experimental pKa, look up calculated value
            for (pid, chain, resname_std, resid), expt_pka in exp_pkas.items():
                if pid != pdb_id:
                    continue

                # Try to find the calculated pKa
                # KB3 uses internal residue names; try multiple variants
                calc_pka = None
                key_found = None
                for kb3_name in [resname_std] + [k for k, v in KB3_TO_STANDARD.items() if v == resname_std]:
                    key = (chain, kb3_name, resid)
                    if key in calc_pkas:
                        calc_pka = calc_pkas[key]
                        key_found = key
                        break
                # Also try without chain filter (single-chain proteins may list all in chain A)
                if calc_pka is None:
                    for (c, rn, rid), val in calc_pkas.items():
                        if rid == resid and (norm_resname(rn) == resname_std or rn == resname_std):
                            calc_pka = val
                            key_found = (c, rn, rid)
                            break

                delta = round(calc_pka - expt_pka, 3) if calc_pka is not None else None
                status = "matched" if calc_pka is not None else "missing"

                # Full metadata from super_table
                st_key = (pdb_id, chain, resname_std, resid)
                st_row = super_table.get(st_key, {})

                overview_rows.append({
                    "pdb_id": pdb_id,
                    "chain": chain,
                    "resname": resname_std,
                    "resid": resid,
                    "run_name": run_name,
                    "run_date": run_date,
                    "structure_type": structure_type,
                    "datasets": datasets,
                    # Experimental values
                    "expt_pka": expt_pka,
                    "expt_pka_pkad2": st_row.get("expt_pka_pkad2", ""),
                    "expt_pka_pkad3": st_row.get("expt_pka_pkad3", ""),
                    "expt_pka_pkad_r": st_row.get("expt_pka_pkadR", "") or st_row.get("expt_pka_pkad-r", "") or st_row.get("expt_pka_pkad_r", ""),
                    "pka_conflict": st_row.get("pka_conflict", ""),
                    "expt_uncertainty": st_row.get("expt_uncertainty", ""),
                    "expt_method": st_row.get("expt_method", ""),
                    "reference": st_row.get("reference", ""),
                    "notes": st_row.get("notes", ""),
                    "in_pkad2": st_row.get("in_pkad2", ""),
                    "in_pkad3": st_row.get("in_pkad3", ""),
                    "in_pkad_r": st_row.get("in_pkadR", "") or st_row.get("in_pkad-r", "") or st_row.get("in_pkad_r", ""),
                    "uniprot_id": st_row.get("uniprot_id", "") or manifest.get("uniprot_id", ""),
                    "uniprot_resid": st_row.get("uniprot_resid", ""),
                    # Calculated values
                    "calc_pka": calc_pka if calc_pka is not None else "",
                    "delta_pka": delta if delta is not None else "",
                    "status": status,
                    "calc_resname_kb3": key_found[1] if key_found else "",
                    # Traceback
                    "source_results_file": str(rf.relative_to(REPO_ROOT)),
                    "karlsberg_dir_exists": karlsberg_dir.is_dir(),
                })

            # Also include calculated pKas that have NO experimental match
            matched_keys = set()
            for row in overview_rows:
                if row["pdb_id"] == pdb_id and row["run_name"] == run_name and row["status"] == "matched":
                    matched_keys.add((row["chain"], row["resname"], row["resid"]))

            for (chain, kb3_name, resid), calc_pka in calc_pkas.items():
                std_name = norm_resname(kb3_name)
                if (chain, std_name, resid) not in matched_keys:
                    # Check if there is an experimental entry for this residue
                    exp_key = (pdb_id, chain, std_name, resid)
                    expt_pka = exp_pkas.get(exp_key)
                    if expt_pka is None:
                        status = "no_expt"
                    else:
                        status = "matched"

                    st_key = (pdb_id, chain, std_name, resid)
                    st_row = super_table.get(st_key, {})

                    overview_rows.append({
                        "pdb_id": pdb_id,
                        "chain": chain,
                        "resname": std_name,
                        "resid": resid,
                        "run_name": run_name,
                        "run_date": run_date,
                        "structure_type": structure_type,
                        "datasets": datasets,
                        "expt_pka": expt_pka if expt_pka is not None else "",
                        "expt_pka_pkad2": st_row.get("expt_pka_pkad2", ""),
                        "expt_pka_pkad3": st_row.get("expt_pka_pkad3", ""),
                        "expt_pka_pkad_r": "",
                        "pka_conflict": st_row.get("pka_conflict", ""),
                        "expt_uncertainty": st_row.get("expt_uncertainty", ""),
                        "expt_method": st_row.get("expt_method", ""),
                        "reference": st_row.get("reference", ""),
                        "notes": st_row.get("notes", ""),
                        "in_pkad2": st_row.get("in_pkad2", ""),
                        "in_pkad3": st_row.get("in_pkad3", ""),
                        "in_pkad_r": "",
                        "uniprot_id": st_row.get("uniprot_id", "") or manifest.get("uniprot_id", ""),
                        "uniprot_resid": st_row.get("uniprot_resid", ""),
                        "calc_pka": calc_pka,
                        "delta_pka": round(calc_pka - expt_pka, 3) if expt_pka is not None else "",
                        "status": status,
                        "calc_resname_kb3": kb3_name,
                        "source_results_file": str(rf.relative_to(REPO_ROOT)),
                        "karlsberg_dir_exists": karlsberg_dir.is_dir(),
                    })

    # Also mark "not_run" structures (in manifest but no results file found)
    run_pdb_ids = set(row["pdb_id"] for row in run_meta_rows)
    for pdb_id, mrow in job_manifest.items():
        if pdb_id not in run_pdb_ids:
            # No results at all
            for (pid, chain, resname, resid), expt_pka in exp_pkas.items():
                if pid != pdb_id:
                    continue
                st_key = (pdb_id, chain, resname, resid)
                st_row = super_table.get(st_key, {})
                overview_rows.append({
                    "pdb_id": pdb_id,
                    "chain": chain,
                    "resname": resname,
                    "resid": resid,
                    "run_name": "not_run",
                    "run_date": "",
                    "structure_type": mrow.get("structure_type", "crystal"),
                    "datasets": mrow.get("datasets", ""),
                    "expt_pka": expt_pka,
                    "expt_pka_pkad2": st_row.get("expt_pka_pkad2", ""),
                    "expt_pka_pkad3": st_row.get("expt_pka_pkad3", ""),
                    "expt_pka_pkad_r": "",
                    "pka_conflict": st_row.get("pka_conflict", ""),
                    "expt_uncertainty": st_row.get("expt_uncertainty", ""),
                    "expt_method": st_row.get("expt_method", ""),
                    "reference": st_row.get("reference", ""),
                    "notes": st_row.get("notes", ""),
                    "in_pkad2": st_row.get("in_pkad2", ""),
                    "in_pkad3": st_row.get("in_pkad3", ""),
                    "in_pkad_r": "",
                    "uniprot_id": st_row.get("uniprot_id", "") or mrow.get("uniprot_id", ""),
                    "uniprot_resid": st_row.get("uniprot_resid", ""),
                    "calc_pka": "",
                    "delta_pka": "",
                    "status": "not_run",
                    "calc_resname_kb3": "",
                    "source_results_file": "",
                    "karlsberg_dir_exists": False,
                })

    # ---------------------------------------------------------------------------
    # Write outputs
    # ---------------------------------------------------------------------------
    output_dir.mkdir(parents=True, exist_ok=True)

    # Overview table
    overview_path = output_dir / "pkad_overview_table.csv"
    if overview_rows:
        fieldnames = list(overview_rows[0].keys())
        with open(overview_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(overview_rows)
        print(f"Written: {overview_path}  ({len(overview_rows)} rows)")

    # Run metadata
    meta_path = output_dir / "pkad_run_metadata.csv"
    if run_meta_rows:
        fieldnames = list(run_meta_rows[0].keys())
        with open(meta_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(run_meta_rows)
        print(f"Written: {meta_path}  ({len(run_meta_rows)} rows)")

    # Statistics
    stats_path = output_dir / "pkad_statistics.txt"
    _write_statistics(overview_rows, stats_path)
    print(f"Written: {stats_path}")


def _write_statistics(rows: list, output_path: Path) -> None:
    """Compute and write per-dataset and overall RMSE/MAE/N statistics."""
    # Group deltas by (run_name, dataset, resname)
    all_deltas = []
    by_run = defaultdict(list)
    by_dataset = defaultdict(list)
    by_resname = defaultdict(list)
    by_run_dataset = defaultdict(list)

    for row in rows:
        if row["status"] != "matched" or row["delta_pka"] == "":
            continue
        d = float(row["delta_pka"])
        rn = row["run_name"]
        ds_str = row.get("datasets", "") or ""
        resname = row["resname"]

        all_deltas.append(d)
        by_run[rn].append(d)
        by_resname[resname].append(d)
        for ds in [x.strip() for x in ds_str.split(",") if x.strip()]:
            by_dataset[ds].append(d)
            by_run_dataset[(rn, ds)].append(d)

    lines = ["=" * 60, "PKAD Benchmark Statistics", "=" * 60, ""]

    lines.append("--- Overall ---")
    s = compute_stats(all_deltas)
    lines.append(f"  N={s['n']}  RMSE={s['rmse']}  MAE={s['mae']}  MeanDelta={s['mean_delta']}")
    lines.append("")

    lines.append("--- By run ---")
    for run, deltas in sorted(by_run.items()):
        s = compute_stats(deltas)
        lines.append(f"  {run}: N={s['n']}  RMSE={s['rmse']}  MAE={s['mae']}  MeanDelta={s['mean_delta']}")
    lines.append("")

    lines.append("--- By dataset ---")
    for ds, deltas in sorted(by_dataset.items()):
        s = compute_stats(deltas)
        lines.append(f"  {ds}: N={s['n']}  RMSE={s['rmse']}  MAE={s['mae']}  MeanDelta={s['mean_delta']}")
    lines.append("")

    lines.append("--- By residue type ---")
    for resname, deltas in sorted(by_resname.items()):
        s = compute_stats(deltas)
        lines.append(f"  {resname}: N={s['n']}  RMSE={s['rmse']}  MAE={s['mae']}  MeanDelta={s['mean_delta']}")
    lines.append("")

    lines.append("--- By run × dataset ---")
    for (run, ds), deltas in sorted(by_run_dataset.items()):
        s = compute_stats(deltas)
        lines.append(f"  {run} | {ds}: N={s['n']}  RMSE={s['rmse']}  MAE={s['mae']}  MeanDelta={s['mean_delta']}")

    output_path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--results-root",
        type=Path,
        default=_repo("results"),
        help="Root directory containing run subdirectories (default: repo_root/results/)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent,
        help="Directory to write output files (default: same dir as this script)",
    )
    parser.add_argument("--verbose", "-v", action="store_true", default=True)
    args = parser.parse_args()

    build_overview_table(
        results_root=args.results_root,
        output_dir=args.output_dir,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
