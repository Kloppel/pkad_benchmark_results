"""
Build protonated protein structures from KB3 pKa predictions using CHARMM.

Given a heavy-atom PDB and a KB3 *_results.txt pKa file, this script:
  1. Parses predicted pKa values
  2. Assigns protonation states at the requested pH (default 7.0)
  3. Builds a protonated PDB (with H atoms) via CHARMM initial modelling
  4. Writes the protonated structure to an output PDB file

The CHARMM topology/parameter files used are the same ones used during
the benchmark calculations and are bundled in the protonation/toppar/ directory.
CHARMM itself must be installed separately and available either on PATH or via
the --charmm flag.

Dependencies:
  - CHARMM (c46 or later, or charmm-lite c46)
  - biopython (optional, for PDB parsing; falls back to manual parser)

Usage examples:
    # Single structure at pH 7.0
    python protonate_from_pkas.py \\
        --pdb structures/1ao6.pdb \\
        --pkas results/pkad_benchmark_20260619/1ao6_results.txt \\
        --output protonated/1ao6_pH7.pdb

    # Batch: all structures in the benchmark
    python protonate_from_pkas.py \\
        --batch \\
        --pdb-dir structures/ \\
        --results-dir results/pkad_benchmark_20260619/ \\
        --output-dir protonated/ \\
        --ph 7.0

    # Specify pH and CHARMM binary
    python protonate_from_pkas.py \\
        --pdb structures/1bni.pdb \\
        --pkas results/pkad_benchmark_20260619/1bni_results.txt \\
        --ph 5.0 \\
        --charmm /projects/biomodeling/charmm-lite/c46b1/ser/bin/charmm \\
        --output protonated/1bni_pH5.pdb
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).parent.parent
TOPPAR_DIR = Path(__file__).parent / "toppar"

# CHARMM topology / parameter file loading order
# Must match the order used in the benchmark calculations
_TOPOLOGY_FILES = [
    TOPPAR_DIR / "top.inp",
    TOPPAR_DIR / "converted_0.rtf",
    TOPPAR_DIR / "converted_1.rtf",
]
_PARAMETER_FILES = [
    TOPPAR_DIR / "par.inp",
    TOPPAR_DIR / "converted_0.prm",
    TOPPAR_DIR / "converted_1.prm",
    TOPPAR_DIR / "toppar_water_ions.str",
    TOPPAR_DIR / "par_all36m_prot_additions.prm",
    TOPPAR_DIR / "backfill.prm",
]

# ---------------------------------------------------------------------------
# Residue name mapping
# ---------------------------------------------------------------------------

# KB3 internal names (as appear in *_results.txt) → CHARMM/PDB standard names
KB3_TO_CHARMM = {
    "EPP": "GLU",
    "DPP": "ASP",
    "HSP": "HSP",   # protonated HIS (positive)
    "HSE": "HSE",   # epsilon-protonated HIS (neutral)
    "HSD": "HSD",   # delta-protonated HIS (neutral)
    "NTE": "NTE",
    "CTE": "CTE",
}

# Residue types with standard pKa reference values (approximate)
STD_PKA = {
    "ASP": 3.8,
    "DPP": 3.8,
    "GLU": 4.2,
    "EPP": 4.2,
    "HIS": 6.5,
    "HSP": 6.5,
    "HSE": 6.5,
    "HSD": 6.5,
    "CYS": 8.3,
    "TYR": 10.5,
    "LYS": 10.5,
    "ARG": 12.5,
    "NTE": 8.0,
    "CTE": 3.3,
}

# For each residue type, what CHARMM patch to apply for the protonated form
# (at pH << pKa → protonated; at pH >> pKa → deprotonated)
#
# Convention used by KB3/CHARMM:
#   - ASP/EPP at low pH: ASPP patch (adds H to one carboxylate O)
#   - GLU/DPP at low pH: GLUP patch
#   - HIS: HSD (delta-H only) / HSE (epsilon-H only) / HSP (both H, +1)
#   - CYS: default is thiol (SH); deprotonated CYS has no special patch
#   - TYR: default is phenol (OH); deprotonated TYR → TYRO patch
#   - LYS: default is ammonium (NH3+); deprotonated → LYSO (rare)
#   - N-terminus: NTER (default) or PROP (proline), or specific charged forms
#   - C-terminus: CTER (default)

PROTONATED_PATCH = {
    "ASP": "ASPP",
    "DPP": "ASPP",
    "GLU": "GLUP",
    "EPP": "GLUP",
    "HIS": "HSP",   # use HSP patch (both H, charged)
}
DEPROTONATED_PATCH = {
    "TYR": "TYRO",  # deprotonated tyrosine (TYR-O-)
    "CYS": "CYSD",  # deprotonated cysteine (CYS-S-)
    "LYS": "LYSO",  # deprotonated lysine (neutral amine; very rare at pH 7)
}


# ---------------------------------------------------------------------------
# pKa file parser
# ---------------------------------------------------------------------------

def parse_pka_results(path: Path) -> dict:
    """
    Parse KB3 *_results.txt.

    Returns dict: (chain, resname, resid) → pKa (float).
    Line format: '       LYS-1_A: 11.72 '
    """
    pkas = {}
    pattern = re.compile(r"^\s*(\w+)-(\d+)_([A-Z]):\s*([\d\.\-]+)\s*$")
    with open(path) as fh:
        for line in fh:
            m = pattern.match(line)
            if m:
                resname, resid, chain, value = m.groups()
                try:
                    pkas[(chain, resname.upper(), int(resid))] = float(value)
                except ValueError:
                    pass
    return pkas


# ---------------------------------------------------------------------------
# Protonation state determination
# ---------------------------------------------------------------------------

def determine_protonation_states(pkas: dict, ph: float) -> dict:
    """
    For each residue, decide which CHARMM residue name / patch to use.

    Args:
        pkas: dict from parse_pka_results
        ph:   target pH

    Returns:
        dict: (chain, resid) → {
            'resname_kb3': str,
            'resname_charmm': str,
            'pka': float,
            'is_protonated': bool,
            'patch': str or None
        }
    """
    states = {}

    for (chain, kb3_name, resid), pka in pkas.items():
        is_protonated = ph < pka  # True → more H than typical deprotonated state

        # Determine CHARMM residue name
        resname_charmm = KB3_TO_CHARMM.get(kb3_name, kb3_name)
        patch = None

        if kb3_name in ("DPP", "EPP", "ASP", "GLU"):
            # Acid residues: protonated = uncharged form (pH < pKa → ASPP/GLUP)
            # At pH 7, most ASP/GLU are deprotonated (charged); only patch if pKa > pH
            resname_charmm = "ASP" if kb3_name in ("DPP", "ASP") else "GLU"
            if is_protonated:
                patch = "ASPP" if kb3_name in ("DPP", "ASP") else "GLUP"

        elif kb3_name in ("HSP", "HSE", "HSD", "HIS"):
            # Histidine: protonated = HSP (doubly protonated, +1)
            # KB3 naming: HSP = protonated, HSE = epsilon-H only, HSD = delta-H only
            if kb3_name == "HSP":
                # Already the protonated form
                resname_charmm = "HSP"
            elif is_protonated:
                resname_charmm = "HSP"
            else:
                # Default to HSE (epsilon tautomer) for neutral HIS
                resname_charmm = "HSE"

        elif kb3_name == "TYR":
            resname_charmm = "TYR"
            if not is_protonated:
                patch = "TYRO"

        elif kb3_name == "CYS":
            resname_charmm = "CYS"
            if not is_protonated:
                patch = "CYSD"

        elif kb3_name == "LYS":
            resname_charmm = "LYS"
            if not is_protonated:
                patch = "LYSO"

        elif kb3_name == "ARG":
            resname_charmm = "ARG"
            # ARG is almost never deprotonated at normal pH; pKa ~ 12-13

        elif kb3_name == "NTE":
            # N-terminal: protonated = NTER (default), deprotonated = NTRD (rare)
            resname_charmm = "NTER"

        elif kb3_name == "CTE":
            # C-terminal: protonated = CTER (default), deprotonated = CTRD (rare)
            resname_charmm = "CTER"

        states[(chain, resid)] = {
            "resname_kb3": kb3_name,
            "resname_charmm": resname_charmm,
            "pka": pka,
            "ph": ph,
            "is_protonated": is_protonated,
            "patch": patch,
        }

    return states


# ---------------------------------------------------------------------------
# PDB parser (lightweight, no biopython needed)
# ---------------------------------------------------------------------------

def parse_pdb_segments(pdb_path: Path) -> dict:
    """
    Parse a PDB file and return a dict: chain → list of (resname, resid, atom lines).
    Only considers ATOM/HETATM records for standard amino acids.
    """
    segments = {}
    with open(pdb_path) as fh:
        for line in fh:
            if line[:6] not in ("ATOM  ", "HETATM"):
                continue
            chain = line[21]
            resid = int(line[22:26].strip())
            resname = line[17:20].strip()
            if chain not in segments:
                segments[chain] = []
            segments[chain].append((resname, resid, line))
    return segments


def get_unique_residues(pdb_path: Path) -> dict:
    """Return dict chain → list of (resname, resid) tuples in order."""
    segs = parse_pdb_segments(pdb_path)
    result = {}
    for chain, records in segs.items():
        seen = set()
        residues = []
        for resname, resid, _ in records:
            if (resid, resname) not in seen:
                seen.add((resid, resname))
                residues.append((resname, resid))
        result[chain] = residues
    return result


def detect_disulfide_bonds(pdb_path: Path, cutoff_angstrom: float = 2.5) -> list:
    """
    Detect disulfide bonds by proximity of CYS SG atoms.

    Returns list of ((chainA, resA), (chainB, resB)) tuples.
    """
    sg_atoms = {}
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if line[:6] not in ("ATOM  ", "HETATM"):
                    continue
                resname = line[17:20].strip().upper()
                atom_name = line[12:16].strip().upper()
                if resname == "CYS" and atom_name == "SG":
                    chain = line[21]
                    resid = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    sg_atoms[(chain, resid)] = (x, y, z)
    except (OSError, ValueError):
        return []

    bonds = []
    keys = list(sg_atoms.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            k1, k2 = keys[i], keys[j]
            x1, y1, z1 = sg_atoms[k1]
            x2, y2, z2 = sg_atoms[k2]
            dist = ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5
            if dist <= cutoff_angstrom:
                bonds.append((k1, k2))
    return bonds


# ---------------------------------------------------------------------------
# CHARMM script generation
# ---------------------------------------------------------------------------

def build_charmm_crd(pdb_path: Path, workdir: Path, protonation_states: dict) -> dict:
    """
    Convert PDB to per-chain CHARMM CRD coordinate files.

    Writes {chain}_in.crd files in workdir.
    Returns dict chain → crd_path.
    """
    # Try biopython first
    try:
        from Bio.PDB import PDBParser, PDBIO
        _build_crd_biopython(pdb_path, workdir, protonation_states)
    except ImportError:
        pass

    # Manual fallback: write PDB-format coordinate sections
    crd_paths = {}
    residues = get_unique_residues(pdb_path)
    segments = parse_pdb_segments(pdb_path)

    for chain, res_list in residues.items():
        crd_path = workdir / f"{chain}_in.crd"
        states = {(chain, resid): st for (c, resid), st in protonation_states.items() if c == chain}

        lines = []
        atom_serial = 0
        resid_seq = 0
        prev_resid = None

        for resname_orig, resid in res_list:
            charmm_resname = resname_orig
            if (chain, resid) in states:
                charmm_resname = states[(chain, resid)]["resname_charmm"]
            if resid != prev_resid:
                resid_seq += 1
                prev_resid = resid

        # Write a simplified CRD (CHARMM CARD format)
        # CHARMM CRD: NATOM header then lines of:
        # serial  resid  resname  atname  x  y  z  segid  resid_orig  weighting
        atom_records = []
        for resname_orig, resid in res_list:
            charmm_resname = resname_orig
            if (chain, resid) in states:
                charmm_resname = states[(chain, resid)]["resname_charmm"]

            for rn, rid, atom_line in segments[chain]:
                if rid != resid or rn != resname_orig:
                    continue
                atom_name = atom_line[12:16].strip()
                try:
                    x = float(atom_line[30:38])
                    y = float(atom_line[38:46])
                    z = float(atom_line[46:54])
                except (ValueError, IndexError):
                    continue
                atom_serial += 1
                # CRD format: 5i,1x,5i,1x,4a4,2f10.5,2x,4a4,1x,5i,2f10.5
                # Simplified plain format (CHARMM accepts this):
                atom_records.append(
                    f"{atom_serial:5d}{resid_seq:5d} {charmm_resname:<4s} {atom_name:<4s}"
                    f"{x:10.5f}{y:10.5f}{z:10.5f}  {chain:<4s}{resid:5d}{0.0:10.5f}\n"
                )

        crd_content = f"* Coordinates for chain {chain} from {pdb_path.name}\n"
        crd_content += f"*\n"
        crd_content += f"{atom_serial:5d}\n"
        crd_content += "".join(atom_records)
        crd_path.write_text(crd_content)
        crd_paths[chain] = crd_path

    return crd_paths


def generate_charmm_inp(
    pdb_path: Path,
    workdir: Path,
    output_pdb: Path,
    protonation_states: dict,
    disulfide_bonds: list,
    chains: list,
) -> Path:
    """
    Generate a CHARMM .inp script for building a protonated structure.

    Returns the path to the generated .inp file.
    """
    pdb_stem = pdb_path.stem
    inp_path = workdir / f"{pdb_stem}_build.inp"

    def topo_block() -> str:
        lines = []
        for i, tf in enumerate(_TOPOLOGY_FILES):
            if tf.exists():
                if i == 0:
                    lines.append(f'read rtf card name "{tf}"')
                else:
                    lines.append(f'read rtf card name "{tf}" append')
        return "\n".join(lines)

    def param_block() -> str:
        lines = []
        for i, pf in enumerate(_PARAMETER_FILES):
            if pf.exists():
                if i == 0:
                    lines.append(f'read para card flex name "{pf}"')
                else:
                    lines.append(f'read para card flex name "{pf}" append')
        return "\n".join(lines)

    # Build per-chain generate + patch blocks
    seg_blocks = []
    patch_blocks = []

    for chain in chains:
        crd_file = workdir / f"{chain}_in.crd"
        seg_block = f"""
read sequence coor name "{crd_file}" resid
generate setup {chain} first ACE last CT1

read coor card name "{crd_file}" resid

hbuild sele segid {chain} .and. hydrogen end
ic para
ic build
hbuild sele segid {chain} .and. hydrogen end
ic para
ic build

coor print sele segid {chain} .and. .not. INIT end
"""
        seg_blocks.append(seg_block)

        # Patches for this chain
        for (ch, resid), state in protonation_states.items():
            if ch != chain:
                continue
            patch = state.get("patch")
            if patch:
                patch_blocks.append(f"patch {patch} {chain} {resid}")

    # Disulfide bond patches
    disu_lines = []
    for (c1, r1), (c2, r2) in disulfide_bonds:
        disu_lines.append(f"patch disu {c1} {r1} {c2} {r2}")

    # Write coordinate after patching
    write_block = f"""
! Write protonated structure
open write unit 1 card name "{output_pdb}"
write coor pdb unit 1
close unit 1
"""

    inp_content = f"""! CHARMM build script generated by protonate_from_pkas.py
! PDB: {pdb_path}
! pH at which protonation states were assigned: see patches below

dimension chsize 999999
bomlev -2

{topo_block()}

{param_block()}

{"".join(seg_blocks)}

! Protonation state patches (from KB3 pKa prediction)
{chr(10).join(patch_blocks)}

! Disulfide bond patches (detected by SG-SG distance < 2.5 Å)
{chr(10).join(disu_lines)}

! Rebuild after patching
ic para
ic build

{write_block}

stop
"""

    inp_path.write_text(inp_content)
    return inp_path


# ---------------------------------------------------------------------------
# CHARMM execution
# ---------------------------------------------------------------------------

def find_charmm(hint: str = None) -> str:
    """Locate CHARMM binary. Returns the command string or raises RuntimeError."""
    candidates = []
    if hint:
        candidates.append(hint)
    # Common cluster locations
    candidates += [
        "charmm",
        "/projects/biomodeling/charmm-lite/c46b1/ser/bin/charmm",
        "/usr/local/bin/charmm",
    ]
    # Also check CHARMM env var
    env_charmm = os.environ.get("CHARMM")
    if env_charmm:
        candidates.insert(0, env_charmm)

    for c in candidates:
        if shutil.which(c) or Path(c).is_file():
            return c
    raise RuntimeError(
        "CHARMM binary not found. Use --charmm /path/to/charmm or set CHARMM env var."
    )


def run_charmm(charmm_cmd: str, inp_path: Path, workdir: Path) -> tuple:
    """
    Run CHARMM with the given input script.

    Returns (returncode, stdout, stderr).
    """
    out_path = inp_path.with_suffix(".out")
    err_path = inp_path.with_suffix(".err")

    with open(inp_path) as inp_fh, open(out_path, "w") as out_fh, open(err_path, "w") as err_fh:
        result = subprocess.run(
            [charmm_cmd],
            stdin=inp_fh,
            stdout=out_fh,
            stderr=err_fh,
            cwd=workdir,
            timeout=600,  # 10 minutes
        )
    stdout = out_path.read_text(errors="replace") if out_path.exists() else ""
    stderr = err_path.read_text(errors="replace") if err_path.exists() else ""
    return result.returncode, stdout, stderr


# ---------------------------------------------------------------------------
# Single-structure protonation
# ---------------------------------------------------------------------------

def protonate_structure(
    pdb_path: Path,
    pka_path: Path,
    output_pdb: Path,
    ph: float = 7.0,
    charmm_cmd: str = None,
    workdir: Path = None,
    keep_workdir: bool = False,
    verbose: bool = True,
) -> bool:
    """
    Protonate a single structure.

    Returns True on success, False on failure.
    """
    if verbose:
        print(f"\n{'=' * 60}")
        print(f"Structure : {pdb_path.name}")
        print(f"pKa file  : {pka_path.name}")
        print(f"pH        : {ph}")
        print(f"Output    : {output_pdb}")

    # Find CHARMM
    try:
        charmm = find_charmm(charmm_cmd)
        if verbose:
            print(f"CHARMM    : {charmm}")
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return False

    # Parse pKa values
    pkas = parse_pka_results(pka_path)
    if verbose:
        print(f"Parsed pKa values: {len(pkas)} titratable residues")

    # Determine protonation states
    states = determine_protonation_states(pkas, ph)
    n_protonated = sum(1 for s in states.values() if s["is_protonated"])
    n_patches = sum(1 for s in states.values() if s["patch"])
    if verbose:
        print(f"Protonation states: {n_protonated}/{len(states)} protonated at pH {ph}")
        print(f"Non-default patches: {n_patches}")

    # Detect chains and disulfide bonds
    residues_by_chain = get_unique_residues(pdb_path)
    chains = sorted(residues_by_chain.keys())
    disulfide_bonds = detect_disulfide_bonds(pdb_path)
    if verbose and disulfide_bonds:
        print(f"Disulfide bonds detected: {len(disulfide_bonds)}")

    # Set up working directory
    cleanup = False
    if workdir is None:
        workdir = Path(tempfile.mkdtemp(prefix="kb3_prot_"))
        cleanup = not keep_workdir

    workdir.mkdir(parents=True, exist_ok=True)

    try:
        # Build CRD files
        crd_paths = build_charmm_crd(pdb_path, workdir, states)

        # Generate CHARMM input
        inp_path = generate_charmm_inp(
            pdb_path=pdb_path,
            workdir=workdir,
            output_pdb=output_pdb.resolve(),
            protonation_states=states,
            disulfide_bonds=disulfide_bonds,
            chains=chains,
        )

        if verbose:
            print(f"CHARMM input: {inp_path}")

        # Run CHARMM
        rc, stdout, stderr = run_charmm(charmm, inp_path, workdir)

        if rc != 0:
            print(f"ERROR: CHARMM exited with code {rc}", file=sys.stderr)
            # Print last 20 lines of CHARMM output for diagnosis
            print("--- Last 20 lines of CHARMM output ---", file=sys.stderr)
            for line in stdout.splitlines()[-20:]:
                print(f"  {line}", file=sys.stderr)
            return False

        if output_pdb.exists():
            if verbose:
                print(f"SUCCESS: protonated structure written to {output_pdb}")
            return True
        else:
            print(f"ERROR: CHARMM ran without error but output PDB was not written.", file=sys.stderr)
            return False

    finally:
        if cleanup and workdir.exists():
            shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Print protonation state report (without running CHARMM)
# ---------------------------------------------------------------------------

def print_protonation_report(pka_path: Path, ph: float, output_path: Path = None) -> None:
    """Print a human-readable protonation state assignment table."""
    pkas = parse_pka_results(pka_path)
    states = determine_protonation_states(pkas, ph)

    lines = [
        f"Protonation states at pH {ph}",
        f"{'Chain':>5} {'ResID':>6} {'ResName':>8} {'KB3Name':>8} {'pKa':>7} {'Protonated':>11} {'Patch':>8}",
        "-" * 60,
    ]

    for (chain, resid), st in sorted(states.items()):
        prot = "YES" if st["is_protonated"] else "no"
        patch = st["patch"] or ""
        lines.append(
            f"{chain:>5} {resid:>6} {st['resname_charmm']:>8} {st['resname_kb3']:>8} "
            f"{st['pka']:>7.2f} {prot:>11} {patch:>8}"
        )

    report = "\n".join(lines) + "\n"

    if output_path:
        output_path.write_text(report)
    else:
        print(report)


# ---------------------------------------------------------------------------
# Batch mode
# ---------------------------------------------------------------------------

def batch_protonate(
    pdb_dir: Path,
    results_dir: Path,
    output_dir: Path,
    ph: float = 7.0,
    charmm_cmd: str = None,
    verbose: bool = True,
) -> None:
    """Protonate all structures found in pdb_dir that have matching *_results.txt files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    n_success = 0
    n_fail = 0
    n_skip = 0

    for pdb_file in pdb_files:
        pdb_id = pdb_file.stem.lower()
        pka_file = results_dir / f"{pdb_id}_results.txt"
        out_pdb = output_dir / f"{pdb_id}_pH{ph:.1f}.pdb"

        if out_pdb.exists():
            if verbose:
                print(f"Skipping {pdb_id} (output already exists)")
            n_skip += 1
            continue

        if not pka_file.exists():
            if verbose:
                print(f"No pKa file for {pdb_id}, skipping")
            n_skip += 1
            continue

        ok = protonate_structure(
            pdb_path=pdb_file,
            pka_path=pka_file,
            output_pdb=out_pdb,
            ph=ph,
            charmm_cmd=charmm_cmd,
            verbose=verbose,
        )
        if ok:
            n_success += 1
        else:
            n_fail += 1

    print(f"\nBatch complete: {n_success} success, {n_fail} failed, {n_skip} skipped")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--batch", action="store_true", help="Batch-protonate all structures")
    mode.add_argument(
        "--report-only",
        action="store_true",
        help="Print protonation state table without running CHARMM",
    )

    # Single-structure options
    parser.add_argument("--pdb", type=Path, help="Input heavy-atom PDB file")
    parser.add_argument("--pkas", type=Path, help="KB3 *_results.txt pKa file")
    parser.add_argument("--output", type=Path, help="Output protonated PDB")

    # Batch options
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=REPO_ROOT / "structures",
        help="Directory of input PDB files (batch mode)",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=REPO_ROOT / "results" / "pkad_benchmark_20260623_combined",
        help="Directory with *_results.txt files (batch mode)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent / "protonated",
        help="Output directory (batch mode)",
    )

    # Shared options
    parser.add_argument("--ph", type=float, default=7.0, help="Target pH (default: 7.0)")
    parser.add_argument("--charmm", type=str, help="Path to CHARMM binary")
    parser.add_argument("--workdir", type=Path, help="Working directory for CHARMM temp files")
    parser.add_argument("--keep-workdir", action="store_true", help="Do not delete working dir after run")
    parser.add_argument("--quiet", "-q", action="store_true")

    args = parser.parse_args()
    verbose = not args.quiet

    if args.batch:
        batch_protonate(
            pdb_dir=args.pdb_dir,
            results_dir=args.results_dir,
            output_dir=args.output_dir,
            ph=args.ph,
            charmm_cmd=args.charmm,
            verbose=verbose,
        )

    elif args.report_only:
        if not args.pkas:
            parser.error("--report-only requires --pkas")
        print_protonation_report(args.pkas, args.ph)

    else:
        # Single structure mode
        if not args.pdb or not args.pkas:
            parser.error("Single-structure mode requires --pdb and --pkas")
        out_pdb = args.output or args.pdb.with_stem(args.pdb.stem + f"_pH{args.ph:.1f}")
        protonate_structure(
            pdb_path=args.pdb,
            pka_path=args.pkas,
            output_pdb=out_pdb,
            ph=args.ph,
            charmm_cmd=args.charmm,
            workdir=args.workdir,
            keep_workdir=args.keep_workdir,
            verbose=verbose,
        )


if __name__ == "__main__":
    main()
