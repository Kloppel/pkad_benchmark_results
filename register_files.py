#!/usr/bin/env python3
"""register_files.py — Register all PDB files available for the benchmark.

Generates two registers:
  pkad_data/pdb_register.csv   — crystal heavy-atom + protonated files
  pkad_data/af3_register.csv   — AlphaFold DB prediction + protonated files

For each structure, the register records:
  - Whether the raw heavy-atom PDB exists and its atom count
  - Whether the KB3-generated protonated PDB exists and its atom count
    (written by CHARMM during initial modelling)

Protonated PDB path pattern (from KB3 result tree):
  {results_root}/{run_name}/{run_name}_karlsberg/
      initial_modelling/{run_name}_stripped_out.pdb

Usage
-----
    python3 register_files.py [--data-dir pkad_data]
                              [--results-root results/pkad_benchmark_YYYYMMDD]
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import pandas as pd

HERE = Path(__file__).parent
DEFAULT_DATA_DIR = HERE / "pkad_data"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def count_atoms(pdb_path: Path) -> int:
    """Count ATOM record lines in a PDB file."""
    if not pdb_path.exists():
        return 0
    n = 0
    for line in pdb_path.read_text(errors="replace").splitlines():
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            n += 1
    return n


def find_protonated_pdb(results_root: Path, run_name: str) -> Optional[Path]:
    """Locate the CHARMM-output protonated PDB for a given run_name.

    KB3 writes the protonated structure during initial_modelling.
    Checks two common filename patterns.
    """
    if results_root is None:
        return None
    base = results_root / run_name / f"{run_name}_karlsberg" / "initial_modelling"
    for stem in (f"{run_name}_stripped_out", f"{run_name}_out"):
        p = base / f"{stem}.pdb"
        if p.exists():
            return p
    return None


def _row_dict(
    pdb_id: str,
    heavy_path: Path,
    prot_path: Optional[Path],
) -> dict:
    heavy_exists = heavy_path.exists()
    prot_exists = prot_path is not None and prot_path.exists()
    return {
        "pdb_id":               pdb_id,
        "heavy_atom_path":      str(heavy_path),
        "heavy_atom_exists":    heavy_exists,
        "heavy_atom_n_atoms":   count_atoms(heavy_path) if heavy_exists else 0,
        "protonated_path":      str(prot_path) if prot_path else "",
        "protonated_exists":    prot_exists,
        "protonated_n_atoms":   count_atoms(prot_path) if prot_exists else 0,
    }


# ---------------------------------------------------------------------------
# Crystal register
# ---------------------------------------------------------------------------

def build_crystal_register(
    crystal_pdb_dir: Path,
    results_root: Optional[Path],
    pdb_ids: list[str],
) -> pd.DataFrame:
    """One row per pdb_id: heavy-atom + protonated file info."""
    rows = []
    for pdb_id in sorted(pdb_ids):
        heavy = crystal_pdb_dir / f"{pdb_id}.pdb"
        prot = find_protonated_pdb(results_root, pdb_id) if results_root else None
        rows.append(_row_dict(pdb_id, heavy, prot))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# AF3 register
# ---------------------------------------------------------------------------

def build_af3_register(
    af3_pdb_dir: Path,
    results_root: Optional[Path],
    uniprot_pdb_pairs: list[tuple[str, str]],
) -> pd.DataFrame:
    """One row per (pdb_id, uniprot_id): AF3 prediction + protonated file info."""
    rows = []
    seen: set[str] = set()
    for uniprot_id, pdb_id in sorted(uniprot_pdb_pairs, key=lambda x: x[1]):
        if not uniprot_id or uniprot_id in seen:
            continue
        seen.add(uniprot_id)
        run_name = f"{pdb_id}_af3"
        heavy = af3_pdb_dir / f"AF-{uniprot_id}-F1-model_v4.pdb"
        prot = find_protonated_pdb(results_root, run_name) if results_root else None

        heavy_exists = heavy.exists()
        prot_exists = prot is not None and prot.exists()
        rows.append({
            "pdb_id":               pdb_id,
            "uniprot_id":           uniprot_id,
            "af3_path":             str(heavy),
            "af3_exists":           heavy_exists,
            "af3_n_atoms":          count_atoms(heavy) if heavy_exists else 0,
            "protonated_path":      str(prot) if prot else "",
            "protonated_exists":    prot_exists,
            "protonated_n_atoms":   count_atoms(prot) if prot_exists else 0,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR,
                        help="Data directory with super_table.csv and PDB subdirs.")
    parser.add_argument("--results-root", type=Path, default=None,
                        help="Root of completed KB3 results for protonated PDB lookup.")
    args = parser.parse_args()

    data_dir: Path = args.data_dir
    crystal_dir = data_dir / "pdbs" / "heavy_atom"
    af3_dir = data_dir / "af3_pdbs"
    results_root = args.results_root

    super_table_path = data_dir / "super_table.csv"
    if not super_table_path.exists():
        import sys
        print(f"[ERROR] super_table.csv not found at {super_table_path}.",
              file=sys.stderr)
        sys.exit(1)

    super_table = pd.read_csv(super_table_path, dtype=str)
    pdb_ids = sorted(super_table["pdb_id"].dropna().unique())

    # Collect (uniprot_id, pdb_id) pairs for AF3 register
    has_uniprot = super_table["uniprot_id"].fillna("").str.strip() != ""
    af_rows = (
        super_table[has_uniprot][["uniprot_id", "pdb_id"]]
        .drop_duplicates(subset="uniprot_id")
    )
    uniprot_pairs = list(zip(af_rows["uniprot_id"], af_rows["pdb_id"]))

    # Build registers
    print(f"Building crystal register for {len(pdb_ids)} PDB IDs …")
    crystal_reg = build_crystal_register(crystal_dir, results_root, pdb_ids)
    crys_heavy_ok = crystal_reg["heavy_atom_exists"].sum()
    crys_prot_ok = crystal_reg["protonated_exists"].sum()
    print(f"  heavy_atom: {crys_heavy_ok}/{len(crystal_reg)} present")
    print(f"  protonated: {crys_prot_ok}/{len(crystal_reg)} present")

    print(f"\nBuilding AF3 register for {len(uniprot_pairs)} UniProt IDs …")
    af3_reg = build_af3_register(af3_dir, results_root, uniprot_pairs)
    af3_heavy_ok = af3_reg["af3_exists"].sum()
    af3_prot_ok = af3_reg["protonated_exists"].sum()
    print(f"  af3 prediction: {af3_heavy_ok}/{len(af3_reg)} present")
    print(f"  protonated:     {af3_prot_ok}/{len(af3_reg)} present")

    # Write
    crys_out = data_dir / "pdb_register.csv"
    af3_out = data_dir / "af3_register.csv"
    crystal_reg.to_csv(crys_out, index=False)
    af3_reg.to_csv(af3_out, index=False)
    print(f"\n[write] {crys_out}")
    print(f"[write] {af3_out}")


if __name__ == "__main__":
    main()
