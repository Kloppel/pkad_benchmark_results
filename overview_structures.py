#!/usr/bin/env python3
"""overview_structures.py — Collaborator-facing overview of all benchmark structures.

Reads pkad_data/super_table.csv and checks which crystal / AF3 PDB files are
present on disk.  Prints a formatted ASCII table and writes
pkad_data/overview.csv.

The overview shows, for each unique PDB ID:
  - Which datasets include it (PKAD2 / PKAD3 / PKAD-R)
  - Number of experimental pKa values
  - Residue types with experimental data (comma-separated)
  - Whether the crystal PDB and AF3 prediction are available

Usage
-----
    python3 overview_structures.py [--data-dir pkad_data] [--no-print]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
DEFAULT_DATA_DIR = HERE / "pkad_data"


# ---------------------------------------------------------------------------
# Build overview
# ---------------------------------------------------------------------------

def build_overview(
    super_table: pd.DataFrame,
    crystal_pdb_dir: Path,
    af3_pdb_dir: Path,
) -> pd.DataFrame:
    """Return one row per unique pdb_id with summary columns."""
    rows = []

    for pdb_id, grp in super_table.groupby("pdb_id", sort=True):
        in_pkad2 = grp["in_pkad2"].astype(str).isin(["True", "true", "1"]).any()
        in_pkad3 = grp["in_pkad3"].astype(str).isin(["True", "true", "1"]).any()
        in_pkadR = grp["in_pkadR"].astype(str).isin(["True", "true", "1"]).any()

        datasets = []
        if in_pkad2:
            datasets.append("PKAD2")
        if in_pkad3:
            datasets.append("PKAD3")
        if in_pkadR:
            datasets.append("PKAD-R")

        n_pkas = grp["expt_pka"].dropna().count()

        residue_types = sorted(
            grp["resname"].dropna().str.strip().str.upper().unique()
        )

        uniprot_ids = (
            grp["uniprot_id"].dropna().str.strip()
            .replace("", pd.NA).dropna().unique()
        )
        uniprot_id = uniprot_ids[0] if len(uniprot_ids) > 0 else ""

        # Check PDB file existence
        crystal_path = crystal_pdb_dir / f"{pdb_id}.pdb"
        crystal_exists = crystal_path.exists()

        af3_path = af3_pdb_dir / f"AF-{uniprot_id}-F1-model_v4.pdb" if uniprot_id else None
        af3_exists = af3_path.exists() if af3_path else False

        rows.append({
            "pdb_id":              pdb_id,
            "datasets":            ",".join(datasets),
            "n_pkas":              int(n_pkas),
            "residue_types":       ",".join(residue_types),
            "uniprot_id":          uniprot_id,
            "crystal_pdb_exists":  crystal_exists,
            "crystal_pdb_path":    str(crystal_path) if crystal_exists else "",
            "af3_pdb_exists":      af3_exists,
            "af3_pdb_path":        str(af3_path) if af3_exists else "",
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Formatted ASCII table
# ---------------------------------------------------------------------------

def print_overview_table(df: pd.DataFrame) -> None:
    """Print a readable table of the overview to stdout."""
    # Column widths
    w_pdb = 6
    w_ds = 18
    w_n = 5
    w_res = 28
    w_up = 12
    w_crys = 8
    w_af3 = 8

    hdr = (
        f"{'PDB':<{w_pdb}}  {'Datasets':<{w_ds}}  {'#pKa':>{w_n}}  "
        f"{'Residue types':<{w_res}}  {'UniProt':<{w_up}}  "
        f"{'Crystal':>{w_crys}}  {'AF3':>{w_af3}}"
    )
    sep = "-" * len(hdr)

    total = len(df)
    crystal_ok = df["crystal_pdb_exists"].sum()
    af3_ok = df["af3_pdb_exists"].sum()
    n_pkad2 = df["datasets"].str.contains("PKAD2").sum()
    n_pkad3 = df["datasets"].str.contains("PKAD3").sum()
    n_pkadR = df["datasets"].str.contains("PKAD-R").sum()

    print(sep)
    print("PKAD BENCHMARK — STRUCTURE OVERVIEW")
    print(sep)
    print(f"Total unique PDB structures : {total}")
    print(f"  In PKAD2                  : {n_pkad2}")
    print(f"  In PKAD3                  : {n_pkad3}")
    print(f"  In PKAD-R                 : {n_pkadR}")
    print(f"Crystal PDB available       : {crystal_ok} / {total}")
    print(f"AF3 prediction available    : {af3_ok} / {total}")
    print(sep)
    print(hdr)
    print(sep)

    for _, row in df.iterrows():
        res_str = row["residue_types"]
        if len(res_str) > w_res:
            res_str = res_str[: w_res - 1] + "…"
        crys_mark = "YES" if row["crystal_pdb_exists"] else "---"
        af3_mark = "YES" if row["af3_pdb_exists"] else "---"
        print(
            f"{row['pdb_id']:<{w_pdb}}  "
            f"{row['datasets']:<{w_ds}}  "
            f"{row['n_pkas']:>{w_n}}  "
            f"{res_str:<{w_res}}  "
            f"{row['uniprot_id']:<{w_up}}  "
            f"{crys_mark:>{w_crys}}  "
            f"{af3_mark:>{w_af3}}"
        )
    print(sep)
    print(f"Total: {total} structures  |  "
          f"Crystal: {crystal_ok}  |  AF3: {af3_ok}")
    print(sep)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR,
                        help="Data directory containing super_table.csv.")
    parser.add_argument("--no-print", action="store_true",
                        help="Suppress ASCII table output.")
    args = parser.parse_args()

    data_dir: Path = args.data_dir
    super_table_path = data_dir / "super_table.csv"

    if not super_table_path.exists():
        import sys
        print(f"[ERROR] super_table.csv not found at {super_table_path}.",
              file=sys.stderr)
        print("        Run download_pkad_data.py first.", file=sys.stderr)
        sys.exit(1)

    super_table = pd.read_csv(super_table_path, dtype=str)
    crystal_dir = data_dir / "pdbs" / "heavy_atom"
    af3_dir = data_dir / "af3_pdbs"

    overview = build_overview(super_table, crystal_dir, af3_dir)

    out_path = data_dir / "overview.csv"
    overview.to_csv(out_path, index=False)
    print(f"[write] {out_path}  ({len(overview)} structures)")

    if not args.no_print:
        print()
        print_overview_table(overview)


if __name__ == "__main__":
    main()
