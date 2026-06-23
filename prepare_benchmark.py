#!/usr/bin/env python3
"""prepare_benchmark.py — Build the job manifest for the PKAD2/PKAD3/PKAD-R benchmark.

Reads:
  pkad_data/super_table.csv          — produced by download_pkad_data.py
  pkad_data/pdbs/heavy_atom/         — crystal PDB files
  pkad_data/af3_pdbs/                — AlphaFold DB prediction PDB files

Writes:
  pkad_data/job_manifest.tsv         — one row per (pdb_id, structure_type)
      pdb_id  structure_type  pdb_path  datasets  n_experimental_pkas  uniprot_id

  pkad_data/experimental_pkas.tsv    — full flat table of all experimental pKa values
      source_dataset  pdb_id  chain  resname  resid  expt_pka

Usage
-----
    python3 prepare_benchmark.py [--data-dir pkad_data]

Run after download_pkad_data.py.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
DEFAULT_DATA_DIR = HERE / "pkad_data"


# ---------------------------------------------------------------------------
# Manifest building
# ---------------------------------------------------------------------------

def build_manifest(
    super_table: pd.DataFrame,
    crystal_pdb_dir: Path,
    af3_pdb_dir: Path,
) -> pd.DataFrame:
    """Build the job manifest with one row per (pdb_id, structure_type).

    structure_type is 'crystal' for original crystal structures and 'af3' for
    AlphaFold DB prediction structures.  Only pdb_ids with a UniProt ID (from
    PKAD3) get an AF3 row.  Each pair appears exactly once.
    """
    rows = []

    # --- Group by pdb_id ---
    grouped = super_table.groupby("pdb_id", sort=True)

    for pdb_id, grp in grouped:
        datasets = []
        if grp["in_pkad2"].astype(str).isin(["True", "true", "1"]).any():
            datasets.append("PKAD2")
        if grp["in_pkad3"].astype(str).isin(["True", "true", "1"]).any():
            datasets.append("PKAD3")
        if grp["in_pkadR"].astype(str).isin(["True", "true", "1"]).any():
            datasets.append("PKAD-R")

        n_pkas = grp["expt_pka"].dropna().count()
        datasets_str = ",".join(datasets)

        # Unique UniProt ID for this pdb_id (first non-empty)
        uniprot_ids = grp["uniprot_id"].dropna().str.strip()
        uniprot_ids = uniprot_ids[uniprot_ids != ""]
        uniprot_id = uniprot_ids.iloc[0] if len(uniprot_ids) > 0 else ""

        # Crystal row
        crystal_path = crystal_pdb_dir / f"{pdb_id}.pdb"
        rows.append({
            "pdb_id":               pdb_id,
            "structure_type":       "crystal",
            "pdb_path":             str(crystal_path),
            "pdb_exists":           crystal_path.exists(),
            "datasets":             datasets_str,
            "n_experimental_pkas":  n_pkas,
            "uniprot_id":           "",
        })

        # AF3 row — only if this pdb_id has a UniProt ID
        if uniprot_id:
            af3_path = af3_pdb_dir / f"AF-{uniprot_id}-F1-model_v4.pdb"
            rows.append({
                "pdb_id":               pdb_id,
                "structure_type":       "af3",
                "pdb_path":             str(af3_path),
                "pdb_exists":           af3_path.exists(),
                "datasets":             datasets_str,
                "n_experimental_pkas":  n_pkas,
                "uniprot_id":           uniprot_id,
            })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR,
                        help="Directory containing super_table.csv and PDB subdirs.")
    parser.add_argument("--include-missing-pdbs", action="store_true",
                        help="Include manifest rows even when the PDB file does not exist.")
    args = parser.parse_args()

    data_dir: Path = args.data_dir
    crystal_dir = data_dir / "pdbs" / "heavy_atom"
    af3_dir = data_dir / "af3_pdbs"

    super_table_path = data_dir / "super_table.csv"
    if not super_table_path.exists():
        print(f"[ERROR] super_table.csv not found at {super_table_path}.",
              file=sys.stderr)
        print("        Run download_pkad_data.py first.", file=sys.stderr)
        sys.exit(1)

    super_table = pd.read_csv(super_table_path, dtype=str)
    print(f"[super_table] {len(super_table)} rows, "
          f"{super_table['pdb_id'].nunique()} unique PDB IDs")

    manifest_full = build_manifest(super_table, crystal_dir, af3_dir)

    # Filter rows where PDB file does not exist (unless --include-missing-pdbs)
    missing_crystal = manifest_full[
        (manifest_full["structure_type"] == "crystal") &
        (~manifest_full["pdb_exists"])
    ]
    missing_af3 = manifest_full[
        (manifest_full["structure_type"] == "af3") &
        (~manifest_full["pdb_exists"])
    ]

    if not args.include_missing_pdbs:
        manifest = manifest_full[manifest_full["pdb_exists"]].copy()
    else:
        manifest = manifest_full.copy()

    n_crystal = (manifest["structure_type"] == "crystal").sum()
    n_af3 = (manifest["structure_type"] == "af3").sum()
    print(f"Manifest rows: {len(manifest)} total  "
          f"({n_crystal} crystal, {n_af3} AF3)")
    if not args.include_missing_pdbs:
        print(f"Excluded (PDB file missing): "
              f"{len(missing_crystal)} crystal, {len(missing_af3)} AF3")

    # Write manifest (without pdb_exists helper column)
    manifest_out = manifest.drop(columns=["pdb_exists"])
    out_path = data_dir / "job_manifest.tsv"
    manifest_out.to_csv(out_path, sep="\t", index=False)
    print(f"[write] {out_path}")

    # Write experimental_pkas.tsv (flat table for collect_results.py)
    exp_cols = ["source_dataset", "pdb_id", "chain", "resname", "resid", "expt_pka"]
    available = [c for c in exp_cols if c in super_table.columns]
    exp_df = (
        super_table[available]
        .rename(columns={"source_dataset": "dataset"})
        .dropna(subset=["expt_pka"])
        .copy()
    )
    exp_path = data_dir / "experimental_pkas.tsv"
    exp_df.to_csv(exp_path, sep="\t", index=False)
    print(f"[write] {exp_path}  ({len(exp_df)} rows)")

    print("\nNext step:")
    print(f"  python3 queue_all.py")
    print(f"  # or dry-run:")
    print(f"  python3 queue_all.py --dry-run")


if __name__ == "__main__":
    main()
