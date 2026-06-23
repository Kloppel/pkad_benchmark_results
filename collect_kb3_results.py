#!/usr/bin/env python3
"""collect_kb3_results.py — Collect KB3 pKa results and compare to experiment.

Reads:
  pkad_data/experimental_pkas.tsv           — pdb_id / chain / resname / resid / expt_pka
  results/<run_dir>/<pdbid>_results.txt     — KB3 output per structure

Format of *_results.txt:
    RESNAME-RESID_CHAIN: pKa

Writes:
  pkad_data/benchmark_summary_<run_dir>.tsv — one row per matched/missing residue

Prints per-structure and overall RMSE/MAE/N to stdout.

Usage:
    python3 collect_kb3_results.py --results-dir results/pkad_benchmark_20260619
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).parent
DEFAULT_DATA_DIR = HERE / "pkad_data"
DEFAULT_RESULTS_ROOT = HERE / "results"

# Residue names Karlsberg reports that map to the experimental TSV names
RESNAME_MAP = {
    "NTE": None,   # N-terminus — skip, no experimental entry
    "CTE": None,   # C-terminus — skip
}

# Titratable residue names kept for matching
TITRATABLE = {"ASP", "GLU", "HIS", "LYS", "ARG", "CYS", "TYR", "NTE", "CTE"}

_RESULT_LINE = re.compile(
    r"^\s*([A-Z]+)-(\d+)_([A-Z]):\s*([-\d.]+)\s*$"
)


def parse_results_txt(txt_path: Path) -> list[dict]:
    """Parse a KB3 *_results.txt file.

    Each line:  RESNAME-RESID_CHAIN: pKa
    Returns list of dicts: resname, resid, chain, calc_pka.
    """
    records = []
    for line in txt_path.read_text(errors="replace").splitlines():
        m = _RESULT_LINE.match(line)
        if not m:
            continue
        resname, resid, chain, pka_str = m.groups()
        records.append({
            "resname": resname,
            "resid": resid,
            "chain": chain,
            "calc_pka": float(pka_str),
        })
    return records


def load_experimental(tsv_path: Path) -> dict[tuple, dict]:
    """Return dict keyed by (pdb_id, resname, resid) → row dict."""
    data = {}
    with open(tsv_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (
                row["pdb_id"].strip().lower(),
                row["resname"].strip().upper(),
                row["resid"].strip(),
            )
            data[key] = row
    return data


def compute_stats(deltas: list[float]) -> dict:
    n = len(deltas)
    if n == 0:
        return {"N": 0, "RMSE": None, "MAE": None}
    mae = sum(abs(d) for d in deltas) / n
    rmse = math.sqrt(sum(d * d for d in deltas) / n)
    return {"N": n, "RMSE": round(rmse, 3), "MAE": round(mae, 3)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results-dir", type=Path, required=True,
                        help="Run directory (contains <pdbid>_results.txt files).")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    args = parser.parse_args()

    results_dir: Path = args.results_dir
    data_dir: Path = args.data_dir

    if not results_dir.exists():
        sys.exit(f"[ERROR] Results directory not found: {results_dir}")

    exp_tsv = data_dir / "experimental_pkas.tsv"
    if not exp_tsv.exists():
        sys.exit(f"[ERROR] experimental_pkas.tsv not found at {exp_tsv}")

    experimental = load_experimental(exp_tsv)

    # Collect all pdb_ids from experimental data
    all_pdb_ids = sorted({k[0] for k in experimental})

    summary_rows: list[dict] = []
    all_deltas: list[float] = []
    per_pdb_stats: list[tuple] = []

    for pdb_id in all_pdb_ids:
        txt_path = results_dir / f"{pdb_id}_results.txt"
        if not txt_path.exists():
            # Structure was not run
            pdb_keys = [k for k in experimental if k[0] == pdb_id]
            for key in pdb_keys:
                row = experimental[key]
                summary_rows.append({
                    "pdb_id": pdb_id,
                    "chain": row.get("chain", ""),
                    "resname": key[1],
                    "resid": key[2],
                    "expt_pka": row["expt_pka"],
                    "calc_pka": "",
                    "delta": "",
                    "status": "not_run",
                })
            continue

        calc_records = parse_results_txt(txt_path)
        # Index by (resname, resid) — chain may differ between experimental and calc
        calc_index: dict[tuple, float] = {}
        for rec in calc_records:
            calc_index[(rec["resname"], rec["resid"])] = rec["calc_pka"]

        pdb_keys = [k for k in experimental if k[0] == pdb_id]
        pdb_deltas: list[float] = []

        for key in sorted(pdb_keys, key=lambda k: (k[1], int(k[2]) if k[2].isdigit() else 0)):
            row = experimental[key]
            _, resname, resid = key
            expt_pka_str = row["expt_pka"].strip()
            try:
                expt_pka = float(expt_pka_str)
            except ValueError:
                continue

            calc_pka = calc_index.get((resname, resid))
            if calc_pka is None:
                status = "missing"
                delta_str = ""
            else:
                delta = calc_pka - expt_pka
                delta_str = f"{delta:+.2f}"
                pdb_deltas.append(delta)
                all_deltas.append(delta)
                status = "matched"

            summary_rows.append({
                "pdb_id": pdb_id,
                "chain": row.get("chain", ""),
                "resname": resname,
                "resid": resid,
                "expt_pka": expt_pka_str,
                "calc_pka": f"{calc_pka:.2f}" if calc_pka is not None else "",
                "delta": delta_str,
                "status": status,
            })

        if pdb_deltas:
            stats = compute_stats(pdb_deltas)
            per_pdb_stats.append((pdb_id, stats))

    # Write TSV
    run_name = results_dir.name
    out_path = data_dir / f"benchmark_summary_{run_name}.tsv"
    fieldnames = ["pdb_id", "chain", "resname", "resid",
                  "expt_pka", "calc_pka", "delta", "status"]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"[write] {out_path}  ({len(summary_rows)} rows)")

    # Count run status
    n_matched = sum(1 for r in summary_rows if r["status"] == "matched")
    n_missing = sum(1 for r in summary_rows if r["status"] == "missing")
    n_not_run = sum(1 for r in summary_rows if r["status"] == "not_run")
    n_structs_run = sum(1 for _, s in per_pdb_stats)

    print(f"\nStructures: {n_structs_run}/{len(all_pdb_ids)} run")
    print(f"Residues:   {n_matched} matched  |  {n_missing} missing  |  {n_not_run} not run")

    print(f"\n{'PDB':<8}  {'N':>4}  {'RMSE':>6}  {'MAE':>6}")
    print("-" * 30)
    for pdb_id, stats in per_pdb_stats:
        rmse = f"{stats['RMSE']:.3f}" if stats["RMSE"] is not None else "—"
        mae  = f"{stats['MAE']:.3f}"  if stats["MAE"]  is not None else "—"
        print(f"{pdb_id:<8}  {stats['N']:>4}  {rmse:>6}  {mae:>6}")

    overall = compute_stats(all_deltas)
    rmse = f"{overall['RMSE']:.3f}" if overall["RMSE"] is not None else "—"
    mae  = f"{overall['MAE']:.3f}"  if overall["MAE"]  is not None else "—"
    print("-" * 30)
    print(f"{'OVERALL':<8}  {overall['N']:>4}  {rmse:>6}  {mae:>6}")


if __name__ == "__main__":
    main()
