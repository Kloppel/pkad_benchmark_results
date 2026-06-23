#!/usr/bin/env python3
"""collect_results.py — Aggregate per-structure KB3 results and compare to experiment.

Reads:
  pkad_data/experimental_pkas.tsv    — produced by prepare_benchmark.py
  results/<run_dir>/<pdbid>/         — per-structure output directories from run_kb3.py

Writes to pkad_data/benchmark_summary_<run_dir>.tsv:
  pdb_id  chain  resname  resid  dataset  expt_pka  calc_pka  delta  status

And prints per-dataset RMSE / MAE / N to stdout.

Usage
-----
    python3 collect_results.py --results-dir results/pkad_benchmark_20260620

The script is non-destructive: it only reads result files.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

HERE = Path(__file__).parent
DEFAULT_DATA_DIR = HERE / "pkad_data"


# ---------------------------------------------------------------------------
# Load experimental pKas
# ---------------------------------------------------------------------------

def load_experimental(tsv_path: Path) -> dict[tuple[str, str, str, str, str], list[dict]]:
    """Return dict keyed by (pdb_id, chain, resname, resid, dataset)."""
    data: dict = defaultdict(list)
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = (
                row["pdb_id"].lower(),
                row.get("chain", "").strip(),
                row.get("resname", "").strip(),
                row.get("resid", "").strip(),
                row["dataset"],
            )
            data[key].append(row)
    return data


# ---------------------------------------------------------------------------
# Load KB3 results
# ---------------------------------------------------------------------------

def _parse_results_txt(txt_path: Path) -> list[dict]:
    """Parse a *_results.txt file written by run_kb3.py / JobHelper.write_results.

    The format is not rigidly defined, so we accept two common layouts:
      1. Lines like:  CHAIN  RESNAME  RESID  pKa=X.XX
      2. Lines like:  RESNAME RESID CHAIN pKa X.XX
    Returns list of dicts with keys: chain, resname, resid, calc_pka.
    """
    results = []
    for line in txt_path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Look for a pKa value in the line
        pka: Optional[float] = None
        for token in line.split():
            if token.lower().startswith("pka="):
                try:
                    pka = float(token.split("=", 1)[1])
                except ValueError:
                    pass
            elif token.replace(".", "", 1).lstrip("-").isdigit():
                try:
                    pka = float(token)
                except ValueError:
                    pass
        if pka is None:
            continue
        parts = line.split()
        # Try to extract resname / resid / chain from the line tokens
        chain = resname = resid = ""
        for i, tok in enumerate(parts):
            if len(tok) == 1 and tok.isalpha() and tok.isupper():
                chain = tok
            elif len(tok) == 3 and tok.isalpha() and tok.isupper():
                resname = tok
            elif tok.isdigit():
                resid = tok
        results.append({"chain": chain, "resname": resname, "resid": resid,
                        "calc_pka": pka})
    return results


def load_kb3_results(results_dir: Path, pdb_id: str) -> list[dict]:
    """Load all calculated pKas for one structure from its results directory."""
    pdb_dir = results_dir / pdb_id
    if not pdb_dir.exists():
        return []
    calc: list[dict] = []
    # Prefer the pkl-derived titration analysis CSV if present
    for csv_path in sorted(pdb_dir.glob("*.csv")):
        try:
            with open(csv_path, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    if "pka" in (row.get("metric", "") or "").lower():
                        calc.append({
                            "chain": row.get("chain", ""),
                            "resname": row.get("resname", ""),
                            "resid": row.get("resid", ""),
                            "calc_pka": float(row.get("value", "nan")),
                        })
        except Exception:
            pass
    if calc:
        return calc
    # Fallback: parse *_results.txt files
    for txt_path in sorted(pdb_dir.glob("*_results.txt")):
        calc.extend(_parse_results_txt(txt_path))
    return calc


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def compute_stats(deltas: list[float]) -> dict:
    if not deltas:
        return {"N": 0, "RMSE": None, "MAE": None}
    n = len(deltas)
    mae = sum(abs(d) for d in deltas) / n
    rmse = math.sqrt(sum(d * d for d in deltas) / n)
    return {"N": n, "RMSE": round(rmse, 3), "MAE": round(mae, 3)}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True,
                        help="Root directory of a pkad_benchmark run "
                             "(contains one sub-directory per PDB ID).")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR,
                        help="Directory containing experimental_pkas.tsv.")
    args = parser.parse_args()

    results_dir: Path = args.results_dir
    data_dir: Path = args.data_dir

    exp_tsv = data_dir / "experimental_pkas.tsv"
    if not exp_tsv.exists():
        print(f"[ERROR] experimental_pkas.tsv not found at {exp_tsv}.", file=sys.stderr)
        print("        Run prepare_benchmark.py first.", file=sys.stderr)
        sys.exit(1)

    experimental = load_experimental(exp_tsv)
    pdb_ids = sorted({k[0] for k in experimental})

    summary_rows: list[dict] = []
    by_dataset: dict[str, list[float]] = defaultdict(list)

    for pdb_id in pdb_ids:
        calc_list = load_kb3_results(results_dir, pdb_id)
        calc_by_resid: dict[tuple[str, str], float] = {}
        for c in calc_list:
            key = (c.get("resname", ""), c.get("resid", ""))
            calc_by_resid[key] = c["calc_pka"]

        for key, rows in experimental.items():
            if key[0] != pdb_id:
                continue
            _, chain, resname, resid, dataset = key
            expt_pka_str = rows[0].get("expt_pka", "")
            try:
                expt_pka = float(expt_pka_str)
            except ValueError:
                continue

            calc_key = (resname, resid)
            calc_pka = calc_by_resid.get(calc_key)
            delta = (calc_pka - expt_pka) if calc_pka is not None else None
            status = "matched" if calc_pka is not None else "missing"

            summary_rows.append({
                "pdb_id": pdb_id,
                "chain": chain,
                "resname": resname,
                "resid": resid,
                "dataset": dataset,
                "expt_pka": expt_pka_str,
                "calc_pka": f"{calc_pka:.2f}" if calc_pka is not None else "",
                "delta": f"{delta:.2f}" if delta is not None else "",
                "status": status,
            })
            if delta is not None:
                by_dataset[dataset].append(delta)

    # Write summary TSV
    run_name = results_dir.name
    out_path = data_dir / f"benchmark_summary_{run_name}.tsv"
    fieldnames = ["pdb_id", "chain", "resname", "resid", "dataset",
                  "expt_pka", "calc_pka", "delta", "status"]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"[write] {out_path}  ({len(summary_rows)} rows)")

    # Print per-dataset statistics
    print("\nPer-dataset statistics (matched residues only):")
    print(f"{'Dataset':<12}  {'N':>6}  {'RMSE':>6}  {'MAE':>6}")
    print("-" * 36)
    all_deltas: list[float] = []
    for ds in sorted(by_dataset):
        stats = compute_stats(by_dataset[ds])
        rmse_str = f"{stats['RMSE']:.3f}" if stats["RMSE"] is not None else "—"
        mae_str = f"{stats['MAE']:.3f}" if stats["MAE"] is not None else "—"
        print(f"{ds:<12}  {stats['N']:>6}  {rmse_str:>6}  {mae_str:>6}")
        all_deltas.extend(by_dataset[ds])

    if all_deltas:
        combined = compute_stats(all_deltas)
        rmse_str = f"{combined['RMSE']:.3f}" if combined["RMSE"] is not None else "—"
        mae_str = f"{combined['MAE']:.3f}" if combined["MAE"] is not None else "—"
        print("-" * 36)
        print(f"{'COMBINED':<12}  {combined['N']:>6}  {rmse_str:>6}  {mae_str:>6}")


if __name__ == "__main__":
    main()
