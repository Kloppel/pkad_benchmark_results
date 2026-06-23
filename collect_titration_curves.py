#!/usr/bin/env python3
"""collect_titration_curves.py — Collect titration curves from KB3 results and plot them.

Scans a results directory for KB3 output pickle files, runs the built-in
kbp3.titration_analysis on each, and aggregates all titration data into:

  {output_dir}/all_titration_curves.csv     — one row per (run, residue, pH)
  {output_dir}/titration_overlay_{res}.png  — pKa distribution per residue type

Usage
-----
    python3 collect_titration_curves.py --results-dir results/pkad_benchmark_20260620
                                        [--output-dir DIR]
                                        [--skip-existing]
                                        [--no-plots]

The script is non-destructive: it only reads result files.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

HERE = Path(__file__).parent


# ---------------------------------------------------------------------------
# PKL discovery
# ---------------------------------------------------------------------------

def find_results_pkls(results_root: Path) -> list[tuple[str, Path]]:
    """Scan results_root for KB3 top-level results.pkl files.

    Pattern: results_root/{run_name}/{run_name}_karlsberg/results.pkl
    (The pac_ph_*/kbp2_res.pkl and snapshots/ files are excluded automatically
    because they don't match the *_karlsberg/results.pkl pattern.)

    Returns sorted list of (run_name, pkl_path).
    """
    found = []
    for pkl_path in sorted(results_root.glob("*/*_karlsberg/results.pkl")):
        karlsberg_dir_name = pkl_path.parent.name   # e.g. "4pti_karlsberg"
        if not karlsberg_dir_name.endswith("_karlsberg"):
            continue
        run_name = karlsberg_dir_name[: -len("_karlsberg")]
        found.append((run_name, pkl_path))
    return found


# ---------------------------------------------------------------------------
# Per-PKL processing
# ---------------------------------------------------------------------------

def process_one_pkl(
    run_name: str,
    pkl_path: Path,
    output_dir: Path,
    skip_existing: bool = True,
) -> Optional[pd.DataFrame]:
    """Run kbp3.titration_analysis on one results.pkl and return a tagged DataFrame.

    Returns None on failure (non-fatal).
    """
    run_out = output_dir / run_name
    done_marker = run_out / ".done"

    if skip_existing and done_marker.exists():
        # Try to read back the existing CSV
        csvs = list(run_out.glob("*.csv"))
        if csvs:
            try:
                df = pd.read_csv(csvs[0])
                df["run_name"] = run_name
                return df
            except Exception:
                pass

    run_out.mkdir(parents=True, exist_ok=True)
    print(f"  [{run_name}] running titration_analysis …", end=" ")

    try:
        from kbp3.titration_analysis import run_analysis
        run_analysis(
            source=str(pkl_path),
            output_dir=str(run_out),
            no_per_residue_plots=False,
        )
    except Exception as exc:
        print(f"FAILED ({exc})")
        return None

    # Find the output CSV (filename is determined by titration_analysis internals)
    csvs = sorted(run_out.glob("*.csv"))
    if not csvs:
        print("FAILED (no CSV output)")
        return None

    try:
        df = pd.read_csv(csvs[0])
        df["run_name"] = run_name
        done_marker.touch()
        print(f"ok ({len(df)} rows, {csvs[0].name})")
        return df
    except Exception as exc:
        print(f"FAILED reading CSV: {exc}")
        return None


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate_all(records: list[pd.DataFrame], out_csv: Path) -> pd.DataFrame:
    """Concatenate per-structure DataFrames and write the combined CSV."""
    combined = pd.concat(records, ignore_index=True)
    combined.to_csv(out_csv, index=False)
    print(f"[write] {out_csv}  ({len(combined)} rows, {combined['run_name'].nunique()} runs)")
    return combined


# ---------------------------------------------------------------------------
# Overlay plots
# ---------------------------------------------------------------------------

def plot_overlay(df: pd.DataFrame, out_dir: Path) -> None:
    """Plot pKa distribution per residue type (boxplot + strip).

    Looks for a 'pka' or 'pKa' column in df.  Also looks for 'resname' or
    'residue_name' or 'residue'.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[plots] matplotlib not available — skipping plots.")
        return

    # Find pKa and residue-name columns (titration_analysis output may vary)
    pka_col = next((c for c in df.columns if c.lower() in ("pka", "pka_value", "value")), None)
    res_col = next(
        (c for c in df.columns if c.lower() in ("resname", "residue_name", "residue")),
        None,
    )
    if pka_col is None or res_col is None:
        print(f"[plots] Cannot find pKa/residue columns in combined CSV "
              f"(columns: {list(df.columns)}). Skipping.")
        return

    pka_data = pd.to_numeric(df[pka_col], errors="coerce")
    df = df.copy()
    df["_pka"] = pka_data
    df["_resname"] = df[res_col].str.strip().str.upper()

    for resname, grp in df.groupby("_resname"):
        vals = grp["_pka"].dropna()
        if len(vals) < 2:
            continue

        fig, ax = plt.subplots(figsize=(4, 5))
        ax.boxplot(vals, positions=[0], widths=0.4,
                   patch_artist=True,
                   boxprops=dict(facecolor="lightblue"))
        # Strip plot (jittered)
        import random
        jitter = [random.gauss(0, 0.05) for _ in range(len(vals))]
        ax.scatter(jitter, vals, alpha=0.5, s=10, color="navy", zorder=3)

        ax.set_title(f"{resname}  (n={len(vals)})")
        ax.set_ylabel("pKa")
        ax.set_xticks([])
        ax.set_xlabel(resname)

        out_path = out_dir / f"titration_overlay_{resname}.png"
        fig.savefig(out_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"[plot] {out_path.name}  (n={len(vals)})")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True,
                        help="Root directory of a pkad_benchmark KB3 run.")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Output directory for CSV and plots "
                             "(default: {results-dir}/titration_analysis/).")
    parser.add_argument("--skip-existing", action="store_true", default=True,
                        help="Skip runs that already have a .done marker.")
    parser.add_argument("--no-skip-existing", dest="skip_existing", action="store_false")
    parser.add_argument("--no-plots", action="store_true",
                        help="Do not generate overlay plots.")
    args = parser.parse_args()

    results_root: Path = args.results_dir
    if not results_root.exists():
        print(f"[ERROR] results-dir not found: {results_root}", file=sys.stderr)
        sys.exit(1)

    output_dir = args.output_dir or (results_root / "titration_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Discover PKL files
    pkls = find_results_pkls(results_root)
    if not pkls:
        print(f"[WARN] No results.pkl files found in {results_root}")
        print("       Pattern: {results_root}/*/{run_name}_karlsberg/results.pkl")
        sys.exit(0)

    print(f"Found {len(pkls)} result PKL files in {results_root}")

    # Process each PKL
    records = []
    for run_name, pkl_path in pkls:
        df = process_one_pkl(run_name, pkl_path, output_dir, args.skip_existing)
        if df is not None:
            records.append(df)

    if not records:
        print("[WARN] No titration data extracted.")
        sys.exit(0)

    # Aggregate
    all_csv = output_dir / "all_titration_curves.csv"
    combined = aggregate_all(records, all_csv)

    # Plots
    if not args.no_plots:
        print("\nGenerating pKa overlay plots …")
        plot_overlay(combined, output_dir)

    print(f"\nDone. Outputs in {output_dir}/")


if __name__ == "__main__":
    main()
