#!/usr/bin/env python3
"""run_pkad_benchmark.py — Top-level entry point for the PKAD benchmark.

Orchestrates the full pipeline:
  Step 1  download_pkad_data.py    — read PKAD tables, download crystal + AF3 PDBs
  Step 2  prepare_benchmark.py     — build job manifest (crystal + AF3 rows)
  Step 3  queue_all.py             — submit SLURM job array
  Step 4  (post-run) collect_results.py    — compute RMSE/MAE vs. experiment
  Step 5  (post-run) collect_titration_curves.py — aggregate titration curves + plots

Run all setup steps (1 + 2) and print submission instructions:
    python3 run_pkad_benchmark.py

Run a specific step:
    python3 run_pkad_benchmark.py --step 1    # download only
    python3 run_pkad_benchmark.py --step 2    # prepare manifest only
    python3 run_pkad_benchmark.py --step 3    # queue to cluster (dry-run: add --dry-run)
    python3 run_pkad_benchmark.py --step 4 --results-dir results/pkad_benchmark_YYYYMMDD
    python3 run_pkad_benchmark.py --step 5 --results-dir results/pkad_benchmark_YYYYMMDD

Generate CHARMM build scripts for completed runs:
    python3 run_pkad_benchmark.py --step 3 --generate-charmm-scripts \
        --results-root results/pkad_benchmark_YYYYMMDD
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent


def run_step(script: str, extra_args: list[str] | None = None) -> None:
    cmd = [sys.executable, str(HERE / script)] + (extra_args or [])
    print(f"\n{'='*60}")
    print(f"Running: {' '.join(str(a) for a in cmd)}")
    print("=" * 60)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"[ERROR] {script} exited with code {result.returncode}", file=sys.stderr)
        sys.exit(result.returncode)


def print_submission_instructions() -> None:
    manifest = HERE / "pkad_data" / "job_manifest.tsv"
    if not manifest.exists():
        return
    n = sum(1 for _ in open(manifest)) - 1  # subtract header
    n_crystal = 0
    n_af3 = 0
    with open(manifest) as fh:
        next(fh)  # skip header
        for line in fh:
            if "\tcrystal\t" in line:
                n_crystal += 1
            elif "\taf3\t" in line:
                n_af3 += 1

    print(f"""
{'='*60}
Step 3: Submit to cluster
{'='*60}
  cd {HERE}
  mkdir -p logs

  # Total jobs: {n} ({n_crystal} crystal + {n_af3} AF3)
  python3 queue_all.py --dry-run   # preview
  python3 queue_all.py             # submit

  # Monitor progress
  squeue -u $USER -n pkad_kb3

  # Or set a specific results root (date is auto-generated otherwise):
  KB3_RESULTS_ROOT={HERE}/results/pkad_benchmark_$(date -u +%Y%m%d) \\
      python3 queue_all.py

Step 4: Collect pKa comparison statistics (after all jobs finish)
  python3 collect_results.py --results-dir results/pkad_benchmark_<date>

Step 5: Collect titration curves and plots
  python3 collect_titration_curves.py --results-dir results/pkad_benchmark_<date>

Step 6 (optional): Generate CHARMM build scripts for completed runs
  python3 queue_all.py --generate-charmm-scripts \\
      --results-root results/pkad_benchmark_<date>
""")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=int, choices=[1, 2, 3, 4, 5],
                        help="Run only this step.")
    parser.add_argument("--results-dir", type=Path,
                        help="Results directory (required for steps 4 and 5).")
    parser.add_argument("--results-root", type=Path,
                        help="Results root for CHARMM script generation (step 3).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Passed to queue_all.py (step 3).")
    parser.add_argument("--generate-charmm-scripts", action="store_true",
                        help="Passed to queue_all.py (step 3).")
    parser.add_argument("--no-af3", action="store_true",
                        help="Passed to download_pkad_data.py (step 1).")
    parser.add_argument("--skip-existing", action="store_true", default=True)
    parser.add_argument("--no-skip-existing", dest="skip_existing", action="store_false")
    args = parser.parse_args()

    dl_args: list[str] = []
    if args.no_af3:
        dl_args.append("--no-af3")
    if not args.skip_existing:
        dl_args.append("--no-skip-existing")

    if args.step == 1:
        run_step("download_pkad_data.py", dl_args)

    elif args.step == 2:
        run_step("prepare_benchmark.py")

    elif args.step == 3:
        q_args: list[str] = []
        if args.dry_run:
            q_args.append("--dry-run")
        if args.generate_charmm_scripts:
            q_args.append("--generate-charmm-scripts")
        if args.results_root:
            q_args += ["--results-root", str(args.results_root)]
        run_step("queue_all.py", q_args)

    elif args.step == 4:
        if not args.results_dir:
            parser.error("--results-dir is required for --step 4")
        run_step("collect_results.py",
                 ["--results-dir", str(args.results_dir)])

    elif args.step == 5:
        if not args.results_dir:
            parser.error("--results-dir is required for --step 5")
        run_step("collect_titration_curves.py",
                 ["--results-dir", str(args.results_dir)])

    else:
        # Default: steps 1 + 2
        run_step("download_pkad_data.py", dl_args)
        run_step("prepare_benchmark.py")
        print_submission_instructions()


if __name__ == "__main__":
    main()
