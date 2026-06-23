#!/usr/bin/env bash
# delete_synced_results.sh
#
# Deletes local *_karlsberg/ folders where:
#   1. The calculation succeeded (matching *_results.txt is non-empty)
#   2. The same folder exists on flashnas2 at the mirror path
#
# Usage:
#   ./delete_synced_results.sh            # dry-run (no deletions)
#   ./delete_synced_results.sh --delete   # actually delete

set -euo pipefail

LOCAL_RESULTS="$(cd "$(dirname "$0")/results" && pwd)"
REMOTE_RESULTS="/media/flashnas2/jones/Karlsberg3/benchmarking/pkad_benchmark/results"
DRY_RUN=true

if [[ "${1:-}" == "--delete" ]]; then
    DRY_RUN=false
fi

if [[ ! -d "$REMOTE_RESULTS" ]]; then
    echo "[ERROR] Remote path not mounted or missing: $REMOTE_RESULTS" >&2
    exit 1
fi

total=0
deleted=0
skipped_no_remote=0
skipped_failed=0

for run_dir in "$LOCAL_RESULTS"/pkad_benchmark_*/; do
    run_name="$(basename "$run_dir")"

    # skip the combined summary dir
    [[ "$run_name" == *_combined ]] && continue

    for karlsberg_dir in "$run_dir"*_karlsberg/; do
        [[ -d "$karlsberg_dir" ]] || continue
        total=$((total + 1))

        # derive the pdb_id from the folder name (strip _karlsberg suffix)
        folder_name="$(basename "$karlsberg_dir")"
        pdb_id="${folder_name%_karlsberg}"

        results_txt="$run_dir/${pdb_id}_results.txt"

        # check success: results.txt must exist and be non-empty
        if [[ ! -s "$results_txt" ]]; then
            skipped_failed=$((skipped_failed + 1))
            echo "[skip-failed ]  $run_name/$folder_name"
            continue
        fi

        # check remote mirror exists
        remote_dir="$REMOTE_RESULTS/$run_name/$folder_name"
        if [[ ! -d "$remote_dir" ]]; then
            skipped_no_remote=$((skipped_no_remote + 1))
            echo "[skip-no-sync]  $run_name/$folder_name"
            continue
        fi

        deleted=$((deleted + 1))
        if $DRY_RUN; then
            echo "[would delete]  $run_name/$folder_name"
        else
            echo "[deleting    ]  $run_name/$folder_name"
            rm -rf "$karlsberg_dir"
        fi
    done
done

echo ""
echo "Summary: $total folders scanned"
echo "  $deleted $(if $DRY_RUN; then echo 'would be deleted'; else echo 'deleted'; fi)"
echo "  $skipped_failed skipped (calculation failed / no pKa output)"
echo "  $skipped_no_remote skipped (not found on flashnas2)"

if $DRY_RUN && [[ $deleted -gt 0 ]]; then
    echo ""
    echo "Dry run complete. Re-run with --delete to actually remove the folders."
fi
