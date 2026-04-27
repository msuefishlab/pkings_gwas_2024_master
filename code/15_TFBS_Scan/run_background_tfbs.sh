#!/bin/bash
## run_background_tfbs.sh
## Runs get_tfbs_background.sh for every window in the background manifest.
##
## Usage:
##   bash code/15_TFBS_Scan/run_background_tfbs.sh [MANIFEST] [N_PARALLEL]
##
##   MANIFEST    Path to background_windows_manifest.tsv (default: see below)
##   N_PARALLEL  Number of parallel jobs (default: 4)
##
## Example:
##   bash code/15_TFBS_Scan/run_background_tfbs.sh \
##     output_data/15_TFBS_Scan/background/background_windows_manifest.tsv 4

set -eo pipefail

root="$(git rev-parse --show-toplevel)"

MANIFEST="${1:-${root}/output_data/15_TFBS_Scan/background/background_windows_manifest.tsv}"
N_PARALLEL="${2:-4}"

if [[ ! -f "$MANIFEST" ]]; then
    echo "Error: manifest not found: $MANIFEST" >&2
    echo "Run sample_background_windows.py first." >&2
    exit 1
fi

LOG_DIR="${root}/output_data/15_TFBS_Scan/background/logs"
mkdir -p "$LOG_DIR"

SCRIPT="${root}/code/15_TFBS_Scan/get_tfbs_background.sh"

# Count windows to run (minus header)
N_TOTAL=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')
echo "[INFO] Running $N_TOTAL background windows with N_PARALLEL=$N_PARALLEL"
echo "[INFO] Logs: $LOG_DIR"

n_queued=0
n_skipped=0
n_done=0
n_running=0

while IFS=$'\t' read -r window_id locus chr start end size_bp sample_id region; do
    region="${region%$'\r'}"  # strip any \r from Windows-style line endings
    if [[ "$sample_id" == "UNSET" ]]; then
        echo "[WARN] Skipping $window_id — sample_id is UNSET in manifest" >&2
        n_skipped=$(( n_skipped + 1 ))
        continue
    fi
    n_queued=$(( n_queued + 1 ))
    logfile="${LOG_DIR}/${window_id}.log"

    bash "$SCRIPT" "$region" "$window_id" "$sample_id" > "$logfile" 2>&1 &
    n_running=$(( n_running + 1 ))

    # When we hit N_PARALLEL running jobs, wait for the batch to finish
    if (( n_running >= N_PARALLEL )); then
        wait
        n_done=$(( n_done + n_running ))
        echo "[INFO] Progress: $n_done / $n_queued windows complete"
        n_running=0
    fi
done < <(tail -n +2 "$MANIFEST")

# Wait for any remaining jobs
wait
n_done=$(( n_done + n_running ))

echo "[INFO] $n_skipped windows skipped (UNSET sample_id)"
echo "[INFO] $n_done / $n_queued windows complete."
echo "[INFO] Results in: ${root}/output_data/15_TFBS_Scan/background/"
