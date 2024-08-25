#!/bin/bash

# Directory containing the SLURM log files
root="$(git rev-parse --show-toplevel)"
log_dir="${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/"

# Timeout message to search for
timeout_message="DUE TO TIME LIMIT"

# Loop through each log file in the directory
for logfile in "$log_dir"/phylo_run_iqtree_loci_*.log; do
    if grep -q "$timeout_message" "$logfile"; then
        # Extract the array job ID from the filename
        array_id=$(echo "$logfile" | sed -E 's/.*_loci_([0-9]+)\.log/\1/')
        echo "Job with array ID $array_id timed out."
    fi
done
