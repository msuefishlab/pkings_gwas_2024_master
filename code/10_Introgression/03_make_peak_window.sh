#!/bin/bash
set -euo pipefail

# Usage: ./make_window_bed.sh peaks.bed 100 peaks_100kb.bed
#   $1 = input BED file
#   $2 = window size (in kb)
#   $3 = output BED file

if [ $# -ne 3 ]; then
  echo "Usage: $0 <input.bed> <window_size_kb> <output.bed>"
  exit 1
fi

input="$1"
window_kb="$2"
output="$3"

# Convert kb to bp
half_window=$(( (window_kb * 1000) / 2 ))

awk -v half="$half_window" '{
  mid = int(($2 + $3) / 2);
  new_start = (mid - half < 0 ? 0 : mid - half);
  new_end   = mid + half;
  print $1, new_start, new_end, $4
}' OFS="\t" "$input" > "$output"

