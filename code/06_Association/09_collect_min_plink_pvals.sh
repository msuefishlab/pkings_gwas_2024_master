#!/usr/bin/env bash
set -euo pipefail

# Collect minimum p-values per TEST from PLINK .model permutation outputs
# Usage:
#   . code/06_Association/04_collect_min_values_plink_model.sh APA_BAM_BEN_TPBP0_WOB1
#
# Output files (inside the OUTNAME dir):
#   - ${OUTNAME}_perm_min_by_test.long.tsv  (perm, TEST, minP)
#   - ${OUTNAME}_perm_min_by_test.wide.tsv  (perm, TEST columns as headers)

root="$(git rev-parse --show-toplevel)"
OUTNAME="$1"

cd "$root"
source "$root/pkings_gwas.env"

basedir="${root}/output_data/06_Association/${OUTNAME}"
# Pattern like: ${OUTNAME}_PLINK_perm_94.model
pattern="${basedir}/${OUTNAME}_PLINK_perm_"'*.model'

long_out="${basedir}/${OUTNAME}_perm_min_by_test.long.tsv"
wide_out="${basedir}/${OUTNAME}_perm_min_by_test.wide.tsv"

# Clean outputs
rm -f "$long_out" "$wide_out"
mkdir -p "$basedir"

# Header for long (tidy) output
echo -e "perm\tTEST\tminP" > "$long_out"

shopt -s nullglob
files=( $pattern )
if [[ ${#files[@]} -eq 0 ]]; then
  echo "No files matched: $pattern" >&2
  exit 1
fi

# Loop through permutation files
for f in "${files[@]}"; do
  # extract permutation id from filename
  # expects ..._perm_<N>.model
  permid="$(basename "$f" | sed -E 's/.*_perm_([0-9]+)\.model/\1/')"

  # For this file, compute min p per TEST (ignoring NA/nan/-/blank)
  awk -v perm="$permid" '
    BEGIN {
      OFS = "\t"
    }
    NR==1 { 
      # header: CHR SNP A1 A2 TEST AFF UNAFF P
      # skip header row
      next 
    }
    {
      test  = $5
      p_raw = $8
      # sanitize: skip NA/nan/-/empty
      if (p_raw ~ /^(NA|NaN|nan|-)?$/) next
      # convert to numeric
      p = p_raw + 0
      if (!(test in minP) || p < minP[test]) {
        minP[test] = p
      }
      tests_seen[test] = 1
    }
    END {
      for (t in tests_seen) {
        if (t in minP) {
          printf "%s\t%s\t%.12g\n", perm, t, minP[t]
        } else {
          # if somehow seen but no valid p, write NA
          printf "%s\t%s\tNA\n", perm, t
        }
      }
    }
  ' "$f" >> "$long_out"
done

# Make a wide table: one row per perm, columns per TEST
# 1) discover all TEST names present, sorted naturally
tests_list="$(cut -f2 "$long_out" | tail -n +2 | sort -u)"

# 2) print header
{
  printf "perm"
  while IFS= read -r t; do
    printf "\t%s", "$t"
  done <<< "$tests_list"
  printf "\n"
} > "$wide_out"

# 3) pivot from long → wide
# sort long by perm (numeric) then TEST to get deterministic output
sort -k1,1n -k2,2 "$long_out" | awk -v tests="$tests_list" '
  BEGIN {
    OFS="\t"
    # split tests into an array to preserve column order
    n=split(tests, T, "\n")
  }
  NR==1 { next }  # skip header in long
  {
    perm=$1; test=$2; p=$3
    key=perm SUBSEP test
    val[key]=p
    perms[perm]=1
  }
  END {
    # For each perm, print in the header test order
    # (Use NA if a test is missing)
    # Gather and sort perms numerically
    nperm=0
    for (p in perms) {
      list[++nperm]=p
    }
    # numeric sort
    asort(list)
    for (i=1; i<=nperm; i++) {
      p = list[i]
      printf "%s", p
      for (j=1; j<=n; j++) {
        t = T[j]
        k = p SUBSEP t
        v = (k in val ? val[k] : "NA")
        printf OFS "%s", v
      }
      printf "\n"
    }
  }
' >> "$wide_out"

echo "Wrote:"
echo "  - $long_out"
echo "  - $wide_out"
