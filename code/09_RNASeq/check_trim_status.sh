## check_trim_status.sh
## 09-30-18
## JRG
## Checks Status of Trimmomatic Jobs

#!/bin/bash


printf 'file\tstatus\tboth_surviving\tboth_surviving_percent\tforward_only\tforward_only_percent\treverse_only\treverese_only_percent\tdropped\tdropped_percent\n'
for f in $SCRATCH/pkings_trimmed/*/*results.log;
  do
    status=$(if grep -q "TrimmomaticPE: Completed successfully" $f; then echo "complete"; else echo "incomplete"; fi);
    trimdata=$(grep "Input Read Pairs:" $f || echo "";)
    parsed_trimdata=$(echo -e "$trimdata"| grep -Eo -e "[0-9]+\.[0-9]+" -e "[0-9]+"| tr '\n' '\t')
    line=$(printf '%b\t%b\t%b\n' $f $status $parsed_trimdata)
    echo -e $line
  done
