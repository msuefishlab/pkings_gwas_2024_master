source pkings_gwas.env
sc_rmats_filter=${root}/code/09_RNAseq/rmats_filtering.py
declare -a event_array=("SE" "MXE" "RI" "A5SS" "A3SS")
declare -a counttype_array=("JC" "JCEC")

mkdir -p ${root}/output_data/09_RNASeq/rmats_filtered/gw2/
cd ${root}/output_data/09_RNASeq/rmats_filtered/gw2/

for event in "${event_array[@]}"; do
    for counttype in "${counttype_array[@]}"; do
        python $sc_rmats_filter ../../rmats/gw2/post/${event}.MATS.${counttype}.txt 10,0.05,0.95,0.01,0.05,0.5,0.2
    done
done