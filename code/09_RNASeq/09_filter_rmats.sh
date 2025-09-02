source pkings_gwas.env

module purge; module load Biopython/1.83-foss-2023a

sc_rmats_filter=${root}/code/09_RNASeq/rmats_filtering.py
declare -a event_array=("SE" "MXE" "RI" "A5SS" "A3SS")
declare -a counttype_array=("JC" "JCEC")

mkdir -p ${root}/output_data/09_RNASeq/rmats_filtered/gw1/
cd ${root}/output_data/09_RNASeq/rmats_filtered/gw1/

for event in "${event_array[@]}"; do
    for counttype in "${counttype_array[@]}"; do
        python $sc_rmats_filter ${root}/output_data/09_RNASeq/rmats/gw1/post/${event}.MATS.${counttype}.txt 10,0.05,0.95,0.01,0.05,0.5,0.2
    done
done