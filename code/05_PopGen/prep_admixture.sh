#!/bin/bash --login

########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=24:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=16           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=30G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name PREP_ADMIXTURE      # you can give your job a name for easier identification (same as -J)
#SBATCH -A data-machine
########## Command Lines to Run ##########
cd $root

source $root/pkings_gwas.env

MYNAME=$(whoami)

indir=${root}/output_data/03_QC/
outdir=${root}/output_data/05_PopGen/
tmpoutdir=/mnt/local/${MYNAME}/${SLURM_JOB_ID}/

mkdir -p ${outdir}
mkdir -p ${tmpoutdir}

OUTROOT=${snps_only_vcf%".vcf.gz"}


#prefilter VCF for MAF < 0.1 and Missingness > 10% and DP < 5
cat ${root}/input_data/00_Reference_Genome/${chr_name_map} | cut -f2,3 > ${root}/input_data/00_Reference_Genome/${chr_name_map}.chromosome_number_map.txt

echo singularity exec --bind $root:/project_root --bind $outdir:/out_dir --bind $indir:/in_dir --bind $tmpoutdir:/tmp_dir ${gwas_tools_image} bash -c "bcftools view --threads 16 /in_dir/${snps_only_vcf} | bcftools annotate --rename-chrs /project_root/input_data/00_Reference_Genome/${chr_name_map}.chromosome_number_map.txt  | bcftools filter -e \"F_MISSING > 0.20 || MAF < 0.05 || INFO/DP < 5\" -o /tmp_dir/${OUTROOT}.renamed.maf5.miss20.dp5.vcf"

singularity exec --bind $root:/project_root --bind $outdir:/out_dir --bind $indir:/in_dir --bind $tmpoutdir:/tmp_dir ${gwas_tools_image} bash -c "bcftools view --threads 16 /in_dir/${snps_only_vcf} | bcftools annotate --rename-chrs /project_root/input_data/00_Reference_Genome/${chr_name_map}.chromosome_number_map.txt  | bcftools filter -e \"F_MISSING > 0.20 || MAF < 0.05 || INFO/DP < 5\" -o /tmp_dir/${OUTROOT}.renamed.maf5.miss20.dp5.vcf"

echo singularity exec --bind $root:/project_root --bind $outdir:/out_dir --bind $indir:/in_dir --bind $tmpoutdir:/tmp_dir ${gwas_tools_image} bash -c "bcftools +prune -n 1 -w 20000bp --nsites-per-win-mode rand --random-seed 42 /tmp_dir/${OUTROOT}.renamed.maf5.miss20.dp5.vcf -o /tmp_dir/${OUTROOT}.renamed.maf5.miss20.dp5.r20kb.vcf"

singularity exec --bind $root:/project_root --bind $outdir:/out_dir --bind $indir:/in_dir --bind $tmpoutdir:/tmp_dir ${gwas_tools_image} bash -c "bcftools +prune -n 1 -w 20000bp --nsites-per-win-mode rand --random-seed 42 /tmp_dir/${OUTROOT}.renamed.maf5.miss20.dp5.vcf -o /tmp_dir/${OUTROOT}.renamed.maf5.miss20.dp5.r20kb.vcf"

singularity exec --bind $root:/project_root --bind $outdir:/out_dir --bind $tmpoutdir:/tmp_dir ${gwas_tools_image} /plink/plink \
--vcf  /tmp_dir/${OUTROOT}.renamed.maf5.miss20.dp5.r20kb.vcf \
--allow-extra-chr \
--allow-no-sex \
--double-id \
--make-bed \
--pheno /project_root/input_data/06_Association/${OUTNAME}.pheno.txt \
--out /out_dir/$OUTROOT.renamed.maf5.miss20.dp5.r20kb.for_admixture