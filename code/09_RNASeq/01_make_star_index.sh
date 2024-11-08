#!/bin/bash --login

########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=04:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=17           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=50G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name STAR_CREATE_INDEX      # you can give your job a name for easier identification (same as -J)
#SBATCH --output=slurm-%x.%j.out
########## Command Lines to Run ##########

cd /mnt/research/efish/pkings_gwas_slurm/data/reference

module purge; module load GCC/8.3.0  OpenMPI/3.1.4 RSEM/1.3.3

rsem-prepare-reference -gff3 ../../data/annotation_files/ragoo.pk_to_sfor1.1.gff \
../../data/reference/ragoo.pk_0.2_to_sfor1.1.fasta \
../../data/reference/star_index/ \
--star --star-sjdboverhang 149 -p 16


