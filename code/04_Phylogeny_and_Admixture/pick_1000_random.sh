
mkdir -p /mnt/research/efish/pkings_gwas_2024_master/output_data/04_Phylogeny_And_Admixture/1000_rand_filtered_gene_alignments
cd /mnt/research/efish/pkings_gwas_2024_master/output_data/04_Phylogeny_And_Admixture/filtered_gene_alignments
find . -type f -print0 | shuf -z -n 1000 | xargs -0 -I {} cp {} /mnt/research/efish/pkings_gwas_2024_master/output_data/04_Phylogeny_And_Admixture/1000_rand_filtered_gene_alignments

