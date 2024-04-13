root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/04_Phylogeny_And_Admixture/1000_rand_gene_trees

cd ${root}/output_data/0/04_Phylogeny_And_Admixture/gene_trees/

find . -type f -name "*.treefile" -print0 | shuf -z -n 1000 | xargs -0 -I {} cp {} ${root}/output_data/04_Phylogeny_And_Admixture/1000_rand_gene_trees

cat ${root}/output_data/04_Phylogeny_And_Admixture/1000_rand_gene_trees/*.treefile > ${root}/output_data/04_Phylogeny_And_Admixture/1000_rand_gene_trees.treefile