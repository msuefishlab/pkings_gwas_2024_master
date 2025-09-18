# Association Analysis Pipeline

Here the pipeline forks into several comparisons between TP and BP P. kingsleyae groups.  Based on the phylogenetic tree, BP originated multiple times within P. kingsleaye...

Mengono Creek (Cocobeach)
Biroundou Creek
Bikagala + Bavavela Creek
Bambomo Creek

We don't have enough samples from Mengono Creek (in this analysis... perhaps the poolseq data could be repurposed?), and I want to compare to the nearest TP population to minimize effects of genetic drift.  These are the comparisions we want to do:

APA+BENGUE to BAM (no wobbles)
DOV/DOUENGI to BIK/BAV (no wobbles)
DOV/DOUENGI to BIR (no wobbles)

We want to demonstrate the fundamental problem of "no peaks" by comparing all P. kingsleyae.

Finally, what are wobbles?  Let's focus on APA/BAM, and score TP=1, WOB=2, BP=3


## Create Keep Files and Phenotype Files for Each Comparision

	source pkings_gwas.env

	kfiledir=${root}/input_data/06_Association/

	mkdir -p ${kfiledir}

	# All P. kingsleyae, no wobbles
	singularity exec --bind $root:/project_root --bind $kfiledir:/out_dir ${gwas_tools_image} python3 /project_root/code/06_Association/generate_phenotype_and_keepfiles.py \
	/project_root/input_data/01_Terra/data_model/fish_data_2023.txt \
	--output_directory /out_dir \
	--include_species "Paramormyrops kingsleyae"
	# Enter TP=1 BP=2, WOB=-9

	# Ipassa, Bengue, Bambomo, no wobbles
	singularity exec --bind $root:/project_root --bind $kfiledir:/out_dir ${gwas_tools_image} python3 /project_root/code/06_Association/generate_phenotype_and_keepfiles.py \
	/project_root/input_data/01_Terra/data_model/fish_data_2023.txt \
	--output_directory /out_dir \
	--include_pop APA BEN BAM
	# Enter TP=1 BP=2, WOB=-9

	#Douvalou, Douengi, Bikagala, Bavavela, BIR
	singularity exec --bind $root:/project_root --bind $kfiledir:/out_dir ${gwas_tools_image} python3 /project_root/code/06_Association/generate_phenotype_and_keepfiles.py \
	/project_root/input_data/01_Terra/data_model/fish_data_2023.txt \
	--output_directory /out_dir \
	--include_pop DOV DOG BIK BAVA BIR
	# Enter TP=1 BP=2, WOB=-9

	# Ipassa, Bengue,IVI,MOV, DOV,DOG Bambomo, no wobbles
	singularity exec --bind $root:/project_root --bind $kfiledir:/out_dir ${gwas_tools_image} python3 /project_root/code/06_Association/generate_phenotype_and_keepfiles.py \
	/project_root/input_data/01_Terra/data_model/fish_data_2023.txt \
	--output_directory /out_dir \
	--include_pop APA BEN BAM IVI,MOV,DOV,DOG
	# Enter TP=1 BP=2, WOB=-9

	#Douvalou, Douengi, Bikagala, Bavavela, BIR IVI MOV DOV DOG, no wobbles
	singularity exec --bind $root:/project_root --bind $kfiledir:/out_dir ${gwas_tools_image} python3 /project_root/code/06_Association/generate_phenotype_and_keepfiles.py \
	/project_root/input_data/01_Terra/data_model/fish_data_2023.txt \
	--output_directory /out_dir \
	--include_pop DOV DOG BIK BAVA BIR IVI MOV DOV DOG
	# Enter TP=1 BP=2, WOB=-9

	#Douvalou, Douengi, Biroundou
	singularity exec --bind $root:/project_root --bind $kfiledir:/out_dir ${gwas_tools_image} python3 /project_root/code/06_Association/generate_phenotype_and_keepfiles.py \
	/project_root/input_data/01_Terra/data_model/fish_data_2023.txt \
	--output_directory /out_dir \
	--include_pop  DOV DOG BIR
	# Enter TP=1 BP=2, WOB=-9

	#APA BAM TP, BP, WOBBLES (RETURN TO LATER)
	singularity exec --bind $root:/project_root --bind $kfiledir:/out_dir ${gwas_tools_image} python3 /project_root/code/06_Association/generate_phenotype_and_keepfiles.py \
	/project_root/input_data/01_Terra/data_model/fish_data_2023.txt \
	--output_directory /out_dir \
	--include_pop APA BAM
	# Enter TP=1 BP=3, WOB=2



*** YOU ARE HERE ****
	#[DOV,DOG,BEN,APA],[BAM], NO WOBBLES
	singularity exec --bind $root:/project_root --bind $kfiledir:/out_dir ${gwas_tools_image} python3 /project_root/code/06_Association/generate_phenotype_and_keepfiles.py \
	/project_root/input_data/01_Terra/data_model/fish_data_2023.txt \
	--output_directory /out_dir \
	--include_pop DOV DOG BEN APA BAM
	# Enter TP=1 BP=2, WOB=-9

	#[DOV,DOG,BEN,APA],[BIK,BAVA,BIR] NO WOBBLES
	singularity exec --bind $root:/project_root --bind $kfiledir:/out_dir ${gwas_tools_image} python3 /project_root/code/06_Association/generate_phenotype_and_keepfiles.py \
	/project_root/input_data/01_Terra/data_model/fish_data_2023.txt \
	--output_directory /out_dir \
	--include_pop DOV DOG BEN APA BIK BAVA BIR
	# Enter TP=1 BP=2, WOB=-9




## Prepare Data
		bash code/06_Association/01_run_prep_data.sh PKINGS_TP1_BP2_WOB9
		bash code/06_Association/01_run_prep_data.sh DOV_DOG_BIK_BAV_BIR_TP1_BP2_WOB9
		bash code/06_Association/01_run_prep_data.sh DOV_DOG_BIK_BAVA_BIR_IVI_MOV_DOV_DOG_TP1_BP2_WOB9
		bash code/06_Association/01_run_prep_data.sh APA_BEN_BAM_IVI_MOV_DOV_DOG_TP1_BP2_WOB9

		bash code/06_Association/01_run_prep_data.sh DOV_DOG_BEN_APA_BAM_TP1_BP2_WOB9
		bash code/06_Association/01_run_prep_data.sh DOV_DOG_BEN_APA_BIK_BAVA_BIR_TP1_BP2_WOB9


## Make RelMat
	bash code/06_Association/02_run_rel_mat.sh APA_BEN_BAM_TP1_BP2_WOB9
	bash code/06_Association/02_run_rel_mat.sh DOV_DOG_BIK_BAV_BIR_TP1_BP2_WOB9
	bash code/06_Association/02_run_rel_mat.sh DOV_DOG_BIR_TP1_BP2_WOB9
	bash code/06_Association/02_run_rel_mat.sh DOV_DOG_BIK_BAVA_BIR_IVI_MOV_DOV_DOG_TP1_BP2_WOB9
	bash code/06_Association/02_run_rel_mat.sh APA_BEN_BAM_IVI_MOV_DOV_DOG_TP1_BP2_WOB9

	bash code/06_Association/02_run_rel_mat.sh DOV_DOG_BEN_APA_BAM_TP1_BP2_WOB9
	bash code/06_Association/02_run_rel_mat.sh DOV_DOG_BEN_APA_BIK_BAVA_BIR_TP1_BP2_WOB9


## Run Assoc
	bash code/06_Association/03_run_gemma.sh APA_BEN_BAM_TP1_BP2_WOB9
	bash code/06_Association/03_run_gemma.sh DOV_DOG_BIK_BAV_BIR_TP1_BP2_WOB9
	bash code/06_Association/03_run_gemma.sh DOV_DOG_BIR_TP1_BP2_WOB9
	bash code/06_Association/03_run_gemma.sh DOV_DOG_BIK_BAVA_BIR_IVI_MOV_DOV_DOG_TP1_BP2_WOB9
	bash code/06_Association/03_run_gemma.sh APA_BEN_BAM_IVI_MOV_DOV_DOG_TP1_BP2_WOB9

	bash code/06_Association/03_run_gemma.sh DOV_DOG_BEN_APA_BAM_TP1_BP2_WOB9
	bash code/06_Association/03_run_gemma.sh DOV_DOG_BEN_APA_BIK_BAVA_BIR_TP1_BP2_WOB9


## Run Permution
	. code/06_Association/04_gemma_permutation_local.sh APA_BEN_BAM_TP1_BP2_WOB9
	. code/06_Association/04_gemma_permutation_local.sh DOV_DOG_BIK_BAV_BIR_TP1_BP2_WOB9
	. code/06_Association/04_gemma_permutation_local.sh DOV_DOG_BIR_TP1_BP2_WOB9
	. code/06_Association/04_gemma_permutation_local.sh DOV_DOG_BIK_BAVA_BIR_IVI_MOV_DOV_DOG_TP1_BP2_WOB9
	. code/06_Association/04_gemma_permutation_local.sh APA_BEN_BAM_IVI_MOV_DOV_DOG_TP1_BP2_WOB9

	. code/06_Association/04_gemma_permutation_local.sh DOV_DOG_BEN_APA_BAM_TP1_BP2_WOB9
	. code/06_Association/04_gemma_permutation_local.sh DOV_DOG_BEN_APA_BIK_BAVA_BIR_TP1_BP2_WOB9


## Collect Min Values from Permution
	. code/06_Association/05_collect_min_values.sh APA_BEN_BAM_TP1_BP2_WOB9
	. code/06_Association/05_collect_min_values.sh DOV_DOG_BIK_BAV_BIR_TP1_BP2_WOB9
	. code/06_Association/05_collect_min_values.sh DOV_DOG_BIR_TP1_BP2_WOB9
	. code/06_Association/05_collect_min_values.sh DOV_DOG_BIK_BAVA_BIR_IVI_MOV_DOV_DOG_TP1_BP2_WOB9
	. code/06_Association/05_collect_min_values.sh APA_BEN_BAM_IVI_MOV_DOV_DOG_TP1_BP2_WOB9

	. code/06_Association/05_collect_min_values.sh DOV_DOG_BEN_APA_BAM_TP1_BP2_WOB9
	. code/06_Association/05_collect_min_values.sh DOV_DOG_BEN_APA_BIK_BAVA_BIR_TP1_BP2_WOB9


## Render First GWAS Report
	bash code/06_Association/06_Render_GWAS_report.sh APA_BEN_BAM_TP1_BP2_WOB9
	bash code/06_Association/06_Render_GWAS_report.sh DOV_DOG_BIK_BAV_BIR_TP1_BP2_WOB9
	bash code/06_Association/06_Render_GWAS_report.sh DOV_DOG_BIR_TP1_BP2_WOB9
	bash code/06_Association/06_Render_GWAS_report.sh DOV_DOG_BIK_BAVA_BIR_IVI_MOV_DOV_DOG_TP1_BP2_WOB9
	bash code/06_Association/06_Render_GWAS_report.sh APA_BEN_BAM_IVI_MOV_DOV_DOG_TP1_BP2_WOB9

	bash code/06_Association/06_Render_GWAS_report.sh DOV_DOG_BEN_APA_BAM_TP1_BP2_WOB9
	bash code/06_Association/06_Render_GWAS_report.sh DOV_DOG_BEN_APA_BIK_BAVA_BIR_TP1_BP2_WOB9


## Get LD SNPs For Round 2 GWAS Report
	bash code/06_Association/07_get_ld_snps.sh APA_BEN_BAM_TP1_BP2_WOB9
	bash code/06_Association/07_get_ld_snps.sh DOV_DOG_BIK_BAV_BIR_TP1_BP2_WOB9
	bash code/06_Association/07_get_ld_snps.sh DOV_DOG_BIR_TP1_BP2_WOB9
	bash code/06_Association/07_get_ld_snps.sh DOV_DOG_BIK_BAVA_BIR_IVI_MOV_DOV_DOG_TP1_BP2_WOB9
	bash code/06_Association/07_get_ld_snps.sh APA_BEN_BAM_IVI_MOV_DOV_DOG_TP1_BP2_WOB9

	bash code/06_Association/07_get_ld_snps.sh DOV_DOG_BEN_APA_BAM_TP1_BP2_WOB9
	bash code/06_Association/07_get_ld_snps.sh DOV_DOG_BEN_APA_BIK_BAVA_BIR_TP1_BP2_WOB9
