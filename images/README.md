# PKINGS GWAS Images

	singularity pull docker://broadinstitute/gatk:4.3.0.0
	singularity pull https://depot.galaxyproject.org/singularity/plink2:2.00a3.3--hb2a7ceb_0
	singularity pull https://depot.galaxyproject.org/singularity/gemma:0.98--h9dd4a16_0

## GWAS Tools

	docker build --platform linux/amd64 -t gwas_tools .
	docker tag gwas_tools jasongallant/gwas_tools
	docker push jasongallant/gwas_tools

## Convert to Singularity
	cd images/
	singularity build gwas_tools.sif docker://jasongallant/gwas_tools:latest