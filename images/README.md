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

## RNASeq Tools

    docker build --platform linux/amd64 -f RNAseqDockerfile -t rnaseq_tools .
    docker tag rnaseq_tools jasongallant/rnaseq_tools
    docker push jasongallant/rnaseq_tools

## Convert to Singularity

    cd images/
    singularity build rnaseq_tools.sif docker://jasongallant/rnaseq_tools:latest

## RMATs

    docker build --platform linux/amd64 -t rmats_turbo -f ./RMATS_Dockerfile .
    docker tag rmats_turbo jasongallant/rmats_turbo
    docker push jasongallant/rmats_turbo

    cd images/
    singularity build rmats_turbo.sif docker://jasongallant/rmats_turbo:latest

## Liftoff

    docker build --platform linux/amd64 -t liftoff -f ./LIFTOFF_Dockerfile .
    docker tag liftoff jasongallant/liftoff
    docker push jasongallant/liftoff

    cd images/
    singularity build liftoff.sif docker://jasongallant/liftoff:latest

## TWISST

    docker build --platform linux/amd64 -f TwisstDockerFile -t twisst_tools .
    docker tag twisst_tools jasongallant/twisst_tools
    docker push jasongallant/twisst_tools

    cd images
    singularity build twisst_tools.sif docker://jasongallant/twisst_tools:latest

## DSuite

    docker build --platform linux/amd64 -f DsuiteDocker -t dsuite_tools .
    docker tag dsuite_tools jasongallant/dsuite_tools
    docker push jasongallant/dsuite_tools

    cd images
    singularity build dsuite_tools.sif docker://jasongallant/dsuite_tools:latest

## Dadi

    docker build --platform linux/amd64 -f DadiDocker -t dadi_tools .
    docker tag dadi_tools jasongallant/dadi_tools
    docker push jasongallant/dadi_tools

    cd images
    singularity build dadi_tools.sif docker://jasongallant/dadi_tools:latest
