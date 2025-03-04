FROM bioconductor/bioconductor_docker:3.18

# Install R packages
RUN install2.r --error \
	tidyverse \
	ggthemes \
	patchwork \
	getopt \
	yaml \
	tools \
	data.table \
	reshape2 \
	stringr \
	ape \
	quadprog \
	igraph \
	fastmatch \
	phangorn \
	zoo \
	kableExtra \
	gggenes \
	ggrepel

RUN R -e 'BiocManager::install("gdsfmt")'
RUN R -e 'BiocManager::install("GenomicRanges")'
RUN R -e 'BiocManager::install("rtracklayer")'
RUN R -e 'BiocManager::install("SNPRelate")'
RUN R -e 'BiocManager::install("ggtree")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("biomaRt")'
RUN R -e 'BiocManager::install("BioMartGOGeneSets")'

RUN mkdir -p /snphylo/bin

RUN apt-get update && \
	apt-get install -y muscle phylip curl git python2 bcftools vcftools samtools tabix wget build-essential libgsl-dev libblas-dev liblapack-dev zlib1g-dev ruby time parallel coreutils binutils openmpi-bin openmpi-doc libopenmpi-dev 

RUN git clone https://github.com/jasongallant/SNPhylo.git

RUN mv SNPhylo /snphylo

WORKDIR /snphylo/SNPhylo

RUN ln -s $(which /usr/bin/python2.7) /usr/bin/python

COPY dnaml /usr/bin/dnaml
RUN chmod +x /usr/bin/dnaml

RUN printf "yyymusclenydnamln" | bash setup.sh

RUN chmod -R 777 /snphylo

RUN mkdir -p /pixy

WORKDIR /pixy

RUN git clone https://github.com/ksamuk/pixy.git

WORKDIR /pixy/pixy

RUN pip install .

RUN mkdir -p /seqkit

WORKDIR /seqkit

RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.6.1/seqkit_linux_amd64.tar.gz

RUN tar -zxvf seqkit_linux_amd64.tar.gz

RUN mkdir -p /qctool

WORKDIR /qctool

RUN wget https://www.well.ox.ac.uk/~gav/resources/qctool_v2.2.0-CentOS_Linux7.8.2003-x86_64.tgz

RUN tar -zxvf qctool_v2.2.0-CentOS_Linux7.8.2003-x86_64.tgz

RUN mv ./*/qctool /qctool

RUN mkdir -p /plink

WORKDIR /plink

RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip

RUN unzip plink_linux_x86_64_20231211.zip

RUN mkdir -p /admixture

WORKDIR /admixture

RUN wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz

RUN tar -zxvf admixture_linux-1.3.0.tar.gz

RUN mv dist/*/admixture .

RUN mkdir -p /gemma

WORKDIR /gemma

RUN git clone https://github.com/genetics-statistics/GEMMA.git

WORKDIR /gemma/GEMMA

RUN git checkout 8132837

RUN make

RUN mkdir -p /vcf2phylip

WORKDIR /vcf2phylip

RUN git clone https://github.com/edgardomortiz/vcf2phylip.git

RUN mkdir -p /gemma-wrapper

WORKDIR /gemma-wrapper

RUN git clone -b slurm_patched --single-branch https://github.com/jasongallant/gemma-wrapper.git

ENV GEMMA_COMMAND=/gemma/GEMMA/bin/gemma

RUN mkdir -p iqtree

WORKDIR /iqtree/

RUN wget https://github.com/iqtree/iqtree2/releases/download/v2.3.1/iqtree-2.3.1-Linux.tar.gz

RUN tar -zxvf iqtree-2.3.1-Linux.tar.gz

RUN mv iqtree*/bin/iqtree2 .

RUN mkdir -p /gff3sort

WORKDIR /gff3sort

RUN git clone https://github.com/billzt/gff3sort.git

WORKDIR /cflib

RUN git clone https://github.com/pomo-dev/cflib.git

RUN pip3 install setuptools==57.5.0
RUN pip3 install biopython pyvcf argparse pysam cflib-pomo whatshap

RUN mkdir -p /shapeit

WORKDIR /shapeit

RUN wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.17.linux.tar.gz

RUN tar -zxvf shapeit.v2.r904.glibcv2.17.linux.tar.gz

RUN mv shapeit*/bin/shapeit .

CMD ["/bin/bash"]