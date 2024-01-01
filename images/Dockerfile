FROM bioconductor/bioconductor_docker:3.17

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
	phangorn

RUN R -e 'BiocManager::install("gdsfmt")'
RUN R -e 'BiocManager::install("SNPRelate")'

RUN mkdir -p /snphylo/bin

RUN apt-get update && \
    apt-get install -y muscle phylip curl git python2 bcftools vcftools tabix

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

RUN python3 setup.py install





CMD ["/bin/bash"]
