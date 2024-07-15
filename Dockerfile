FROM ubuntu:latest

FROM python:3.9

RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    git 

RUN wget https://github.com/lh3/minimap2/releases/download/v2.27/minimap2-2.27.tar.bz2 && \
    tar -xvjf minimap2-2.27.tar.bz2 && \
    cd minimap2-2.27 && \
    make && \
    cp minimap2 /usr/bin/

# download and install samtools version 1.19.2
RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 && \
    tar -xvjf samtools-1.19.2.tar.bz2 && \
    cd samtools-1.19.2 && \
    ./configure && \
    make && \
	  make install

# install miniconda (adjust the version as needed)
ENV MINICONDA_VERSION py39_24.5.0-0
ENV CONDA_DIR /miniconda3

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-$MINICONDA_VERSION-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

# make non-activate conda commands available
ENV PATH $CONDA_DIR/bin:$PATH
    
# make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile
    
# make conda activate command available from /bin/bash --interative shells
RUN conda init bash

RUN conda install -c "bioconda/label/cf201901" cd-hit
RUN conda install -c bioconda gffread
RUN conda install -c bioconda gffcompare
RUN conda install -c forge metamorpheus

RUN pip install biopython==1.77
RUN pip install py-cdhit
RUN pip install peptides

FROM rocker/r-ver:4.3.2

FROM bioconductor/bioconductor_docker:devel

# install R packages required 
RUN R -e 'install.packages(c(\
              "shiny", \
              "shinyjs", \
              "shinythemes", \
              "shinydashboard", \
              "data.table", \ 
              "dplyr", \
              "magrittr", \
              "ggplot2", \
              "tidyr", \
              "readr", \
              "tibble", \
              "purrr", \
              "stringr", \
              "forcats", \
              "markdown", \
              "devtools", \
              "ggrepel", \
              "phylotools", \
              "optparse"), \
              repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23")'

# Bioconductor
RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("rtracklayer", \
              "ORFik", \
              "GenomicFeatures", \
              "GenomicRanges", \
              "GenomicAlignments",\
              "BSgenome.Mmusculus.UCSC.mm39",\
              "BSgenome.Hsapiens.UCSC.hg38",\
              "bambu", \
              "Biostrings", \
              "SummarizedExperiment", \
              "Rsamtools", \
              "vsn", \
              "patchwork", ask = F))'

# install ggtranscript from github with devtools
RUN R -e "devtools::install_github('dzhang32/ggtranscript')"

# copy the shiny app directory into the image
COPY ./GenomeProt /srv/shiny-server/

# expose local port
EXPOSE 3838

# run shiny app 
CMD Rscript -e "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=3838)"