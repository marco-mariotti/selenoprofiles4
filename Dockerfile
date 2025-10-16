# Use a Conda base image
FROM continuumio/miniconda3

# Install dependencies
RUN conda install -y python=3.12.0 -c bioconda -c conda-forge

RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    zlib1g-dev \
    gcc \
    g++ \
    libffi-dev \
    libssl-dev \
    && apt-get clean

RUN pip install selenoprofiles4[addons]

ENV PATH=$PATH:/blast-2.2.26/bin

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz && \
    tar -xzvf blast-2.2.26-x64-linux.tar.gz    

# Run selenoprofiles setup and download commands
RUN selenoprofiles -setup && \
    yes "" | selenoprofiles -download 

# Accept Anaconda Terms of Service
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install bioconda::exonerate

RUN conda install -c anaconda gawk

RUN conda install bioconda::wise2

RUN conda install -c mmariotti -c conda-forge -c etetoolkit ncbi_db
