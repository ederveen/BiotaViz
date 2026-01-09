FROM ubuntu:20.04

ARG USER=docker

# set environment without graphical interface
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    autoconf automake \
    make \
    pandoc \
    cargo \
    curl wget pgp \
    libtool \
    libcurl4-openssl-dev \ 
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    libncurses5-dev \
    libgdbm-dev \
    libnss3-dev \
    libreadline-dev \
    libffi-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

#------------------------------------------------------------------------------------#
# 2. Python:3.11 setup
#------------------------------------------------------------------------------------#
# Getting python 3.11 from source
RUN wget https://www.python.org/ftp/python/3.11.3/Python-3.11.3.tgz \
    && tar -xvf Python-3.11.3.tgz \
    && cd Python-3.11.3 \
    && ./configure --enable-optimizations \
    && make altinstall \
    && rm -rf /tmp/* /var/lib/apt/lists/*

# Install pip3 for Python 3.11
RUN wget https://bootstrap.pypa.io/get-pip.py \
    && python3.11 get-pip.py

RUN pip install \
    numpy \
    biom-format \
    h5py

#------------------------------------------------------------------------------------#
# 3. R-base:latest setup
#------------------------------------------------------------------------------------#
# Setting up R-base:latest, includes public key setup and addition of r-base version to apt depository
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" >> /etc/apt/sources.list
# apt update, upgrade and clean after above set-up is required!
# Otherwise it gives an error of broken pipes.
RUN apt-get update && apt-get upgrade -y && apt-get clean
RUN apt-get update && apt-get install -y r-base

# Copy requirements
COPY install2.r .

# Required package for install2.r
RUN R -e "install.packages('docopt', dependencies=TRUE)"

# Essential R dependencies
RUN Rscript install2.r --error --skipinstalled \
    networkD3 \
    tibble \
    htmlwidgets \
    rmarkdown \
    && rm -rf /tmp/downloaded_packages

#------------------------------------------------------------------------------------#
# 3. Adding executables
#------------------------------------------------------------------------------------#
# Python executables
COPY biom2biotaviz.py /usr/local/bin/
COPY clean_biom_txt.py /usr/local/bin/
COPY Biotaviz_counts_to_abundance.py /usr/local/bin/
COPY sankey-file-prep.py /usr/local/bin/

RUN chmod +x /usr/local/bin/*.py

# R executables
COPY sankey-diagram-html-generator.R /usr/local/bin/
RUN echo '#!/bin/bash' > /usr/local/bin/sankey-diagram-html-generator.R \
    && echo 'exec Rscript /usr/local/bin/sankey-diagram-html-generator.R "$@"' >> /usr/local/bin/sankey-diagram-html-generator \
    && chmod +x /usr/local/bin/sankey-diagram-html-generator \
    && chmod +x /usr/local/bin/sankey-diagram-html-generator.R

COPY sankey-diagram-png-generator.R /usr/local/bin/
RUN echo '#!/bin/bash' > /usr/local/bin/sankey-diagram-png-generator.R \
    && echo 'exec Rscript /usr/local/bin/sankey-diagram-png-generator.R "$@"' >> /usr/local/bin/sankey-diagram-png-generator \
    && chmod +x /usr/local/bin/sankey-diagram-png-generator \
    && chmod +x /usr/local/bin/sankey-diagram-png-generator.R

#------------------------------------------------------------------------------------#
# 4. non-root user
#------------------------------------------------------------------------------------#

ARG USERNAME=docker
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME