FROM ubuntu:24.04

LABEL maintainer="lukas.schmidbauer@oth-regensburg.de"

WORKDIR /repro

RUN apt update \
    && DEBIAN_FRONTEND=noninteractive apt install -y build-essential wget cmake gfortran gobjc gobjc++ gnustep gnustep-devel libbz2-dev liblzma-dev libpcre2-dev libcurl4-openssl-dev libcairo2-dev libtiff5-dev libreadline-dev libxml2-dev libharfbuzz-dev libfribidi-dev libglpk-dev libgsl-dev libgmp-dev libmpc-dev libudunits2-dev libgdal-dev libmagick++-dev \
    && apt clean \
    && rm -rf /var/lib/apt/lists/*

# ENV CONDA_DIR /opt/conda
# RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh -O ~/miniconda.sh && \
#     /bin/bash ~/miniconda.sh -b -p /opt/conda

# ENV PATH=$CONDA_DIR/bin:$PATH

RUN apt-get clean && apt-get update
RUN apt-get install -y texlive-full

RUN apt clean && apt update
RUN apt install -y software-properties-common dirmngr
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN DEBIAN_FRONTEND=noninteractive apt install -y r-base r-base-dev

RUN R -e "install.packages('tidyverse', dependencies=TRUE)"
RUN R -e "install.packages('patchwork', dependencies=TRUE)"
RUN R -e "install.packages('tikzDevice', dependencies=TRUE)"
RUN R -e "install.packages('scales', dependencies=TRUE)"
RUN R -e "install.packages('ggh4x', dependencies=TRUE)"
RUN R -e "install.packages('ggpmisc', dependencies=TRUE)"
RUN R -e "install.packages('stringr', dependencies=TRUE)"
RUN R -e "install.packages('ggpubr', dependencies=TRUE)"
RUN R -e "install.packages('gridExtra', dependencies=TRUE)"
RUN R -e "install.packages('rlist', dependencies=TRUE)"
RUN R -e "install.packages('Cairo', dependencies=TRUE)"

RUN apt-get -y install python3-pip python3.12-venv 
RUN apt install 7zip

ENV VENV=/venv
RUN python3 -m venv $VENV
ENV PATH="$VENV/bin:$PATH"

RUN pip install pandas
RUN pip install numpy
RUN pip install dill
RUN pip install matplotlib
RUN pip install dwave-ocean-sdk==7.1.0

RUN pip install networkx

RUN pip install qiskit==1.1.0
RUN pip install qiskit-aer==0.14.2
RUN pip install qiskit-ibm-runtime==0.24.1

COPY ./ ./

#OpenShell
CMD ["/bin/bash"]
