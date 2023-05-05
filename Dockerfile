FROM rocker/tidyverse:4.3.0

# Install prerequisites
RUN apt-get update && apt-get install --yes curl software-properties-common

# Install R packages - the tidyverse is already compiled
RUN Rscript \
    -e "options(repos='https://cloud.r-project.org/')" \
    -e "install.packages('BiocManager')" \
    -e "BiocManager::install('fgsea')" \
    -e "BiocManager::install('biomaRt')" \
    -e "install.packages('ggraph')" \
    -e "install.packages('assertthat')" \
    -e "install.packages('argparser')" \
    -e "install.packages('uuid')" \
    -e "install.packages('grid')" \
    -e "BiocManager::install('DESeq2')"

# Install rust and cargo
RUN curl https://sh.rustup.rs -sSf > rstins && sh rstins -y

ENV PATH="/root/.cargo/bin:${PATH}"

# Install python 3.10
RUN add-apt-repository --yes ppa:deadsnakes/ppa && apt-get update && apt-get install --yes python3.11 python3.11-venv

# Install tectonic - Ubuntu does not have a repo for it, so we download the 
# precompiled from source
RUN curl --proto '=https' --tlsv1.2 -fsSL https://drop-sh.fullyjustified.net |sh && mv ./tectonic /usr/bin/ && apt-get install --yes biber=2.17-2 

RUN curl http://nz2.archive.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.1f-1ubuntu2.18_amd64.deb > libss.deb && dpkg -i ./libss.deb

# Copy the source files in. The ./data/ directory will be mounted at runtime
COPY . .

# Make the python env
RUN make venv

CMD make
