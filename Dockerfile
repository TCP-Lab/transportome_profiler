FROM rocker/tidyverse:4.3.0

# Install rust and cargo
RUN curl https://sh.rustup.rs -sSf | sh && rustup stable

# Install python 3.10
RUN add-apt-repository --yes ppa:deadsnakes/ppa && apt-get update && apt-get install python3.10

# Install tectonic - Ubuntu does not have a repo for it, so we download the 
# precompiled from source
RUN curl --proto '=https' --tlsv1.2 -fsSL https://drop-sh.fullyjustified.net |sh && mv ./tectonic /usr/bin/

# Install R packages - the tidyverse is already compiled
RUN Rscript \
    -e "options(repos='https://cloud.r-project.org/')" \
    -e "install.packages('BiocManager')" \
    -e "BiocManager::install('fgsea')" \
    -e "BiocManager::install('biomaRt')" \
    -e "install.packages('ggraph')" \
    -e "install.packages('assertthat')"

# Install python 3.10
RUN add-apt-repository --yes ppa:deadsnakes/ppa && apt-get update && apt-get install python3.10

# Copy the source files in. The ./data/ directory will be mounted at runtime
COPY . .

RUN make
