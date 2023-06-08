FROM rocker/tidyverse:4.3.0

# Install prerequisites
RUN apt-get update && apt-get install --yes curl software-properties-common

# Copy the required R packages
COPY ./src/helper_scripts/install_r_pkgs.R .
# Install R packages - the tidyverse is already compiled
RUN Rscript ./install_r_pkgs.R

# Install python 3.11
RUN add-apt-repository --yes ppa:deadsnakes/ppa && apt-get update && apt-get install --yes python3.11 python3.11-venv

# Install tree
RUN apt install tree

# Install xsv
RUN XSV_VERSION=$(curl -s "https://api.github.com/repos/BurntSushi/xsv/releases/latest" | grep -Po '"tag_name": "\K[0-9.]+') ; curl -Lo xsv.tar.gz "https://github.com/BurntSushi/xsv/releases/latest/download/xsv-${XSV_VERSION}-x86_64-unknown-linux-musl.tar.gz" ; tar xf xsv.tar.gz -C /usr/local/bin

# Copy the requirements for the python env + the makefile
COPY requirements.txt makefile ./
# Make the python env
RUN make env

# Copy the rest of the files
COPY . .

# Run the analysis
CMD make
