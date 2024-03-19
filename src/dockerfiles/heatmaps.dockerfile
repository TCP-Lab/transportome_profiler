FROM rocker/tidyverse:4.3.0
# Install build prerequisites
RUN apt-get update && apt-get install --yes curl software-properties-common build-essential

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y && \
    echo 'source $HOME/.cargo/env' >> $HOME/.bashrc

# Copy the required R packages
COPY ./src/helper_scripts/install_r_pkgs.R .
# Install R packages - the tidyverse is already compiled
RUN Rscript --vanilla ./install_r_pkgs.R

COPY ./src/requirements.txt /src/
# Install python 3.11 and tree
RUN add-apt-repository --yes ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install --yes python3.11 jq

# Install python packages
RUN curl --proto '=https' --tlsv1.2 -sSf https://bootstrap.pypa.io/get-pip.py | python3.11 && python3.11 -m pip install -r ./src/requirements.txt

# Make sure python is findable as 'python'
RUN ln -s "$(which python3.11)" "/usr/bin/python"

# Install xsv
RUN XSV_VERSION=$(curl -s "https://api.github.com/repos/BurntSushi/xsv/releases/latest" | grep -Po '"tag_name": "\K[0-9.]+') ; curl -Lo xsv.tar.gz "https://github.com/BurntSushi/xsv/releases/latest/download/xsv-${XSV_VERSION}-x86_64-unknown-linux-musl.tar.gz" ; tar xf xsv.tar.gz -C /usr/local/bin

# Copy the rest of the files
COPY . .
