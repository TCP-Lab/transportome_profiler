FROM archlinux:base-devel

## The arch maintainers suggest to do this before anything else
RUN pacman --noconfirm -Syu

## Install R and R packages
RUN pacman --noconfirm -Syu r

COPY ./src/helper_scripts/install_r_pkgs.R .
RUN Rscript --vanilla ./install_r_pkgs.R

## Install Rust
RUN pacman --noconfirm -Syu rust

COPY ./src/requirements.txt /src/
RUN pacman --noconfirm -Syu ttf-fira-code

# Install python and python packages
RUN pacman --noconfirm -Syu python && python -m pip install ./src/requirements.txt

# Install rust dependencies
RUN cargo install xsv
RUN cargo install --git https://github.com/MrHedmad/fast-cohen.git

# Install kerblam
RUN cargo install --locked --git https://github.com/MrHedmad/kerblam.git

# Copy the rest of the files
COPY . .
