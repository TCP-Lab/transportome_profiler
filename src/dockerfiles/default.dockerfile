FROM archlinux:base-devel-20250420.0.338771

## The arch maintainers suggest to do this before anything else
RUN pacman --noconfirm -Syu

## Install R and R packages
RUN pacman --noconfirm -Syu r

COPY ./src/helper_scripts/install_r_pkgs.R .
RUN Rscript --vanilla ./install_r_pkgs.R

## Install Rust
RUN pacman --noconfirm -Syu rust

## Set a timezone file (for the lubridate R package)
RUN ln -sf /usr/share/zoneinfo/Etc/GMT /etc/localtime

# Install python and python packages
COPY ./src/requirements.txt /src/
RUN pacman --noconfirm -Syu git python python-pip && python -m pip install --break-system-packages -r ./src/requirements.txt

# Install rust dependencies
RUN cargo install xsv
RUN cargo install --git https://github.com/MrHedmad/fast-cohen.git

# Install kerblam
RUN cargo install --locked --git https://github.com/MrHedmad/kerblam.git

# Install miscellaneous other packages
RUN pacman --noconfirm -Syu jq ttf-fira-code ripgrep

ENV PATH="$PATH:/root/.cargo/bin"

# Copy the rest of the files
COPY . .
