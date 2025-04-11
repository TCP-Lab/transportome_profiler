FROM rocker/tidyverse:4.4

## Use the r2u framework for faster install of R packages
RUN apt update -qq && apt install --yes --no-install-recommends wget \
        ca-certificates gnupg curl software-properties-common build-essential && \
    wget -q -O- https://eddelbuettel.github.io/r2u/assets/dirk_eddelbuettel_key.asc \
        | tee -a /etc/apt/trusted.gpg.d/cranapt_key.asc && \
    echo "deb [arch=amd64] https://r2u.stat.illinois.edu/ubuntu jammy main" \
        > /etc/apt/sources.list.d/cranapt.list && \
    apt update -qq && \
    echo "Package: *" > /etc/apt/preferences.d/99cranapt && \
    echo "Pin: release o=CRAN-Apt Project" >> /etc/apt/preferences.d/99cranapt && \
    echo "Pin: release l=CRAN-Apt Packages" >> /etc/apt/preferences.d/99cranapt && \
    echo "Pin-Priority: 700"  >> /etc/apt/preferences.d/99cranapt && \
    apt install --yes --no-install-recommends python3-dbus python3-gi python3-apt && \
    Rscript -e 'install.packages("bspm")' && \
    RHOME=$(R RHOME) && \
    echo "suppressMessages(bspm::enable())" >> ${RHOME}/etc/Rprofile.site && \
    echo "options(bspm.version.check=FALSE)" >> ${RHOME}/etc/Rprofile.site
# Copy the required R packages
COPY ./src/helper_scripts/install_r_pkgs.R .
# Install R packages - the tidyverse is already compiled
RUN Rscript --vanilla ./install_r_pkgs.R

# Install Rust -- this is the same as the Rust image for debian.
# We need this uglyness as we're root, so installing rust in the default way does not work so well.
ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH \
    RUST_VERSION=1.78.0

RUN set -eux; \
    dpkgArch="$(dpkg --print-architecture)"; \
    case "${dpkgArch##*-}" in \
        amd64) rustArch='x86_64-unknown-linux-gnu'; rustupSha256='a3d541a5484c8fa2f1c21478a6f6c505a778d473c21d60a18a4df5185d320ef8' ;; \
        armhf) rustArch='armv7-unknown-linux-gnueabihf'; rustupSha256='7cff34808434a28d5a697593cd7a46cefdf59c4670021debccd4c86afde0ff76' ;; \
        arm64) rustArch='aarch64-unknown-linux-gnu'; rustupSha256='76cd420cb8a82e540025c5f97bda3c65ceb0b0661d5843e6ef177479813b0367' ;; \
        i386) rustArch='i686-unknown-linux-gnu'; rustupSha256='cacdd10eb5ec58498cd95dbb7191fdab5fa4343e05daaf0fb7cdcae63be0a272' ;; \
        *) echo >&2 "unsupported architecture: ${dpkgArch}"; exit 1 ;; \
    esac; \
    url="https://static.rust-lang.org/rustup/archive/1.27.0/${rustArch}/rustup-init"; \
    wget "$url"; \
    echo "${rustupSha256} *rustup-init" | sha256sum -c -; \
    chmod +x rustup-init; \
    ./rustup-init -y --no-modify-path --profile minimal --default-toolchain $RUST_VERSION --default-host ${rustArch}; \
    rm rustup-init; \
    chmod -R a+w $RUSTUP_HOME $CARGO_HOME; \
    rustup --version; \
    cargo --version; \
    rustc --version;
# --- End Rust installation

## Install fonts
RUN sudo add-apt-repository --yes universe && \
    apt-get update && \
    apt-get install --yes fonts-firacode

COPY ./src/requirements.txt /src/
# Install python 3.13 and tree
RUN add-apt-repository --yes ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install --yes python3.13 jq

# Install python packages
RUN curl --proto '=https' --tlsv1.2 -sSf https://bootstrap.pypa.io/get-pip.py | python3.13 && python3.13 -m pip install -r ./src/requirements.txt

# Make sure python is findable as 'python'
RUN ln -s "$(which python3.13)" "/usr/bin/python"

# Install rust dependencies
RUN cargo install xsv
RUN cargo install --git https://github.com/MrHedmad/fast-cohen.git

# Install kerblam
RUN cargo install --locked --git https://github.com/MrHedmad/kerblam.git

# Copy the rest of the files
COPY . .
