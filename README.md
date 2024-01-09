# Transportome Profiler analysis

![Kerblam!](https://img.shields.io/badge/Kerblam!-v0.2.3-blue?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAC0AAAAtCAMAAAANxBKoAAABlVBMVEUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADW1tYNDHwcnNLKFBQgIB/ExMS1tbWMjIufDQ3S0tLOzs6srKyioqJRUVFSS0o0MjIBARqPj48MC3pqaWkIB2MtLS1ybm3U1NS6uroXirqpqamYmJiSkpIPZ4yHh4eFhIV8fHwLWnuBe3kMC3cLCnIHBlwGBlgFBU8EBEVPRkICAi4ADRa+EhIAAAwJCQmJiYnQ0NDKysoZkMK2trYWhLOjo6MTeKMTd6KgoKCbm5uKiIaAgIAPDHhubm4JT20KCW0KCWoIS2cHBUxBQUEEAz9IQT4DAz0DKTpFPTgCAjcCASoBASAXFxcgGRa5ERG1ERGzEBCpDw+hDg4fFA2WDAyLCgouAQFaWloFO1MBHStWBATnwMkoAAAAK3RSTlMA7zRmHcOuDQYK52IwJtWZiXJWQgXw39q2jYBgE/j2187JubKjoJNLSvmSt94WZwAAAvlJREFUSMeF1GdXGkEUgOGliIgIorFH0+u7JBIChEgJamyJvWt6783eS8rvzszAusACvp88x4d7hsvsaqdU57h8oQnobGmtb6xMzwbOkV9jJdvWBRwf7e9uLyzs7B3+o7487miC+AjcvZ3rkNZyttolbKxPv2fyPVrKYKcPhp7oIpPv0FkGN5N5rmd7afAFKH0MH99DihrTK2j3RTICF/Pt0trPUr9AxXyXpkJ3xu6o97tgQJDQm+Xlt6E8vs+FfNrg6kQ1pOuREVSPoydf9YjLpg14gMW1X0IInGZ+9PWr0Xl+R43pxzgM3NgCiekvqfE50hFdT7Ly8Jbo2R/xWYNTl8Ptwk6lgsHUD+Ji2NMlBFZ8ntzZRziXW5kLZsaDom/0yH/G+CSkapS3CvfFCWTxJZgMyqbYVLtLMmzoVywrHaPrrNJX4IHCDyCmF+nXhHXRkzhtCncY+PMig3pu0FfzJG900RBNarTTxrTCEwne69miGV5k8cPst3wOHSfrmJmcCH6Y42NEzzXIX8EFXmFE/q4ZXJrKW4VsY13uzqivF74OD39CbT/0HV/1yQW9Xn8e1O0w+WAG0VJS4P4Mzc7CK+2B7jt6XtFYMhl7Kv4YWMKnsJkXZiW3NgQXxTEKamM2fL8EjzwGv1srykZveBULj6bBZX2Bwbs03cXTQ3HAb9FOGNsS4wt5fw9zv0q9oZo54Gf4UQ95PLbJj/E1HFZ9DRgTuMecPgjfUqlF7Jo1B9wX+JFxmMh7mAoGv9B1pkg2tDoVl7i3G8mjH1mUN3PaspJaqM1NH/sJq2L6QJzEZ4FTCRosuKomdxjYSofDs8DcRPZh8hQd5IbE3qt1ih+MveuVeP2DxOMJAlphgSs1mt3GVWO6yMNGUDZDi1uzJLDNqxbZDLab3mqQB5mExtLYrtU45L10qlfMeSbVQ91eFlfRmnclZyR2VcB5y7pOYhouuSvg2rxHCZG/HHZnsVkVtg7NmkdirS6LzbztTq1EPo9dXRWxqtP7D+wL5neoEOq/AAAAAElFTkSuQmCC&link=https%3A%2F%2Fgithub.com%2FMrHedmad%2Fkerblam)

This repository contains the code for the analysis on the expression profile of
the transportome in Cancer based on the [MTP-DB](https://github.com/TCP-Lab/MTP-DB).

> Read the preprint here: [Profiling the Expression of Transportome Genes in cancer: A systematic approach](https://doi.org/10.1101/2023.07.18.549498)

This is a two-step process. The database is queried for information by the
script in `src/geneset_maker`. The algorithm generates gene sets ready for use
by GSEA. We then generate a series of DEG tables with `src/run_dea` based on the
comparison of gene expression from TCGA (cancer) and GTEx (healthy) tissues.
Finally, GSEA is called by `src/gsea_runner` in a pre-ranked manner on all the
DEG tables with all of the genesets, making enrichment tables.
The enrichment tables are then processed by a script in `src/gsea_runner` to
make enrichment plots.

## Running the analysis
You can run the analysis in two ways: locally or in a Docker container.
In both cases, you must first clone and link the repository locally:
```bash
# Clone the repo
git clone git@github.com:TCP-Lab/transportome_profiler.git
cd ./transportome_profiler
# Make housekeeping directories. See the following note.
./link A_FOLDER_THAT_WILL_CONTAIN_THE_DATA
```
The variable `A_FOLDER_THAT_WILL_CONTAIN_THE_DATA` should be the directory where
you want the (rather bulky) data to live.
The directory can be anywhere you want, but
**DO NOT CHOOSE `./transportome_profiler/data`** as a data directory!
It is the place where your actual data directory will be linked to by `link`.

### Running in Docker
You will need to have `docker` installed. Once you do, you can either pull the
remote container or build a new one locally.

The script `run_docker` in the root of the repository helps you do just that:
- There is not (**yet!**) a remote docker container for the analysis. When there will be, run `run_docker --pull` to run the analysis in that.
- To build and run a docker container locally, run `run_docker`. The docker container will be built and executed immediately after.

Currently, docker will save everything as root. Run `chown ./data` after it is finished to re-own the files that it generated. Sorry for the inconvenience!

### Running locally
You need some requirements to be installed before you can run the analysis locally:
- `R` version `4.3.0`.
- `Python` version `3.11`.
- The `tree` utility (`sudo apt install tree` on Debian-like or `sudo pacman -Syu tree` on Arch).
- The `xsv` program, required by [`metasplit`](https://github.com/MrHedmad/metasplit) (`sudo pacman -Syu xsv` on Arch, not packaged by Debian, but [this guide might be useful](https://lindevs.com/install-xsv-on-ubuntu). If you have `cargo` installed, you can simply run `cargo install xsv`).
- A series of R packages that can be installed with `Rscript ./src/helper_scripts/install_R_pkgs.R`
- Quite a bit of RAM (the full analysis takes > 50 Gb of RAM) and time. If you want to run with less ram (but slower), override the `make` variable `split_threads` with `make split_threads=1`. By default this is 3. If you have a ton of memory and want to run more threads in parallel increase this value. 

`make` handles making and using a Python `venv` when appropriate.

If you have all the requirements, you can simply:
```bash
# Run the analysis
make all
# or just `make`, since `all` is the default.
```

#### A note for developers
This is the same workflow you use to start working on the project on a new PC. You can then simply re-run `make all` as you work.

### Other make targets
There are several make targets if you do not want to rerun all the steps of the analysis (like `all` does):
- `make all`: Runs the whole analysis - from the retrieval of the data to the generation of the final output.
- `make clean`: Cleanup every data file and the output.
- `make restart`: Cleanup input and output data and the virtual enviroment. To restart the analysis from scratch.
- `make scrub`: Cleanup just like `restart`, but also remove the soft links made by `./link` and the corresponding data folders. This should leave your disks squeaky clean, just like before you started the analysis.
- `make paper`: Runs the analysis and rebuilds the paper.
- `make thin_paper`: Rebuild the paper, without running the analysis. Figures that are not found are replaced by placeholders. Useful when you just want to write the paper without having to run the whole shebang.

# The Manuscript: Transportome Profiler
To compile the manuscript, you will need a local installation of `texlive` and `biber`.
Biber should be packaged in most distros. You can google how to install it.
`texlive` is a bit harder.
For Arch users, follow the [archwiki article](https://wiki.archlinux.org/title/TeX_Live).
In general, the recommended way to install `texlive` is to follow the [official TeXlive guide](https://tug.org/texlive/quickinstall.html).

The instruction above are copied here, but you should check if they changed:
```bash
cd /tmp/
wget https://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
zcat < install-tl-unx.tar.gz | tar xf -
cd install-tl-*
# Install only the medium image
perl ./install-tl --no-interaction --scheme=medium
```

You will probably need to install extra packages to compile successfully.
The command to do this for `tlmgr` is in `.../transportome_profiler/paper/install_pkgs.sh`:
```bash
sudo ./install_pkgs.sh
```
I am a weird man online. Suffice to say you should not `sudo` anything without inspecting it manually beforehand.

Once everything is installed, you should simply run `make paper`. `make` will run the analysis (to generate the figures for the manuscript), and then generate the paper itself. Keep reading if you want to *just* compile the manuscript, and not run the analysis proper.

## Just edit the manuscript
If you just want to edit the manuscript, but do not care about restarting the analysis, you can simply run `make thin_paper`. The paper will be built (with no figures) and saved to `/paper/src/main.pdf`. 

If you have access to pre-built images (from e.g. a previous run) you can just drop them in `/paper/src/resources/figures/generated` and even `make thin_paper` will pick them up and isert them in the output.

