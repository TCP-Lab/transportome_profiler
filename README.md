# Transportome Profiler analysis

![Kerblam!](https://img.shields.io/badge/Kerblam!-v1.0.0-rc.1-blue?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAC0AAAAtCAMAAAANxBKoAAABlVBMVEUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADW1tYNDHwcnNLKFBQgIB/ExMS1tbWMjIufDQ3S0tLOzs6srKyioqJRUVFSS0o0MjIBARqPj48MC3pqaWkIB2MtLS1ybm3U1NS6uroXirqpqamYmJiSkpIPZ4yHh4eFhIV8fHwLWnuBe3kMC3cLCnIHBlwGBlgFBU8EBEVPRkICAi4ADRa+EhIAAAwJCQmJiYnQ0NDKysoZkMK2trYWhLOjo6MTeKMTd6KgoKCbm5uKiIaAgIAPDHhubm4JT20KCW0KCWoIS2cHBUxBQUEEAz9IQT4DAz0DKTpFPTgCAjcCASoBASAXFxcgGRa5ERG1ERGzEBCpDw+hDg4fFA2WDAyLCgouAQFaWloFO1MBHStWBATnwMkoAAAAK3RSTlMA7zRmHcOuDQYK52IwJtWZiXJWQgXw39q2jYBgE/j2187JubKjoJNLSvmSt94WZwAAAvlJREFUSMeF1GdXGkEUgOGliIgIorFH0+u7JBIChEgJamyJvWt6783eS8rvzszAusACvp88x4d7hsvsaqdU57h8oQnobGmtb6xMzwbOkV9jJdvWBRwf7e9uLyzs7B3+o7487miC+AjcvZ3rkNZyttolbKxPv2fyPVrKYKcPhp7oIpPv0FkGN5N5rmd7afAFKH0MH99DihrTK2j3RTICF/Pt0trPUr9AxXyXpkJ3xu6o97tgQJDQm+Xlt6E8vs+FfNrg6kQ1pOuREVSPoydf9YjLpg14gMW1X0IInGZ+9PWr0Xl+R43pxzgM3NgCiekvqfE50hFdT7Ly8Jbo2R/xWYNTl8Ptwk6lgsHUD+Ji2NMlBFZ8ntzZRziXW5kLZsaDom/0yH/G+CSkapS3CvfFCWTxJZgMyqbYVLtLMmzoVywrHaPrrNJX4IHCDyCmF+nXhHXRkzhtCncY+PMig3pu0FfzJG900RBNarTTxrTCEwne69miGV5k8cPst3wOHSfrmJmcCH6Y42NEzzXIX8EFXmFE/q4ZXJrKW4VsY13uzqivF74OD39CbT/0HV/1yQW9Xn8e1O0w+WAG0VJS4P4Mzc7CK+2B7jt6XtFYMhl7Kv4YWMKnsJkXZiW3NgQXxTEKamM2fL8EjzwGv1srykZveBULj6bBZX2Bwbs03cXTQ3HAb9FOGNsS4wt5fw9zv0q9oZo54Gf4UQ95PLbJj/E1HFZ9DRgTuMecPgjfUqlF7Jo1B9wX+JFxmMh7mAoGv9B1pkg2tDoVl7i3G8mjH1mUN3PaspJaqM1NH/sJq2L6QJzEZ4FTCRosuKomdxjYSofDs8DcRPZh8hQd5IbE3qt1ih+MveuVeP2DxOMJAlphgSs1mt3GVWO6yMNGUDZDi1uzJLDNqxbZDLab3mqQB5mExtLYrtU45L10qlfMeSbVQ91eFlfRmnclZyR2VcB5y7pOYhouuSvg2rxHCZG/HHZnsVkVtg7NmkdirS6LzbztTq1EPo9dXRWxqtP7D+wL5neoEOq/AAAAAElFTkSuQmCC&link=https%3A%2F%2Fgithub.com%2FMrHedmad%2Fkerblam)

This repository contains the code for the analysis on the expression profile of
the transportome in Cancer based on the [MTP-DB](https://github.com/TCP-Lab/MTP-DB).

> [!NOTE]
> Read the preprint here: [Profiling the Expression of Transportome Genes in cancer: A systematic approach](https://doi.org/10.1101/2023.07.18.549498)
> 
> It's potentially out of date.

The project follows the [Kerblam! standard](https://github.com/MrHedmad/kerblam).

## Running the analysis
You can run the analysis pipelines with Kerblam! and `docker`:
```bash
# Clone the repo
git clone git@github.com:TCP-Lab/transportome_profiler.git
cd ./transportome_profiler

kerblam data fetch # Fetch the input data not present in the repository
kerblam run <pipeline>
```

Kerblam! will build docker containers and run the analysis locally.
To run without docker, read below.

### Pipelines
The project currently encompasses the following pipelines:
- `heatmaps`: Create large heatmaps from the expression matrices by using GSEA
  on computed gene rankings, testing all possible gene lists that can be made
  from the `MTP-DB`.
  - The `test` profile makes this pipeline much faster by running on smaller
    (i.e. sampled) input data (~75% reduction in sample number, only 5000 random genes).

### Running locally without docker
You need some requirements to be installed before you can run the analysis locally:
- `R` version `4.3.0+`.
  - Install R requirements with `./src/helper_scripts/install_r_pkgs.R`.
- `Python` version `3.11+`.
  - Install python requirements with `pip install -r ./src/requirements.txt`.
- The `jq` utility (that you can find [here](https://jqlang.github.io/jq/)).
- The `xsv` program, required by [`metasplit`](https://github.com/MrHedmad/metasplit)
  (`sudo pacman -Syu xsv` on Arch, not packaged by Debian, but [this guide might be useful](https://lindevs.com/install-xsv-on-ubuntu).
  If you have `cargo` installed, you can simply run `cargo install xsv`).
- Follow the extra installation guide for [`generanker`](https://github.com/TCP-Lab/gene_ranker) (namely installing [fast-cohen](https://github.com/MrHedmad/fast-cohen))
- The `xls2csv` utility (on arch `yay -Syu perl-xls2csv`)
- A series of R packages that can be installed with `Rscript ./src/helper_scripts/install_R_pkgs.R`
- Quite a bit of RAM (some steps require > 50 Gb of RAM) and time.
  Override `N_THREADS` (with `export N_THREADS=...`) to run with less threads.

If you have all the requirements, you can simply:
```bash
kerblam run <pipeline> --local
```

> [!IMPORTANT]
> The manuscript for this project is also available online in [this repository](https://github.com/TCP-Lab/transportome_profiler_paper).

