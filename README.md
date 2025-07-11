# Transportome Profiler analysis

![Kerblam!](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2FTCP-Lab%2Ftransportome_profiler%2Frefs%2Fheads%2Fmain%2Fkerblam.toml&query=%24.meta.version&prefix=v.&logo=Rocket&logoColor=red&label=Kerblam!)

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

