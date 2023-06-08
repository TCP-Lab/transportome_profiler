# Transportome Profiler

This repository contains the code for the analysis on the expression profile of the transportome in Cancer based on the [MTP-DB](https://github.com/CMA-Lab/MTP-DB).

This is a two-step process. The database is queried for information by the script in `src/geneset_maker`. The algorithm generates gene sets ready for use by GSEA. We then generate a series of DEG tables with `src/run_dea` based on the comparison of gene expression from TCGA (cancer) and GTEx (healthy) tissues. Finally, GSEA is called by `src/gsea_runner` in a pre-ranked manner on all the DEG tables with all of the genesets, making enrichment tables.
The enrichment tables are then processed by a script in `src/gsea_runner` to make enrichment plots.

## Running the analysis
A docker container is coming soon(tm). For now, follow these steps:

You need some requirements to be installed before you can run the analysis:
- `R`
- `Python`
- The `tree` utility (`sudo apt install tree` on Debian-like or `sudo pacman -Syu tree` on Arch).
- The `xsv` program, required by [`metasplit`](https://github.com/MrHedmad/metasplit) (`sudo pacman -Syu xsv` on Arch, not packaged by Debian, but [this guide might be useful](https://lindevs.com/install-xsv-on-ubuntu). If you have `cargo` installed, you can simply run `cargo install xsv`).
- A series of R packages that can be installed with `Rscript ./src/helper_scripts/install_R_pkgs.R`
- Quite a bit of RAM (the full analysis takes > 50 Gb of RAM) and time. If you want to run with less ram (but slower), override the `make` variable `split_threads` with `make split_threads=1`. By default this is 3. If you have a ton of memory and want to run more threads in parallel increase this value. 

`make` handles making and using a Python `venv` when appropriate.

If you have all the requirements, you can:
```bash
# Clone the repo
git clone git@github.com:CMA-Lab/transportome_profiler.git
cd ./transportome_profiler
# Make housekeeping directories. See the following note.
./link A_FOLDER_THAT_WILL_CONTAIN_THE_DATA
# Run the analysis
make all
# or just `make`, since `all` is the default.
```

The variable `A_FOLDER_THAT_WILL_CONTAIN_THE_DATA` should be the directory where you want the (rather bulky) data to live. The directory can be anywhere you want, but care must be taken that **YOU DO NOT CHOOSE `./transportome_profiler/data`** as a data directory! It is the place where your actual data directory will be linked to. Beware of paradoxes!

### A note for developers
This is the same workflow you use to start working on the project on a new PC. You can then simply re-run `make all` as you work.

## Other make targets
There are several make targets if you do not want to rerun all the steps of the analysis (like `all` does):
- `make all`: Runs the whole analysis - from the retrieval of the data to the generation of the final output.
- `make retrieve`: Just download the remote data files.
- `make analyze`: Retrieve + analyze the data, but do not make the latex paper.
- `make clean`: Cleanup every data file and the output.
- `make scrub`: Cleanup just like `clean`, but also remove the soft links made by `./link` and the corresponding data folders. This should leave your disks squeaky clean, just like before you started the analysis.
