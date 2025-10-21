#!/usr/bin/env python

import os
import re
import tarfile
import tempfile
from pathlib import Path
from shutil import copy


def main(path: Path, temp_dir: Path):
    """Point this to the tranportome_profiler_results.tar file to get images

    It creates the 'images' folder with the correct file names to be used
    in the tranportome profiler paper.
    """

    # Extract the tar archive to the temp dir
    with tarfile.open(path) as tar:
        tar.extractall(temp_dir, filter="data")

    # Extract all inner tar archives
    result_regex = re.compile(r"packaged_output_(.*?).tar.gz")
    results = {}
    for file in temp_dir.iterdir():
        match = result_regex.match(file.name)
        if not match:
            continue
        result_name = match.group(1)
        print(f"Looking into {file} ({result_name})")
        results[result_name] = file

    for name, tar_path in results.items():
        extract_images(output_dir=Path(path).parent, target=tar_path, slug=name)

    ## We're done - extract the extra plots
    extra_plots = Path(temp_dir) / "extra_plots.tar.gz"
    tar = tarfile.open(extra_plots, mode="r:gz")
    with tempfile.TemporaryDirectory() as tmp:
        tar.extractall(tmp, filter="data")
        for target, dist in EXTRA_PLOT_TARGETS.items():
            full_dist = Path(path).parent / dist
            os.makedirs(full_dist.parent, exist_ok=True)
            copy(Path(tmp) / target, full_dist)


TARGETS = {
    "data/out/figures/deregulation_heatmap.png": "images/heatmaps/tcga_gtex/deregulation_heatmap_{slug}.png",
    "data/out/figures/geo_deregulation_heatmap.png": "images/heatmaps/geo/geo_deregulation_heatmap_{slug}.png",
    "data/out/figures/geo_upset.png": "images/upset/upset_geo_{slug}.png",
    "data/out/figures/top_disregulation_thr_1.5_set_channels.png": "images/top_disregulation/tcga_gtex/td_1.5_channels_{slug}.png",
    "data/out/figures/top_disregulation_thr_1.5_set_transporters.png": "images/top_disregulation/tcga_gtex/td_1.5_transporters_{slug}.png",
    "data/out/figures/top_disregulation_thr_1.5_set_whole_transportome.png": "images/top_disregulation/tcga_gtex/td_1.5_whole_transportome_{slug}.png",
    "data/out/figures/top_disregulation_thr_1_set_channels.png": "images/top_disregulation/tcga_gtex/td_1_channels_{slug}.png",
    "data/out/figures/top_disregulation_thr_1_set_transporters.png": "images/top_disregulation/tcga_gtex/td_1_transporters_{slug}.png",
    "data/out/figures/top_disregulation_thr_1_set_whole_transportome.png": "images/top_disregulation/tcga_gtex/td_1_whole_transportome_{slug}.png",
    "data/out/figures/top_disregulation_thr_2_set_channels.png": "images/top_disregulation/tcga_gtex/td_2_channels_{slug}.png",
    "data/out/figures/top_disregulation_thr_2_set_transporters.png": "images/top_disregulation/tcga_gtex/td_2_transporters_{slug}.png",
    "data/out/figures/top_disregulation_thr_2_set_whole_transportome.png": "images/top_disregulation/tcga_gtex/td_2_whole_transportome_{slug}.png",
    "data/out/figures/upset.png": "images/upset/upset_{slug}.png",
    "data/out/figures/correlation.png": "images/corr/correlation_{slug}.png",
}

EXTRA_PLOT_TARGETS = {
    "data/out/figures/channels_expression_upset.png": "images/extra_plots/channels_expression_upset.png",
    "data/out/figures/expression_means.png": "images/extra_plots/expression_means.png",
    "data/out/figures/expression_means_GTEX_only.png": "images/extra_plots/expression_means_GTEX_only.png",
    "data/out/figures/expression_means_TCGA_only.png": "images/extra_plots/expression_means_TCGA_only.png",
    "data/out/figures/transporters_expression_upset.png": "images/extra_plots/transporters_expression_upset.png",
    "data/out/figures/whole_transportome_expression_upset.png": "images/extra_plots/whole_transportome_expression_upset.png",
    "data/out/figures/plot_commonality_set_channels.png": "images/extra_plots/plot_commonality_set_channels.png",
    "data/out/figures/plot_commonality_set_transporters.png": "images/extra_plots/plot_commonality_set_transporters.png",
    "data/out/figures/plot_commonality_set_whole_transportome.png": "images/extra_plots/plot_commonality_set_whole_transportome.png",
}


def extract_images(output_dir: Path, target: Path, slug: str):
    with tempfile.TemporaryDirectory() as tmp:
        tar = tarfile.open(target, mode="r:gz")
        tar.extractall(tmp, filter="data")
        for source, dest in TARGETS.items():
            full_target = output_dir / dest.format(slug=slug)
            os.makedirs(full_target.parent, exist_ok=True)
            copy(Path(tmp) / source, full_target)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Path to the output tar archive")

    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as dir:
        main(args.path, Path(dir))
