#!/usr/bin/env python3
from __future__ import annotations

import json
import logging
import random
from collections import Counter
from copy import copy
from enum import Enum
from logging import StreamHandler
from pathlib import Path
from sqlite3 import Connection, connect

import numpy
import pandas as pd
from bonsai import Node, Tree, IdBuilder
from colorama import Back, Fore, Style
from tqdm import tqdm

random.seed(1)
numpy.random.seed(1)

builder = IdBuilder()

## >>>> Logging setup
class ColorFormatter(logging.Formatter):
    # Change this dictionary to suit your coloring needs!
    COLORS = {
        "WARNING": Fore.YELLOW,
        "ERROR": Fore.RED,
        "DEBUG": Style.BRIGHT + Fore.MAGENTA,
        "INFO": Fore.GREEN,
        "CRITICAL": Style.BRIGHT + Fore.RED,
    }

    def format(self, record):
        reset = Fore.RESET + Back.RESET + Style.NORMAL
        color = self.COLORS.get(record.levelname, "")
        if color:
            record.name = Style.BRIGHT + Fore.CYAN + record.name + reset
            if record.levelname != "INFO":
                record.msg = color + record.msg + reset
            record.levelname = color + record.levelname + reset
        return logging.Formatter.format(self, record)


log = logging.getLogger("Ariadne")
log.setLevel(logging.DEBUG)
log.propagate = False
# Keep this at DEBUG - set levels in handlers themselves

format = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
console_formatter = ColorFormatter(format)

stream_h = StreamHandler()
stream_h.setFormatter(console_formatter)
stream_h.setLevel(logging.INFO)

log.addHandler(stream_h)
## <<<< Logging setup


class PruneDirection(Enum):
    TOPDOWN = "topdown"
    BOTTOMUP = "bottomup"


def calc_similarity(node_a: set, node_b: set) -> float:
    # For compatibility
    node_a = set(node_a)
    node_b = set(node_b)

    # This is the Sorens-Dice coefficient. It was chosen because it is very easy to compute
    # and it's really similar to the Jaccard index anyway.
    # It empirically worked pretty well.
    return 2 * (len(node_a.intersection(node_b)) / (len(node_a) + len(node_b)))


def prune(tree: Tree, similarity: float, direction: PruneDirection) -> Tree:
    original_len = len(tree.nodes)
    log.info(f"Pruning {tree}.")

    reverse_sort = direction == PruneDirection.TOPDOWN

    def is_similar(node: Node, nodes: list[Node]) -> bool:
        for other in nodes:
            # This is reversed.
            # If the other node is the root, ignore it.
            if any([other.id == "0", node.id == "0"]):
                continue
            # If the two nodes are NOT similar, we can go on and
            # check the others.
            if calc_similarity(node.data, other.data) < similarity:
                continue
            # If they ARE similar, we can stop, and return TRUE
            return True
        return False

    cycle = 0
    pruned = True
    while pruned:
        pruned = False
        log.info(f"Prune cycle {cycle} -- {len(tree.nodes)} nodes in tree.")
        # Find all leaves
        leaves = tree.leaves()

        # Sort them
        leaves.sort(
            key=lambda x: tree.depth_of(x.id),
            reverse=reverse_sort,
        )

        # Prune
        for node in tqdm(leaves):
            other_nodes = list(copy(list(tree.nodes.values())))
            other_nodes.remove(node)
            if is_similar(node, other_nodes):
                log.debug(f"Pruned {node}")
                pruned = True
                tree.prune(node.id)

        cycle += 1

    len_diff = original_len - len(tree.nodes)
    log.info(
        f"Prune finished. Removed {len_diff} nodes ({round(len_diff / original_len * 100, 3)}% of total)"
    )

    return tree


def main(args: dict) -> None:
    log.info(f"Launching with args: {args}")

    log.info(f"Connecting to {args.database_path}...")
    connection: Connection = connect(args.database_path)

    # 1. Generate large tables
    log.info("Making large tables...")
    with args.basic_gene_lists.open("r") as stream:
        sets = json.load(stream)

    large_tables = make_large_tables(connection, sets)

    log.info(f"Made {len(large_tables)} large tables.")

    # 2. Generate lists from large tables
    log.info("Generating gene trees...")
    trees = {}
    for name, table in large_tables.items():
        log.info(f"Processing table {name}")
        tree = generate_gene_list_trees(
            table,
            name,
            min_pop_score=args.min_pop_score,
            min_set_size=args.min_set_size,
            min_recurse_set_size=args.min_recurse_set_size,
            recurse=not args.no_recurse,
        )
        trees[name] = tree

    # 3. Make the union of the genesets following the structure
    log.info("Pasting trees together...")
    large_tree = Tree(_id_fn=builder)
    for source, sink in tqdm(sets["structure"], desc="Merging"):
        if source == "root":
            large_tree.create_node(sink, None)
        else:
            large_tree.create_node(sink, large_tree.get_one_node_named(source).id)
        large_tree.paste(
            trees[sink],
            large_tree.get_one_node_named(sink).id,
            update_data=True,
        )

    if not args.no_prune:
        log.info("Pruning tree...")
        large_tree = prune(
            large_tree,
            similarity=args.prune_similarity,
            direction=PruneDirection(args.prune_direction),
        )

    large_tree.to_node_json(Path(args.out_json).open("w+"))
    large_tree.to_representation(Path(args.out_repr).open("w+"), force_uuid=True)

    log.info("Finished!")


def make_large_tables(conn: Connection, sets: dict) -> dict[pd.DataFrame]:
    """Generate large tables from a database and a list of genesets

    The database is seen as a series of tables that have to be row-wise joined
    to create larger tables to be queried for gene sets.

    The number and way to combine these tables is given by the `sets` param.

    Args:
        conn (Connection): A connection to the database
        sets (dict): A dictionary with two keys: "genesets" and "queries".
            The "genesets" key has to have a dictionary with as keys the names
            of the large tables, and as values list of table names from the
            database that have to be joined together to form a large table.
            The "queries" key has to have a dictionary with as keys table names
            of the database and as values SQL queries that return the table in
            question for joining.
            Tables in "genesets" that are not in "queries" will be retrieved
            with "SELECT * FROM {table_name};"

    Returns:
        dict[pd.DataFrame]: A dictionary with the same keys as the
            sets["genesets"] dictionary and as values pandas DataFrames with
            the resulting data.
    """
    assert type(sets.get("genesets", None)) is dict, "Genesets dictionary not found."
    queries = sets.get("queries", None)

    large_tables = {}
    for set_name, tables in sets["genesets"].items():
        loaded_tables = []
        for table_name in tables:
            query = (
                queries[table_name]
                if queries and table_name in queries
                else f"SELECT * FROM {table_name};"
            )
            log.debug(f"Loading table {table_name} with query {query}.")
            loaded_table = pd.read_sql(query, conn)
            loaded_table = loaded_table.reindex(sorted(loaded_table.columns), axis=1)
            loaded_tables.append(loaded_table)

        large_table = pd.concat(loaded_tables, ignore_index=True, sort=True)
        large_tables[set_name] = large_table.sort_values("ensg", axis=0)

    return large_tables


def generate_gene_list_trees(
    dataframe: pd.DataFrame,
    name: str,
    id_col: str = "ensg",
    min_pop_score: float = 0.5,
    min_set_size: int = 10,
    min_recurse_set_size: int = 40,
    recurse: bool = True,
) -> Tree:
    """Generate gene lists from a dataframe.

    Takes a dataframe with at least the ID_COL, then uses the other columns to
    generate possible gene lists with.

    Args:
        dataframe (pd.DataFrame): The dataframe to source
        id_col (str, optional): The column to use as IDs. Defaults to "ensg".
        min_pop_score (float, optional): Minimum percentage of non-NA values in
          a column to be considered for gene lists. Defaults to 0.5.
        min_set_size (int, optional): Minimum number of genes to produce a valid
          gene set. Defaults to 10.
        min_recurse_set_size (int, optional): If recurse is true, minimum parent
          gene set size to have before running recursion on it.
        recurse (bool, optional): Recurse of sub-dataframes? Defaults to TRUE
        name (str, optional): The name to give to the overall set of genesets.
          in other words, the name of the parent node for this geneset.

    Returns:
        Tree: A Tree structure of nodes, where each node contains the geneset
    """

    def generate_list(tree: Tree, father_node_id: str, frame: pd.DataFrame, layer: int):
        log.debug(f"Enumerating layer {layer}: {list(frame.columns)}")
        # This is the recursive wrapper

        valid_cols = []
        for current_col in sorted(list(frame.columns)):
            if sum(frame[current_col].isna()) / len(frame.index) > 1 - min_pop_score:
                log.debug(f"Layer {layer} -- col {current_col} ... SKIPPED (too empty)")
                continue
            valid_cols.append(current_col)

        iter_cols = valid_cols  # Remove ID col
        iter_cols.remove("ensg")
        for current_col in iter_cols:
            # Skip processing of id col
            if current_col == id_col:
                continue
            counts = Counter(frame[current_col].dropna())

            for value in sorted(set(counts.elements())):
                if counts[value] < min_set_size:
                    log.debug(
                        f"Layer {layer} -- col {current_col} -- value {value} ... SKIPPED (too small)"
                    )
                    continue

                putative_list = (
                    frame[id_col][frame[current_col] == value]
                    .drop_duplicates()
                    .to_list()
                )

                # Skip if the putative gene set is too small
                if len(putative_list) < min_set_size:
                    log.debug(
                        f"Layer {layer} -- col {current_col} -- value {value} ... SKIPPED (too small pure set)"
                    )
                    continue

                node_name = f"{current_col}::{value}"
                node_id = tree.create_node(
                    node_name, father_node_id, data=putative_list
                )

                if not recurse:
                    log.debug(
                        f"Layer {layer} -- col {current_col} -- value {value} ... ACCEPTED NR (id : {node_id})"
                    )
                    continue

                if len(putative_list) < min_recurse_set_size:
                    log.debug(
                        f"Layer {layer} -- col {current_col} -- value {value} ... ACCEPTED NR (id : {node_id})"
                    )
                    continue

                log.debug(
                    f"Layer {layer} -- col {current_col} -- value {value} ... ACCEPTED RC (id : {node_id})"
                )

                # Add back the ID col
                recurse_cols = list(frame.columns)
                recurse_cols.remove(current_col)
                new_data = frame.loc[frame[current_col] == value, recurse_cols]

                tree = generate_list(tree, node_id, new_data, layer=layer + 1)

        return tree

    all_colnames = list(dataframe.columns)
    # Remove the ID col
    all_colnames.remove(id_col)

    tree = Tree(_id_fn=builder)
    tree_root = tree.create_node(
        name, parent=None, data=dataframe[id_col].drop_duplicates().to_list()
    )

    subtree = generate_list(tree, tree_root, dataframe, 0)

    tree.paste(subtree, tree_root)

    log.debug(f"Generated {len(tree.all_nodes())} gene lists from table.")

    return tree


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "database_path", type=Path, help="Database used as annotation source"
    )
    parser.add_argument(
        "basic_gene_lists", type=Path, help="Gene lists to draw from the database"
    )
    parser.add_argument("out_json", type=Path, help="Output tree JSON representation")
    parser.add_argument("out_repr", type=Path, help="Output tree visual representation")

    parser.add_argument(
        "--min_pop_score",
        type=float,
        default=0.5,
        help="Minimum fraction of non-null values to consider cols",
    )
    parser.add_argument(
        "--min_set_size", type=int, default=10, help="Minimum generated set size"
    )
    parser.add_argument(
        "--min_recurse_set_size",
        type=int,
        default=40,
        help="Minimum size of set to recurse on",
    )
    parser.add_argument("--no_recurse", help="Suppress recursion", action="store_true")
    parser.add_argument(
        "--no_prune", help="Do not run pruning on the gene lists", action="store_true"
    )
    parser.add_argument(
        "--prune_similarity",
        type=float,
        help="Node similarity threshold for pruning",
        default=0.45,
    )
    parser.add_argument(
        "--prune_direction",
        choices=["topdown", "bottomup"],
        help="Direction to prune nodes in",
        default="bottomup",
    )
    parser.add_argument("--verbose", help="Increase verbosity", action="store_true")
    parser.add_argument(
        "--json",
        help="Dump a JSON of the tree instead of nested dirs",
        action="store_true",
    )

    args = parser.parse_args()

    if args.verbose:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    main(args)
