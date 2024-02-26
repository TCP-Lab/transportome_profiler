#!/usr/bin/env zsh

from pathlib import Path
import json
from bonsai.bonsai import Tree
from dataclasses import dataclass

def id_to_name(node_id: str, tree: Tree):
    node = tree.nodes[node_id]

    return node.name

@dataclass
class HashDiff:
    path_original: str
    path_hash: str
    content_original: str
    content_hash: str
    unique_original: str
    unique_hash: str

    id: str

def hash_node(node_id: str, tree: Tree) -> HashDiff:
    path = tree.get_path_of(node_id)
    path = [id_to_name(x, tree) for x in path]
    path = "/".join(path)

    contents = tree.nodes[node_id].data
    contents.sort()

    blob = f"{path}{'-'.join(contents)}"

    print(f"{blob}")

    return HashDiff(
        path_original=path,
        path_hash=hash(path),
        content_original='-'.join(contents),
        content_hash=hash(contents),
        unique_hash=hash(blob),
        unique_original=blob,
        id=node_id
    )

def check_diff(one: HashDiff, two: HashDiff):
    if one.unique_hash == two.unique_hash:
        return

    if one.content_hash != two.content_hash:
        print(f"Content {one.content_original} differes from {two.content_original}")

    if one.path_hash != two.path_hash:
        print(f"Path {one.path_original} differs from {two.path_hash}")

def main(args):

    with args.one.open("r") as one_conn, args.two.open("r") as two_conn:
        one = json.load(one_conn)
        one = Tree.from_node_json(one)
        two = json.load(two_conn)
        two = Tree.from_node_json(two)

    print(f"Loaded bonsai: one: {one} , two: {two}")

    print("Compute hashes")
    one_hashes = {}
    two_hashes = {}

    for node in one.nodes.values():
        t = hash_node(node.id, one)
        one_hashes[t.unique_hash] = t

    for node in two.nodes.values():
        two_hashes[hash_node(node.id, two)] = node.id

    print("Checking for differences...")

    all_hashes = []
    ids = set()
    for x in (list(one_hashes.values()) + list(two_hashes.values())):
        if x.unique_hash not in ids:
            all_hashes.append(x)
            ids.add(x.unique_hash)

    one_ids = [x.unique_hash for x in one_hashes.values()]
    two_ids = [x.unique_hash for x in two_hashes.values()]




if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("one", type=Path, help="First JSON to parse")
    parser.add_argument("two", type=Path, help="Second JSON to parse")

    args = parser.parse_args()

    main(args)
