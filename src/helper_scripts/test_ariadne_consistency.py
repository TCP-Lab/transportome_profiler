import functools as ft
import json
import logging
import subprocess
import tempfile
import uuid

import bonsai
from tqdm import tqdm

log = logging.getLogger(__name__)


class IncongruencyError(Exception):
    def __init__(self, location: str, inconsistency: str):
        self.location = location
        self.inconsistency = inconsistency

    def __str__(self):
        return f"Inconsistency error at {self.location}: {self.inconsistency}"


def launch_ariadne(exec, args, target_file):
    log.debug(f"Running Ariadne to create {target_file}")
    all_args = [exec]
    all_args.extend(args)
    all_args.extend([target_file, "/dev/null"])
    ret = subprocess.run(all_args, capture_output=True)
    if ret.returncode != 0:
        raise RuntimeError(
            f"Ariadne failed to run:\n{ret.stderr.decode("UTF-8")}",
        )

    with open(target_file, "r") as f:
        return json.load(f)


def which(iterable):
    enum = [i for i, x in enumerate(iterable) if x]
    if len(enum) == 1:
        return enum[0]
    return enum


def find(items, fn):
    """Find one and only one item in a list

    The fn is a function that takes each item and returns True or False.
    """
    checks = [fn(x) for x in items]
    if sum(checks) == 1:
        return items[which(checks)]
    elif sum(checks) > 1:
        raise RuntimeError("Found more than one matching item")
    else:
        raise RuntimeError("Item not found")


def get_children_named(tree: bonsai.Tree, parent_id, name):
    childrens = tree.get_direct_children(parent_id)
    return find(childrens, lambda x: x.name == name)


def check_two_way_congruency(list_1, list_2):
    # Could probably forgo one check and check for len but whatever
    if not all([x in list_2 for x in list_1]):
        raise RuntimeError(
            f"Did not find some items in list 1 inside list 2. Length of one: {len(list_1)} vs two {len(list_2)}"
        )
    if not all([x in list_1 for x in list_2]):
        raise RuntimeError(
            f"Did not find some items in list 2 inside list 1. Length of one: {len(list_1)} vs two {len(list_2)}"
        )
    return True


def get_location_of(tree: bonsai.Tree, node_id):
    path = tree.get_path_of(node_id)
    names = [tree.nodes[x].name for x in path]
    spath = ""
    for name in names:
        spath += name + "///"

    return spath.strip("/")


def test_congruency(one: bonsai.Tree, two: bonsai.Tree):
    def compare_childrens(parent_one_id, parent_two_id):
        first_children = one.get_direct_children(parent_one_id)
        second_children = two.get_direct_children(parent_two_id)
        first_children_names = [one.nodes[x.id].name for x in first_children]
        second_children_names = [two.nodes[x.id].name for x in second_children]

        log.debug(
            f"Comparing {one.nodes[parent_one_id].name} with {two.nodes[parent_two_id].name}"
        )

        if not first_children and not second_children:
            # There are no nodes to compare - the two parents are leaves.
            # This is the safe exit condition of the recursion.
            log.debug("Nodes are leaves. Returning.")
            return
        elif not first_children:
            # There are second children but not first ones. Panic!
            raise IncongruencyError(
                get_location_of(one, parent_one_id),
                f"This node in the first tree should have children: {second_children_names}",
            )
        elif not second_children:
            # Same, but reversed
            raise IncongruencyError(
                get_location_of(two, parent_two_id),
                f"This node in the second tree should have children: {first_children_names}",
            )

        if not len(first_children) == len(second_children):
            raise IncongruencyError(
                f"{one.nodes[parent_one_id].name} and {two.nodes[parent_two_id].name}",
                f"Different number of children, {len(first_children)} vs {len(second_children)}: {first_children_names} vs {second_children_names}",
            )

        # It does not matter which of the two lists we check - they have the same length.
        for child in first_children:
            try:
                sibling = find(second_children, lambda x: x.name == child.name)
            except RuntimeError as e:
                raise IncongruencyError(
                    f"{get_location_of(one, child.id)}",
                    f"could not locate sibling node: {str(e)}",
                )

            try:
                check_two_way_congruency(child.data, sibling.data)
            except RuntimeError as e:
                raise IncongruencyError(f"{get_location_of(one, child.id)}", str(e))

            # We are done comparing these children: they are identical.
            # now, check in turn their children.
            log.debug("Moving down...")
            compare_childrens(child.id, sibling.id)

    # Get the roots of the trees and start checking
    compare_childrens(one.root.id, two.root.id)


def main(args, extra):
    ariadne = ft.partial(launch_ariadne, exec=args.exec, args=extra)
    with tempfile.TemporaryDirectory() as scratch:
        # We add the "out_json" and "out_repr" files
        # the latter is just redirected to be deleted
        template = f"{scratch}/{uuid.uuid4()}"
        template_tree = bonsai.Tree.from_node_json(ariadne(target_file=template))

        for _ in tqdm(range(args.trials)):
            response_path = f"{scratch}/{uuid.uuid4()}"
            response_tree = bonsai.Tree.from_node_json(
                ariadne(target_file=response_path)
            )

            test_congruency(template_tree, response_tree)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("exec", help="Path to the executable file")
    parser.add_argument("--trials", help="Number of trials to check", type=int)
    parser.add_argument("--verbose", help="Increase verbosity", action="store_true")

    args, extra = parser.parse_known_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    main(args, extra)
