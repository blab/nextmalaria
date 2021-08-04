"""
Given auspice-ready JSONs (i.e. from `augur export`) produce an
auspice-compatable JSON which hides the root node). Modified from
https://github.com/seattleflu/augur-build/blob/master/scripts/annotate_hidden_nodes.py
written by James Hadfield, July 2019
"""
import argparse
import json

## Magics / hardcoded parameters
STUB_DIV_LENGTH = 0.00001
#STUB_TIME_LENGTH = 0.05 # ~18 days

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract sample sequences by name",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, metavar="JSON", help="Unified JSON from augur export v2")
    parser.add_argument("--output", required=True, metavar="JSON", help="Unified JSON for auspice")
    return parser.parse_args()

def make_stub(node):
    stub = {
        "name": "stub",
        "node_attrs": {
            #"num_date": {
            #    "value": node["node_attrs"]["num_date"]["value"] - STUB_TIME_LENGTH
            #},
            "div": node["node_attrs"]["div"] - STUB_DIV_LENGTH,
            "hidden": "always"

        },
        "children": [node]
    }
    return stub

def hide_root(tree):
    tree["node_attrs"]["hidden"] = "always"
    root = tree["children"][0]
    root["node_attrs"]["hidden"] = "always"
    restOfTree = tree["children"][1]
    new_tree = make_stub(restOfTree)
    return new_tree

def flatten_tree(tree):
    flat = []
    stack = [tree]
    while len(stack) > 0:
        node = stack.pop()
        flat.append(node)
        if "children" in node:
            for child in node["children"]:
                stack.append(child)
    return flat


def shift_hidden_div(flat_tree):
    """
    Shift the num_date and div attributes of the hidden nodes so that the
	non-hidden nodes take up the entire screen in auspice. (Auspice calculates
	the view based on all the nodes, whether they are hidden or not, so a deep
	MRCA -- even if hidden -- doesn't look great.) This uses the trick that both
	flat_tree & tree contain references to the same object, so modifications to
	one can affect the other. Note: plenty of optimisation potential here
    """
    # find the date of the earliest node which has a cluster and modify all
    # nodes earlier than that to have that date
    #min_date = min([n["node_attrs"]["num_date"]["value"] for n in flat_tree if get_cluster(n)])
    #for n in flat_tree:
    #    if n["node_attrs"]["num_date"]["value"] < min_date:
            # setting the `num_date` here will also overwrite any confidences if they are set
    #        n["node_attrs"]["num_date"] = {"value": min_date}
    div = flat_tree[0]["node_attrs"]["div"] + STUB_DIV_LENGTH
    # modify clusters to each have divergence starting from 0
    for n in flat_tree:
        if "hidden" in n["node_attrs"]:
            n["node_attrs"]["div"] = -STUB_DIV_LENGTH
        else:
            n["node_attrs"]["div"] -= div


if __name__ == '__main__':
    args = parse_args()

    with open(args.input, "rU") as fh:
        unified = json.load(fh)

    tree = unified["tree"]
    new_tree = hide_root(tree)
    flat_tree = flatten_tree(new_tree)
    shift_hidden_div(flat_tree)

    unified["tree"] = new_tree
    with open(args.output, "w") as fh:
        json.dump(unified, fh, indent=2, sort_keys=True)
