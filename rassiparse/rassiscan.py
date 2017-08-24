#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys

from natsort import natsorted
import yaml

from rassiparse import parse_rassi, handle_rassi, set_mo_images

def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("rassi",
        help="Scan logfile containing the rassi outputs.")
    parser.add_argument("mofns",
        help=".yaml file created by mopicgen.py --dumpfns.")

    return parser.parse_args(args)


def handle_scan_step(spinfree_states, mo_fns):
    as_dicts = dict()
    # Iterave over all states for a given step
    for i, sfs in enumerate(spinfree_states):
        sfs.set_mo_images(mo_fns)
        as_dicts[i] = sfs.as_dict()

    return as_dicts


def run():
    args = parse_args(sys.argv[1:])

    # Split into single &rassi chunks
    rassi_chunk_regex = "Start Module: rassi(.+?)Stop Module:  rassi"
    with open(args.rassi) as handle:
        text = handle.read()
    rassi_chunks = re.findall(rassi_chunk_regex, text, re.DOTALL)
    sfs_lists = list()
    for rc in rassi_chunks:
        sf_states, trans_dict = parse_rassi(rc)
        handle_rassi(sf_states, trans_dict)
        sfs_lists.append(sf_states)

    # Load mo names
    with open(args.mofns) as handle:
        mo_fns_dict = yaml.load(handle.read())
    mo_fn_lists = [natsorted(mo_fns_dict[i])
                    for i in range(1, len(rassi_chunks)+1)
    ]
    assert(len(sfs_lists) == len(mo_fn_lists))

    # Iterate over all steps in the scan
    faw_dict = dict()
    for i, (spinfree_states, mo_fns) in enumerate(zip(sfs_lists, mo_fn_lists)):
        faw_dict[i] = handle_scan_step(spinfree_states, mo_fns)

    out_fn = "rassiscan.yaml"
    with open(out_fn, "w") as handle:
        handle.write(yaml.dump(faw_dict))


if __name__ == "__main__":
    run()
