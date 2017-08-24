#!/usr/bin/env python3

import argparse
import logging
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HARTREE2EV = 27.211386

def parse_roots_ens(text, regex):
    roots, ens = zip(
        *[(int(root), float(energy))
          for root, energy in re.findall(regex, text)
         ]
    )
    max_root = max(roots)
    ens = np.array(ens).reshape((-1, max_root))
    ens -= ens.min()
    ens *= HARTREE2EV

    return ens, max_root

def parse_reference_weights(text):
    # Reference weight:           0.55837
    regex = "Reference weight:\s*([\d\.]+)"
    ref_weights = np.array(
                    re.findall(regex, text),
                    dtype=np.float
    )
    return ref_weights


def get_rasscf_energies(text):
    # RASSCF root number  4 Total energy:  -1232.07159976
    regex = "RASSCF root number\s*(\d+).+?([\-\d\.]+)"
    return parse_roots_ens(text, regex)


def get_caspt2_energies(text):
    #::    CASPT2 Root  1     Total energy:  -1236.87839955
    regex = "::\s*CASPT2 Root\s*(\d+)\s*Total energy:\s*([\d\-\.]+)"
    ens, roots = parse_roots_ens(text, regex)
    ref_weights = parse_reference_weights(text)
    return ens, ref_weights


def get_mscaspt2_energies(text):
    #::    MS-CASPT2 Root  1     Total energy:  -1236.87839291
    regex = "::\s*MS-CASPT2 Root\s*(\d+)\s*Total energy:\s*([\d\-\.]+)"
    ens, roots = parse_roots_ens(text, regex)
    ref_weights = parse_reference_weights(text)
    return ens, ref_weights


def get_natural_bond_orders(text, roots, bond_pairs):
    regex = "Atom A       Atom B       Bond Order\n(.+?)\-\-"
    nbo_analysis = re.findall(regex, text, re.DOTALL)
    # Delete all "|" characters
    nbo_analysis = [nboa.replace("|", "").split() for nboa in nbo_analysis]
    all_bond_inds = list()
    all_bos = list()
    for nboa in nbo_analysis:
        assert(len(nboa) % 3 == 0)
        from_atoms, to_atoms, bos = zip(
            *[nboa[i:i+3] for i in range(0, len(nboa), 3)]
        )
        splt = lambda lst: [re.match("([A-Z]+)(\d+)", i).groups()
                            for i in lst]
        to_ints = lambda lst: [int(i) for i in lst]
        from_atoms, from_inds = zip(*splt(from_atoms))
        to_atoms, to_inds = zip(*splt(to_atoms))
        from_inds = to_ints(from_inds)
        to_inds = to_ints(to_inds)
        bos = [float(bo) for bo in bos]
        bond_inds = tuple(zip(from_inds, to_inds))
        all_bond_inds.append(bond_inds)
        all_bos.append(bos)

    all_matched_bos = list()
    for bond_inds, bos in zip(all_bond_inds, all_bos):
        matched_bos = [None for bp in bond_pairs]
        inds = [(i, bond_inds.index(indss))
                for i, indss in enumerate(bond_pairs)
                if indss in bond_inds
        ]
        for i, index in inds:
            matched_bos[i] = bos[index]
        all_matched_bos.append(matched_bos)

    # Split into roots
    per_root = list()
    for i in range(roots):
        per_root.append(
            [all_matched_bos[j] for j in range(i, len(all_matched_bos), 4)]
        )
    columns = [f"{bp[0]}-{bp[1]}" for bp in bond_pairs]
    per_root_dfs = [pd.DataFrame(pr, columns=columns) for pr in per_root]

    return per_root_dfs


def check_cmocorr(text):
    logging.warning("--cmocorr tested only with RAS2!")
    regex = "Analyzing Orbital spaces\.(.+?)Analyzing Orbitals\."
    orbital_spaces_blocks = re.findall(regex, text, re.DOTALL)
    all_ras2_overlaps = list()
    for i, osb in enumerate(orbital_spaces_blocks):
        #Ras 2           1   6   6   6.000
        #                2   4   4   4.000
        #                   10  10  10.000
        regex = "Ras 2\s*(.+?)\n\n"
        #regex = "Ras 2(.+?)Secondary"
        ras2_overlaps = re.search(regex, osb, re.DOTALL).groups()[0]
        ras2_overlaps = ras2_overlaps.split("\n")
        # Dropt the last line, that is a sum of the previous lines
        ras2_overlaps = [_.strip().split() for _ in ras2_overlaps][:-1]
        all_ras2_overlaps.extend(ras2_overlaps)

    columns = ("sym", "ref", "chk", "Overlap")
    df = pd.DataFrame(all_ras2_overlaps, columns=columns)
    print("cmocorr - ras2 overlaps")
    print(df)


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("fn")
    parser.add_argument("--nbo", type=int, nargs="+")
    parser.add_argument("--cmocorr", action="store_true")
    parser.add_argument("--caspt2", action="store_true")
    parser.add_argument("--mscaspt2", action="store_true")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    with open(args.fn) as handle:
        text = handle.read()

    rasscf_ens, roots = get_rasscf_energies(text)
    ens_df = pd.DataFrame(rasscf_ens)
    ens_df.to_csv("rasscf_ens.csv", index=False)

    subplots = 1
    dfs_to_plot = [ens_df, ]
    titles = ["SA-RASSCF", ]
    ref_weights = None
    if args.caspt2:
        caspt2_ens, ref_weights = get_caspt2_energies(text)
        caspt2_df = pd.DataFrame(caspt2_ens)
        caspt2_df.to_csv("caspt2_ens.csv", index=False)
        subplots += 1
        dfs_to_plot.append(caspt2_df)
        titles.append("SS-CASPT2")

    if args.mscaspt2:
        mscaspt2_ens, ref_weights = get_mscaspt2_energies(text)
        mscaspt2_df = pd.DataFrame(mscaspt2_ens)
        mscaspt2_df.to_csv("mscaspt2_ens.csv", index=False)
        subplots += 1
        dfs_to_plot.append(mscaspt2_df)
        titles.append("MS-CASPT2")

    if ref_weights is not None:
        rw_df = pd.DataFrame(ref_weights)
        print("Reference weights:")
        print(rw_df.describe())

    fig, axes = plt.subplots(subplots)
    for i, df, title in zip(range(subplots), dfs_to_plot, titles):
        if subplots == 1:
            axes = [axes, ]
        df.plot(ax=axes[i], title=title)

    plt.tight_layout()
    plt.show()


    if args.nbo:
        bonds = args.nbo
        assert((len(bonds) % 2) == 0), "Bonds have to be given in pairs!"

        bond_pairs = list()
        # Sanitize the bond-indces input
        for i in range(0, len(bonds), 2):
            ind1, ind2 = bonds[i:i+2]
            assert(ind1 != ind2)
            # Force ind1 < ind2
            pair = (min(ind1, ind2), max(ind1, ind2))
            bond_pairs.append(pair)

        bo_dfs = get_natural_bond_orders(text, roots, bond_pairs)
        logging.warning("Need check if NBO converged!")

        fig, axes = plt.subplots(roots)
        for root, bo_df in zip(range(roots), bo_dfs):
            bo_df.plot(ax=axes[root], title=f"Root {root+1}")
        plt.tight_layout()
        plt.show()

    if args.cmocorr:
        check_cmocorr(text)

if __name__ == "__main__":
    run()
