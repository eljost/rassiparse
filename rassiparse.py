#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import OrderedDict
import logging
import os.path
import re
import sys

from docx import Document
import numpy as np
import simplejson as json

from helper_funcs import chunks, swap_chars_in_str, print_bt_table
import rex
import td


def make_docx(output, verbose_confs_dict, irreps, docx_fn):
    """Export the supplied excited states into a .docx-document."""
    docx_fn = os.path.splitext(docx_fn)[0] + ".docx"
    # The table header
    header = ("State",
              "Sym.",
              "λ / nm",
              "E / eV",
              "f",
              "natural orbitals",
              "Weight / %")
    trans_fmt = "{} → {}"
    weight_fmt = "{:.0%}"
    if irreps is None:
        # Create a fake 'irreps' dict were all jobiphs belong to
        # the totalsymmetric irrep 'A'
        irreps = {i: "A" for i in range(8)}
        logging.warning("Irrep-Hack used!")

    # Prepare the data to be inserted into the table
    as_lists = [[i, irreps[jobiph], Enm, EeV, f]
                for i, (state, jobiph, root, E, EeV, Enm, f)
                in enumerate(output[1:], 1)]
    id_tpls = [(jobiph, root)
               for state, jobiph, root, *_
               in output[1:]]
    as_fmt_lists = [[
        "S{}".format(id_),
        sym,
        "{:.1f}".format(l),
        "{:.2f}".format(dE),
        "{:.4f}".format(f)] for id_, sym, l, dE, f in as_lists]

    # For one excited state there may be several configurations
    # that contribute. This loop constructs two string, holding
    # the information about the transitions and the contributing
    # weights.
    for i, id_tpl in enumerate(id_tpls):
        try:
            verbose_confs = verbose_confs_dict[id_tpl]
        except KeyError:
            continue
        trans_list = [trans_fmt.format(from_mo, to_mo)
                      for from_mo, to_mo, _
                      in verbose_confs]
        trans_str = "\n".join(trans_list)
        weight_list = [weight_fmt.format(weight)
                       for _, __, weight in verbose_confs]
        weight_str = "\n".join(weight_list)
        as_fmt_lists[i].extend([trans_str, weight_str])

    # Prepare the document and the table
    doc = Document()
    # We need one additional row for the table header
    table = doc.add_table(rows=len(output[1:])+1,
                          cols=len(header))

    # Set header in the first row
    for item, cell in zip(header, table.rows[0].cells):
        cell.text = item

    # Start from the 2nd row (index 1) and fill in all cells
    # with the parsed data.
    for i, fmt_list in enumerate(as_fmt_lists, 1):
        for item, cell in zip(fmt_list, table.rows[i].cells):
            cell.text = item
    # Save the document
    doc.save(docx_fn)


def parse_rassi(text, reverse, swap):
    """Parse rassi file produced by MOLCAS."""
    lines = text.split("\n")

    rassi_re_tpl = (int, "\([0-9\s:/]+\)", "[2ud0\s]+", float, float)
    rassi_re, rassi_conv = rex.join_re(rassi_re_tpl)
    # List of tuples
    rassi_raw = rex.match_lines(lines, rassi_re, rassi_conv)
    # Convert to list of lists
    rassi_raw = [list(line) for line in rassi_raw]

    state_num_re = "READCI called for state"
    state_nums = rex.find_ints(text, state_num_re)

    jobiph_num_re = "This is on JobIph nr."
    jobiph_nums = rex.find_ints(text, jobiph_num_re)

    root_nums_re = "It is root nr."
    root_nums = rex.find_ints(text, root_nums_re)

    nums_combined = zip(state_nums, jobiph_nums, root_nums)

    # Split according to the configuration number
    # Iterate over all found lines holding configurations.
    # If the number of the current configuration is smaller
    # than the number of the last one a new configuration
    # is found and appended to a new root.
    root = list()
    all_roots = [root, ]
    last_conf = -1
    for line in rassi_raw:
        conf_number = line[0]
        if conf_number < last_conf:
            root = list()
            all_roots.append(root)
        root.append(line)
        last_conf = conf_number

    # Only keep unique roots.
    # This is necessary because MOLCAS prints this shit over and over.
    already_added = list()
    unique_roots = list()
    for i, id_tpl in enumerate(nums_combined):
        if id_tpl in already_added:
            continue
        unique_roots.append(all_roots[i])
        already_added.append(id_tpl)
        #print("Added root with {},{}".format(*id_tpl))

    if reverse:
        for root in unique_roots:
            for line in root:
                conf_num, trash, conf, ci, weight = line
                conf_split = conf.split()
                conf = " ".join([irrep[::-1] for irrep in conf_split])
                line[2] = conf

    for i, root in enumerate(unique_roots):
        state, jobiph, root_num = already_added[i]
        for ji, irrep, a, b in swap:
            if ji == jobiph:
                for line in root:
                    conf_num, trash, conf, ci, weight = line
                    conf_split = conf.split()
                    conf_by_irrep = conf_split[irrep]
                    swapped = swap_chars_in_str(conf_by_irrep, a, b)
                    conf_split[irrep] = swapped
                    conf = " ".join(conf_split)
                    line[2] = conf

    # Get dipole transition strengths
    trans_tpl = (int, int, float, float, float, float, float)
    trans_re, trans_conv = rex.join_re(trans_tpl)
    trans = rex.match_lines(lines, trans_re, trans_conv)

    # Get total energies
    energy_re = "::    RASSI State\s+\d+\s+Total energy:"
    energies = rex.find_floats(text, energy_re)

    # Check if RASSI changed the ordering of the states. For this
    # compare the ordering in "HAMILTON MATRIX FOR THE ORIGINAL
    # STATES" to the ordering in the "RASSI State" table and
    # swap roots accordingly
    org_ens_regex = "\s+HAMILTONIAN MATRIX FOR THE ORIGINAL STATES:.+" \
                    "\s+Diagonal, with energies\n" \
                    "(.+)" \
                    "\s+OVERLAP MATRIX FOR THE ORIGINAL STATES:"
    org_ens_match_obj = re.search(
        org_ens_regex,
        text,
        re.DOTALL
    )
    org_ens = org_ens_match_obj.groups()[0]
    org_ens = [float(n) for n in org_ens.strip().split()]
    # Find indices of the original energies in the RASSI State energies
    inds = [energies.index(oe) for oe in org_ens]
    if inds != list(range(len(inds))):
        logging.warning("Swapping of states detected!")

    # Now reorder the roots and their indices
    # trans and energies are already in the right order
    unique_roots = np.array(unique_roots)[inds]
    already_added = np.array(already_added)[inds]

    return unique_roots, already_added, trans, energies


def make_trans_dict(trans):
    """Create a dictionary holding all the oscillator
    strengths of the different electronic transitions."""
    trans_dict = {
        "1,1": 0.0,
    }
    for t in trans:
        to, ffrom, osc_str = t[:3]
        trans_dict["{},{}".format(to, ffrom)] = osc_str

    return trans_dict


def sort_by_ci_coeff(root):
    return sorted(root, key=lambda lst: -lst[-1])


def make_spectrum(transitions, start_l, end_l):
    es = [td.ExcitedState(i, 0, 0, trs[4], trs[5], trs[6], 0)
          for i, trs in enumerate(transitions[1:])]
    td.make_spectrum(es, start_l, end_l, False)


def hartree2eV(hartree):
    """Convert an energy from electronvolt to atomic units."""
    return hartree * 27.2114


def hartree2nm(hartree):
    """Convert an energy in atomic units to electronvolt."""
    with np.errstate(divide="ignore"):
        return 6.63e-34 * 2.9979e8 * 1e9 / (hartree * 27.2114 * 1.6e-19)


def significant_confs(root):
    return [conf for conf in root if conf[-1] >= 0.1]


def one_by_one_diff(str1, str2):
    """Compare two strings one by one. The output is inspired
    by difflib.ndiff(...)"""
    diffs = list()
    for i, c1 in enumerate(str1):
        c2 = str2[i]
        if c1 == c2:
            diffs.append("  {}".format(c1))
        elif c1 != c2:
            diffs.append("- {}".format(c1))
            diffs.append("+ {}".format(c2))
    return diffs


def conf_diff(c1, c2):
    """Compute transitions between two configurations that are outputted
    by MOLCAS &rasscf calculations.
    E.g. 2222000 2200 and 222u000 220d would give [(3, 10)]
    """

    assert ("u" not in c1) and ("d" not in c1), "Only singlet ground" \
            " states supported"
    """
    # Check if there could be several possible transitions and then abort
    c2_split = c2.split()
    for irrep in c2_split:
        assert not (("d" in irrep) and ("u" in irrep)), "Ambiguos transition"
    """

    # Replace whitespace
    c1 = c1.replace(" ", "")
    c2 = c2.replace(" ", "")

    # Prepare lists that hold the information about the transitions
    from_mos = list()
    to_mos = list()
    mo_pairs = list()

    # Determine occupied orbitals
    occ_indices = [i for i, mo in enumerate(c1) if (mo is "2")]

    diffs = one_by_one_diff(c1, c2)
    offset = 0
    for i, d in enumerate(diffs):
        # Modify index i with offset because ndiff also outputs
        # chars that are in c1 and not in c2 but in the end we
        # just want whats in c2 and not in c1.
        index_in_c2 = i - offset
        # Look for differences in second configuration compared
        # to the first one.
        if d.startswith("+"):
            keep = (index_in_c2, d[-1])
            # Check if difference originates from an occupied orbital
            if index_in_c2 in occ_indices:
                from_mos.append(keep)
            # Check for double excitations
            elif d.endswith("2"):
                logging.warning("Double excitations not supported!")
                return list()
            # Otherwise the MO was formerly unoccupied. Then the 
            # transition goes to this orbital.
            else:
                to_mos.append(keep)
        elif d.startswith("-"):
            offset += 1

    spin_pairs = {
        "u": "d",
        "d": "u"
    }

    #print from_mos
    #print to_mos
    #print mo_pairs
    if (len(from_mos) >= 2) or (len(to_mos) >= 2):
        logging.warning("Double excitations not supported!")
        return list()
    # Correlate information in both lists
    for from_tpl in from_mos:
        i, from_spin = from_tpl
        for to_tpl in to_mos:
            j, to_spin = to_tpl
            if to_spin is spin_pairs[from_spin]:
                mo_pairs.append((i, j))
                from_mos.remove(from_tpl)
                to_mos.remove(to_tpl)
                continue

    return mo_pairs

def run(fn, active_spaces, reverse, swap):
    with open(fn) as handle:
        text = handle.read()

    roots, root_ids, trans, energies, = parse_rassi(text, reverse, swap)

    # Create dictionary for easy lookup holding the
    # oscillator strengths with 'from,to' as keys
    trans_dict = make_trans_dict(trans)
    trs_dct = lambda root: trans_dict["1,{}".format(root)]
    # Sort roots by ci coefficient
    roots = [sort_by_ci_coeff(root) for root in roots]
    # Energies relative to the ground state
    energies_rel = np.array(energies) - energies[0]

    # Assume first root is holding the ground state configuration
    ground_state_conf = significant_confs(roots[0])[0][2]
    # Create a dict to hold verbose information about the
    # configurations for later printing
    verbose_confs_dict = dict()
    for i, root in enumerate(roots):
        state_num, jobiph, root_num = root_ids[i]
        # Filter for configurations with CI coefficients >= 0.2 (8%)
        conf_tpls = significant_confs(root)
        for conf_tpl in conf_tpls:
            conf_id, whoot, conf, ci, weight = conf_tpl
            if active_spaces:
                mo_pairs = conf_diff(ground_state_conf, conf)
                for from_index, to_index in mo_pairs:
                    # Convert the key to a string because our dict
                    # we loaded from the .json-file has string-keys
                    jobiph_mos = active_spaces[jobiph]
                    verbose_from = jobiph_mos[from_index]
                    verbose_to = jobiph_mos[to_index]
                    # Save information in verbose_confs_dict for later
                    # printing.
                    # Use a tuple holding the jobiph number and the root
                    # as key.
                    id_tpl = (jobiph, root_num)
                    verbose_tpl = (verbose_from, verbose_to, weight)
                    try:
                        verbose_confs_dict[id_tpl].append(verbose_tpl)
                    except KeyError:
                        verbose_confs_dict[id_tpl] = [verbose_tpl, ]

    # Create output table
    output = list()
    for i in range(len(energies)):
        state, jobiph, root = root_ids[i]
        try:
            output.append((
                state,
                jobiph,
                root,
                energies[i],
                hartree2eV(energies_rel[i]),
                hartree2nm(energies_rel[i]),
                trs_dct(i + 1))
            )
        except KeyError as err:
            from_state, to_state = err.args[0].split(",")
            logging.warning(
                "Oscillator strength below threshold for transition"
                " {} -> {}.".format(
                    from_state,
                    to_state))
    output = sorted(output, key=lambda row: row[3])

    return output, verbose_confs_dict


def print_booktabs(output, symms):
    bt_list = list()
    # Drop HF-configuration
    for line in output:
        state, jobiph, root, E, EeV, Enm, f = line
        if symms:
            sym = symms[jobiph-1]
        else:
            sym = jobiph
        bt_list.append((sym, EeV, Enm, f))

    print_bt_table(bt_list)


def print_output(output, verbose_confs_dict, header):
    header_tpl = ["{}" for item in output_headers]
    header_str = "\t".join(header_tpl)
    print(header_str.format(*output_headers))
    for line in output:
        state, jobiph, root, E, EeV, Enm, f = line
        id_tpl = (jobiph, root)
        # with E
        #out_tpl =("{}", "{}", "{:.5f}", "{:.2f}", "{:.1f}", "{:.5f}")
        # without E
        out_str ="{} {}\t\t{:.2f}\t\t{:.1f}\t\t{:.5f}"
        print(out_str.format(
            jobiph,
            root,
            #E,
            EeV,
            Enm,
            f
        ))
        try:
            verbose_confs = verbose_confs_dict[id_tpl]
            verbose_tpl = ("\t{}", "->",  "{}", "\t{:.2%}")
            verbose_str = "\t".join(verbose_tpl)
            for vc in verbose_confs:
                from_mo, to_mo, weight = vc
                print(verbose_str.format(
                    from_mo,
                    to_mo,
                    weight
                ))
        except KeyError:
            pass
        print


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse a &rassi-output" \
                                     " from MOLCAS.")
    parser.add_argument("fn", help="Filename of the RASSI-file.")
    parser.add_argument("--rev", action="store_true",
                        help="Reverse the string holding the configuration.")
    parser.add_argument("--swap", type=int, nargs="+",
                        help="""Indices to swap in the f***ing scrambled 
                             configuration line. Expecting groups of 4
                             integers. The first determines to JobIph to
                             be used for the swapping. The next determines
                             the irrep. The last two integers are the
                             indices to be swapped in the configuration string
                             of the irrep. The numbering of the last three
                             integers starts at 0.
                             --swap 2 1 0 6
                             swaps the first char with the seventh char
                             in the irrep 1 in the 2nd Jobfile.""")
    parser.add_argument("--spectrum", nargs=2, type=float,
                        help="Output an spectrum.")
    parser.add_argument("--booktabs", action="store_true",
            help="Format output so that it's usable in latex-tables.")
    parser.add_argument("--sym", nargs="+", default=None,
            help="Symmetries of the WF in the JobIph-files.")
    parser.add_argument("--docx", action="store_true",
            help="Export data to a .docx-table.")
    args = parser.parse_args()
    fn = args.fn
    rev = args.rev

    if args.swap:
        swap = list(chunks(args.swap, 4))
    else:
        swap = list()

    try:
        verbose_fn = os.path.splitext(fn)[0] + ".json"
        with open(verbose_fn) as handle:
            json_data = json.load(handle, object_pairs_hook=OrderedDict)
            active_spaces = {int(key): json_data[key]["as"] for key in json_data}
            irreps = {int(key): json_data[key]["irrep"] for key in json_data}
    except IOError:
        active_spaces = None
        irreps = None

    output, verbose_confs_dict = run(fn, active_spaces, rev, swap)

    output_headers = ("State", "JobIph", "Root", "E in a.u.",
        "dErel in eV", "Erel in nm", "f")
    output_headers = ("JobIph", "Root", "dE in eV", "dE in nm", "f")

    if args.spectrum:
        from_l, to_l = args.spectrum
        make_spectrum(output, from_l, to_l)
        sys.exit()
    if args.booktabs:
        print_booktabs(output, args.sym)
        sys.exit()
    if args.docx:
        make_docx(output, verbose_confs_dict, irreps, fn)

    print_output(output, verbose_confs_dict, output_headers)
