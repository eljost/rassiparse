#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import namedtuple, OrderedDict
import configparser
import logging
import os
import re
import sys

from docx import Document
from jinja2 import Environment, FileSystemLoader
import numpy as np
import simplejson as json
from tabulate import tabulate

import rex
from SpinFreeState import SpinFreeState
from ConfDiff import ConfDiff


THIS_DIR = os.path.dirname(os.path.realpath(__file__))
CONFIG = configparser.ConfigParser()
# Try to load a customized config from THIS_DIR
if not CONFIG.read(os.path.join(THIS_DIR, "config.ini")):
    logging.warning("Used fallback config 'templates/config.tpl'")
    CONFIG.read(os.path.join(THIS_DIR, "templates/config.tpl"))


def make_docx(sf_states, attrs, fn_base):
    """Export the supplied excited states into a .docx-document."""
    docx_fn = fn_base + ".docx"

    sfs_attrs = [sfs.as_str_list(attrs, newlines=True)
                 for sfs in sf_states]
    headers = sf_states[0].get_headers(attrs)

    # Prepare the document and the table
    doc = Document()
    # We need one additional row for the table header
    table = doc.add_table(rows=len(sfs_attrs)+1,
                          cols=len(headers))

    # Set header in the first row
    for item, cell in zip(headers, table.rows[0].cells):
        cell.text = item

    # Start from the 2nd row (index 1) and fill in all cells
    # with the parsed data.
    for i, attrs_for_row in enumerate(sfs_attrs, 1):
        for item, cell in zip(attrs_for_row, table.rows[i].cells):
            cell.text = str(item)
    # Save the document
    doc.save(docx_fn)


def make_html(sf_states, fn_base):
    j2_env = Environment(loader=FileSystemLoader(
                                    os.path.join(THIS_DIR, "templates"),
                                    followlinks=True))
    trans_tpl = j2_env.get_template("rassi_trans.tpl")

    trans_rendered = trans_tpl.render(sf_states=sf_states)
    trans_fn = os.path.join(
                os.getcwd(), fn_base + ".trans.html")
    with open(trans_fn, "w") as handle:
        handle.write(trans_rendered)

    single_tpl = j2_env.get_template("rassi_single.tpl")

    single_rendered = single_tpl.render(sf_states=sf_states)
    single_fn = os.path.join(
                os.getcwd(), fn_base + ".single.html")
    with open(single_fn, "w") as handle:
        handle.write(single_rendered)


def make_export(sf_states, fn_base):
    attrs = "state sym dE_global dE_global_eV osc".split()
    as_lists = [sfs.as_list(attrs) for sfs in sf_states]
    header = sf_states[0].get_headers(attrs)
    # Replace whitespace
    header = " ".join([re.sub("\s", "", h) for h in header])
    header = re.sub("Δ", "d", header)
    as_array = np.array(as_lists)
    out_fn = "{}.dat".format(fn_base)
    np.savetxt(out_fn, as_array, header=header)


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


def parse_state(text, irrep_inds):
    state_re = "state\s*(\d+)"
    jobiph_re = "JobIph nr.\s*(\d+)"
    root_re = "It is root nr.\s*(\d+)"
    sym_re = "Its symmetry  =\s*(\d+)"
    mult_re = "Spin multiplic=\s*(\d+)"
    confs_re_raw = (int, "\([0-9\s:/]+\)", "[2ud0\s]+?", float, float)
    conf_re, conf_conv = rex.join_re(confs_re_raw)

    def find(regex, txt):
        return int(re.search(regex, txt).groups()[0])

    state = find(state_re, text)
    jobiph = find(jobiph_re, text)
    root = find(root_re, text)
    sym = find(sym_re, text)
    mult = find(mult_re, text)
    confs_raw = re.findall(conf_re, text)
    confs = [[cv(cf)
             for cv, cf
             in zip(conf_conv, conf)] for conf in confs_raw]
    for conf in confs:
        old_conf_string = conf[2]
        old_conf_list = old_conf_string.split()
        new_conf_list = list()
        for irrep, inds in enumerate(irrep_inds, 1):
            new_conf_list.append("".join([old_conf_list[ind]
                                 for ind in inds])
            )
        new_conf_string = " ".join(new_conf_list)
        conf[2] = new_conf_string
        logging.debug("Reordered configuration {} to {}.".format(
                      old_conf_string, new_conf_string)
        )

    return (state, jobiph, root, sym, mult, confs)


def parse_subspaces(text, irreps):
    # When there is no symmetry (only 1 irrep) we have only
    # one section in the configuration string, even if all
    # subspaces (RAS1-RAS3) are present.
    if irreps == 1:
        return [[0]]

    # Construct regex based on number of irreps
    ras_base_re = ("RAS{}", "(\d+)") + ("(\d+)", )*irreps
    subspaces = (1, 2, 3)
    curr_ind = 0
    # Initialize a nested list
    irrep_inds = list()
    for irrep in range(irreps):
        irrep_inds.append([])

    for subspace in subspaces:
        ras_re = "\s*".join(ras_base_re).format(subspace)
        mobj = re.search(ras_re, text)
        as_ints = [int(i) for i in mobj.groups()]
        # The first item holds the sum of the following entries
        # and right now we dont need it.
        for irrep, per_irrep in enumerate(as_ints[1:], 1):
            # In this case there are no MOs for this subspace and irrep
            if per_irrep == 0:
                continue
            # When there are MOs in this subspace for this irrep
            # keep the index in the configuration string.
            #
            # Irreps: 1 and 2
            # RAS1           2       1   1
            # RAS2           6       4   2
            # RAS3           2       1   1
            # Configuration: 2 2 2200 20 0 0
            #
            # We have two irreps and there are MOs for both irreps in every
            # subspace. Splitting the configuration string would lead to
            # [2, 2, 2200, 20, 0, 0]. This equals the ordering
            # [RAS1(1), RAS1(2), RAS2(1), RAS2(2), RAS3(1), RAS3(2)]. We keep
            # the indices that belong to a specific irrep se we can reorder
            # the configuration string later on.
            # The example above would lead to the following irrep_inds_dict:
            # {1: [0, 2, 4], 2: [1, 3, 5]}. Entries 0, 2 and 4 belong to the
            # first irrep. Entries 1, 3, 5 belong to the second irrep.
            irrep_inds[irrep-1].append(curr_ind)
            curr_ind += 1

    return irrep_inds


def parse_rassi(text):
    nr_of_irreps_re = "NR of irreps:"
    nr_of_irreps_list = rex.find_ints(text, nr_of_irreps_re)
    assert(len(nr_of_irreps_list) == 1)
    irreps = nr_of_irreps_list[0]
    irrep_inds = parse_subspaces(text, irreps)
    logging.debug("Indices per irrep in the configuration "
                  "string: {}".format(irrep_inds)
    )

    states_regex = "READCI called for(.+?)\*"
    raw_states = re.findall(states_regex, text, re.DOTALL)
    parsed_states = [parse_state(raw_state, irrep_inds)
                     for raw_state in raw_states]
    # Only keep unique states
    unique_state_ids = list()
    states = list()
    for pst in parsed_states:
        state = pst[0]
        if state in unique_state_ids:
            continue
        states.append(pst)
        unique_state_ids.append(state)
    logging.info("Found {} states.".format(len(states)))

    lines = text.split("\n")

    # Get dipole transition strengths but only use the spin-free
    # section.
    sf_trans_re = "Dipole transition strengths:(.+?)" + "-"*20 + "\n--"
    sf_trans_mobj = re.search(sf_trans_re, text, re.DOTALL)
    sf_trans_lines = sf_trans_mobj.groups()[0].split("\n")
    trans_tpl = (int, int, float, float, float, float, float)
    trans_re, trans_conv = rex.join_re(trans_tpl)
    trans = rex.match_lines(sf_trans_lines, trans_re, trans_conv)
    from_to_set = set(tuple(t[:2]) for t in trans)
    # Assure that we have every transition only once, e.g. we didn't
    # accidently parsed also SO states.
    assert(len(trans) == len(from_to_set))

    # Create dictionary for easy lookup holding the
    # oscillator strengths with 'from,to' as keys
    trans_dict = dict()
    for t in trans:
        to, from_, osc = t[:3]
        trans_dict[(to, from_)] = osc

    # Get total energies
    energy_re = "::    RASSI State\s+\d+\s+Total energy:"
    energies = rex.find_floats(text, energy_re)

    energy_shifts = rex.find_floats(text, "\(Shifted by EVAC \(a\.u\.\) =")

    # Check if RASSI changed the ordering of the states. For this
    # compare the ordering in "HAMILTON MATRIX FOR THE ORIGINAL
    # STATES" to the ordering in the "RASSI State" table and
    # swap roots accordingly.
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
    if not np.all(np.isclose(energies, org_ens)):
        logging.warning("Swapping of states detected!")
        sys.exit()

    energies = np.array(energies)
    gs_energy = energies.min()
    gs_index_arr = np.where(energies==gs_energy)[0]
    assert (gs_index_arr.size == 1), "Degenerate ground state found!"
    energies_rel = np.array(energies) - gs_energy

    sf_states = list()
    for pfs, en, en_rel in zip(states, energies, energies_rel):
        state, jobiph, root, sym, mult, confs = pfs
        sfs = SpinFreeState(
            state=state,
            jobiph=jobiph,
            root=root,
            sym=sym,
            mult=mult,
            confs=confs,
            energy=en,
            dE_global=en_rel
        )
        sf_states.append(sfs)

    return sf_states, trans_dict


def sort_by_ci_coeff(root):
    return sorted(root, key=lambda lst: -lst[-1])


def significant_confs(root, thresh=0.1):
    return [conf for conf in root if conf[-1] >= thresh]


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


def get_electron_spins(occ_1, occ_2):
    """Takes two occupation labels ("0", "u", "d", "2") as argument
    and returns the spin(s) of the electron(s) that can produce this
    transition."""
    occ_set = set((occ_1, occ_2))
    ud_tpl = (set(("0", "2")), set(("u", "d")))
    u_tpl = (set(("0", "u")), set(("2", "d")))
    d_tpl = (set(("0", "d")), set(("2", "u")))
    # These occupation changes correspond to transitions of two
    # electrons with up and down spin.
    if occ_set in ud_tpl:
        return ("u", "d")
    elif occ_set in u_tpl:
        return ("u", )
    elif occ_set in d_tpl:
        return ("d", )
    else:
        sys.exit("This should never happen ;)")


def conf_diff(c1, c2):
    """Compute transitions between two configurations that are outputted
    by MOLCAS &rasscf calculations.
    E.g. 2222000 2200 and 222u000 220d would return ((3, 10), )
    """

    new_occ_dict = {"u": "d", "d": "u", "0": "u"}

    # Delete whitespace
    c1 = c1.replace(" ", "")
    c2 = c2.replace(" ", "")

    # Prepare lists that hold the information about the transitions
    from_mos = list()
    to_mos = list()
    mo_pairs = list()

    logging.debug("Got c1: {}".format(c1))
    logging.debug("Got c2: {}".format(c2))
    diffs = np.array(one_by_one_diff(c1, c2))
    if ("u" in c1) and not ("d" in c1):
        logging.debug("Open shell configuration found!")
    # We have to correct the indices in diffs because for every
    # difference we introduce a new item. So we can't use the
    # indices of the items in diffs directly.
    offsets = [1 if i.startswith("+") else 0 for i in diffs]
    indices = np.arange(len(diffs))
    cum_offsets = np.array(offsets, dtype=int).cumsum()
    corrected_indices = indices - cum_offsets
    combined = zip(diffs, corrected_indices)
    # Only keep differences and drop items that are the same in both configs.
    diffs = [d for d in zip(diffs, corrected_indices) if d[0][0] in ("+", "-")]
    for d_pair in chunks(diffs, 2):
        # The '+' and '-' could be dropped now because the first item is alwayas
        # '-' and the second is always '+' but we keep them for now because the
        # code should be easier to follow this way.
        # der zweite index in d_pair ist auch sinnlos weil er immer
        # gleich ist
        logging.debug("Checking diff-pair {}.".format(d_pair))
        # Don't consider the second because it always equals the first index
        (first_occ, ind), (sec_occ, _) = d_pair
        assert(ind == _)
        # This tuple describes the electron spin in this transition
        el_spins = get_electron_spins(first_occ[-1], sec_occ[-1])
        trans_list = [(el_spin, ind) for el_spin in el_spins]
        # This is always an excitation FROM this MO
        if first_occ == "- 2":
            from_mos.extend(trans_list)
        # This is always an excitation INTO this MO
        elif sec_occ == "+ 2":
            to_mos.extend(trans_list)
        # This is always an excitation FROM this MO
        elif first_occ in ("- u", "- d") and sec_occ == "+ 0":
            from_mos.extend(trans_list)
        # This is always an excitation INTO this MO
        elif first_occ == "- 0" and sec_occ in ("+ u", "+ d"):
            to_mos.extend(trans_list)
        # This handles cases where the occupation flips from
        # (u -> d) or (d -> u). This equals the transition of
        # two electrons.
        else:
            from_tpl = (first_occ[-1], ind)
            to_tpl = (sec_occ[-1], ind)
            from_mos.append(from_tpl)
            to_mos.append(to_tpl)

    logging.debug("Transitions from MOs: {}".format(from_mos))
    logging.debug("Transitions to MOs: {}".format( to_mos))
    handle_with_spin_flip = list()
    for occ_mo in from_mos:
        occ_spin, occ_ind = occ_mo
        logging.debug("Handling transition from "
                      "MO {} with spin {}.".format(occ_ind, occ_spin)
        )
        virt_mos = [(virt_spin, virt_ind) for virt_spin, virt_ind in to_mos
                     if (occ_spin == virt_spin)]
        if len(virt_mos) == 1:
            mo_pairs.append((occ_ind, virt_mos[0][1]))
            to_mos.remove(virt_mos[0])
        elif (len(from_mos) == 1) and (len(to_mos) == 1):
            logging.warning("Spin flip transition!")
            mo_pairs.append((occ_ind, to_mos[0][1]))
        else:
            # Don't return anything because this state can't be described
            # completly without considering all MO transitions.
            logging.warning("Ambiguous transitions. Can't handle this.")
            return list()

    return mo_pairs


def set_mo_transitions(sf_states, trans_dict, gs_conf):
    for sfs in sf_states:
        logging.info("Checking state {} from JobIph {}".format(
            sfs.state, sfs.state)
        )
        # Filter for configurations with CI coefficients >= 0.2 (8%)
        conf_tpls = significant_confs(sfs.confs)
        for conf_tpl in conf_tpls:
            conf_id, whoot, conf, ci, weight = conf_tpl
            logging.info("Checking configuration {} with "
                         "weight {:.1%}".format(conf, weight)
            )
            mo_pairs = conf_diff(gs_conf, conf)
            cd = ConfDiff(conf, mo_pairs, weight)
            sfs.add_confdiff(cd)


def set_single_mos(sfs):
    single_mos = list()
    for conf_line in significant_confs(sfs.confs):
        conf = conf_line[2].replace(" ", "")
        images_spins = [(sfs.mo_images[i], spin) for i, spin in enumerate(conf)
                         if spin in ("u", "d")]
        weight = conf_line[-1]
        if images_spins:
            images, spins = zip(*images_spins)
            smo = (weight, images, spins)
        else:
            # This applies to close shell configuration with no singly
            # occupied MOs.
            smo = (weight, [], [])
        single_mos.append(smo)
    sfs.single_mos = single_mos


def get_img_nums(path):
    png_fns = [f for f in os.listdir(path) if f.endswith(".png")]
    mo_mobjs = [re.match("mo_(\d+).+\.png", png) for png in png_fns]
    mo_nums = [int(mobj.groups()[0]) for mobj in mo_mobjs]
    mo_nums = sorted(list(set(mo_nums)))
    return mo_nums


def group_sf_states_by(sf_states, attr):
    attrs = np.array([getattr(sfs, attr) for sfs in sf_states])
    unique_attrs = set(attrs)
    grouped_dict = dict()
    for ua in unique_attrs:
        indices = np.where(attrs==ua)[0]
        grouped_dict[ua] = [sf_states[i] for i in indices]

    return grouped_dict


def print_table_by_attr(objects, attrs, floatfmt=".4f"):
    sfs_attrs = [obj.as_list(attrs) for obj in objects]
    headers = objects[0].get_headers(attrs)
    print(tabulate(sfs_attrs, headers=headers,
                   tablefmt="fancy_grid", floatfmt=floatfmt))


def handle_rassi(sf_states, trans_dict, ground_state=None):
    if not ground_state:
        ground_state = sorted(sf_states, key=lambda sfs: sfs.energy)[0]
    # Determine ground state configuration based on state with
    # minimum energy. Sort the significant confs by weight.
    sig_confs= sorted(significant_confs(ground_state.confs),
                                        key=lambda c: -c[-1])
    gs_conf_line = sig_confs[0]
    gs_conf = gs_conf_line[2]
    gs_weight = gs_conf_line[-1]
    logging.info("Found {} GS configuration ({:.1%}) in state {}.".format(
        gs_conf, gs_weight, ground_state.state))
    for sfs in sf_states:
        trans_tpl = (ground_state.state, sfs.state)
        try:
            osc = trans_dict[trans_tpl]
        except KeyError:
            osc = 0
            logging.warning("Oscillator strength below threshold for "
                            "transition {} -> {}!".format(
                                *trans_tpl)
            )
        sfs.set_ground_state(ground_state, osc)

    set_mo_transitions(sf_states, trans_dict, gs_conf)
    # Sort by energy difference to the ground state
    sf_states = sorted(sf_states, key=lambda sfs: sfs.dE_gs_eV)
    # Add a new numbering according to the energy of the states.
    # The ground state is 0.
    [setattr(sfs, "state_rel", i) for i, sfs in enumerate(sf_states)]


def load_mo_images(path):
    png_fns = [png for png in os.listdir(path)
               if png.endswith(".png")]
    # Filter for MO imgs with matching names
    mo_re = "mo_(\d+)(?:\.job)?(\d+)?\.png"
    img_dict = dict()
    for png in png_fns:
        mobj = re.match(mo_re, png)
        if mobj:
            mo = int(mobj.groups()[0])
            job = mobj.groups()[1]
            # We expect a rassi with only one jobiph
            if not job:
                job = 1
            else:
                job = int(job)
            img_dict.setdefault(job, list()).append((mo, png))

    for job in img_dict:
        img_dict[job] = sorted(img_dict[job], key=lambda tpl: tpl[0])

    return img_dict


def set_mo_images(sf_states, image_dict):
    for sfs in sf_states:
        images = image_dict[sfs.jobiph]
        sfs.set_mo_images(images)


def load_mo_names(fn_base):
    json_fn = fn_base + ".json"
    try:
        with open(json_fn) as handle:
            mo_names = json.load(handle)
    except FileNotFoundError:
        logging.warning("Couldn't find {}.".format(json_fn))
        mo_names = dict()
    return mo_names


def set_mo_names(sf_states, mo_names_dict):
    for sfs in sf_states:
        names = mo_names_dict[str(sfs.jobiph)]
        sfs.set_mo_names(names)


def parse_args(args):
    parser = argparse.ArgumentParser(
                        description="Parse a &rassi-output from MOLCAS.")
    parser.add_argument("fn", help="Filename of the RASSI-output.")
    parser.add_argument("--docx", action="store_true",
            help="Export data to DOCX.")
    parser.add_argument("--html", action="store_true",
            help="Export data to HTML.")
    parser.add_argument("--gs", type=int, default=None,
            help="Set a global ground state configuration that will be used "
                 "to determine differences between configurations. Expects a "
                 "RASSI state (1..number of states). If --gs=0 the ground "
                 "state configuration will be determined from the "
                 "configuration with the highest weight in the ground state.")
    parser.add_argument("--attrs", default="default",
            help="Load a different set of attributes from the config file.")
    parser.add_argument("--info", action="store_true",
            help="Print more information.")
    parser.add_argument("--debug", action="store_true",
            help="Print even more information.")

    return parser.parse_args(args)


def load_attrs(section):
    sf_states_attrs = CONFIG[section]["all_states"].strip().split()
    by_mult_attrs = CONFIG[section]["by_mult"].strip().split()
    docx_attrs = CONFIG[section]["docx"].strip().split()

    return sf_states_attrs, by_mult_attrs, docx_attrs


def handle_spin_free_states(sf_states, trans_dict, mo_names_dict,
                            ground_state=None):
    # Try to load the MO images
    image_dict = load_mo_images(".")

    grouped_by_mult = group_sf_states_by(sf_states, "mult")
    for mult in grouped_by_mult:
        by_mult = grouped_by_mult[mult]
        handle_rassi(grouped_by_mult[mult], trans_dict, ground_state)

        # Precalculate JOB00N string.
        jobiphs = set([sfs.jobiph for sfs in by_mult])
        jobiph_strings = ["JOB{:0>3}".format(j) for j in jobiphs]
        try:
            set_mo_images(by_mult, image_dict)
            [set_single_mos(sfs) for sfs in by_mult]
        except KeyError:
            logging.warning("Couldn't find MO images for "
                            "states from {}".format(jobiph_strings))
        try:
            set_mo_names(by_mult, mo_names_dict)
        except KeyError:
            logging.warning("Couldn't find MO names for "
                            "states from {}".format(jobiph_strings))
    return grouped_by_mult


def run():
    args = parse_args(sys.argv[1:])

    sf_states_attrs, by_mult_attrs, docx_attrs = load_attrs(args.attrs)

    if args.info:
        logging.basicConfig(level=logging.INFO)
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    with open(args.fn) as handle:
        text = handle.read()

    fn_base = os.path.splitext(args.fn)[0]
    fn_base_fmt = "{}.mult{}"

    mo_names_dict = load_mo_names(fn_base)
    sf_states, trans_dict = parse_rassi(text)
    
    if args.gs:
        ground_state = [sfs for sfs in sf_states
                        if sfs.state == args.gs][0]
    else:
        ground_state = None

    grouped_by_mult = handle_spin_free_states(sf_states,
                                              trans_dict,
                                              mo_names_dict,
                                              ground_state
    )
                    
    # Sort by energy difference to the global energy minimum
    sf_states_sorted = sorted(sf_states, key=lambda sfs: sfs.dE_global_eV)
    print_table_by_attr(sf_states_sorted, sf_states_attrs)
    for mult in grouped_by_mult:
        fn_base_mult = fn_base_fmt.format(fn_base, mult)
        # Hmm, this should be handled in a different way
        # Sort by energy difference to the ground state of this multiplicity
        by_mult = sorted(grouped_by_mult[mult], key=lambda sfs: sfs.dE_gs_eV)
        print_table_by_attr(by_mult, by_mult_attrs)
        if args.html:
            make_html(by_mult, fn_base_mult)
        if args.docx:
            make_docx(by_mult, docx_attrs, fn_base_mult)
        make_export(by_mult, fn_base_mult)


if __name__ == "__main__":
    run()
