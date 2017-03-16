#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import namedtuple, OrderedDict
import logging
#logging.basicConfig(level=logging.INFO)
logging.basicConfig(level=logging.DEBUG)
import os
import re
import sys

from docx import Document
from jinja2 import Environment, FileSystemLoader
import numpy as np
import simplejson as json
from tabulate import tabulate

from helper_funcs import chunks
import rex
from SpinFreeState import SpinFreeState
from ConfDiff import ConfDiff


THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def make_docx(sf_states, attrs, fn_base):
    """Export the supplied excited states into a .docx-document."""
    docx_fn = fn_base + ".docx"
    # The table header
    trans_fmt = "{} → {}"
    weight_fmt = "{:.0%}"

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
    j2_env = Environment(loader=FileSystemLoader(THIS_DIR,
                                                 followlinks=True))
    tpl = j2_env.get_template("templates/html.tpl")

    rendered = tpl.render(sf_states=sf_states)
    out_fn = os.path.join(
                os.getcwd(), fn_base + ".html")
    with open(out_fn, "w") as handle:
        handle.write(rendered)


def parse_state(text):
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

    return (state, jobiph, root, sym, mult, confs)


def parse_rassi(text):
    states_regex = "READCI called for(.+?)\*"
    raw_states = re.findall(states_regex, text, re.DOTALL)
    parsed_states = [parse_state(raw_state) for raw_state
                     in raw_states]
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
    # Get dipole transition strengths
    trans_tpl = (int, int, float, float, float, float, float)
    trans_re, trans_conv = rex.join_re(trans_tpl)
    trans = rex.match_lines(lines, trans_re, trans_conv)
    # Create dictionary for easy lookup holding the
    # oscillator strengths with 'from,to' as keys
    trans_dict = dict()
    for t in trans:
        to, from_, osc = t[:3]
        trans_dict[(to, from_)] = osc

    # Get total energies
    energy_re = "::    RASSI State\s+\d+\s+Total energy:"
    energies = rex.find_floats(text, energy_re)

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
    # Find indices of the original energies in the RASSI State energies
    inds = [energies.index(oe) for oe in org_ens]
    if inds != list(range(len(inds))):
        logging.warning("Swapping of states detected!")
        sys.exit()

    """
    # Now reorder the roots and their indices
    # trans and energies are already in the right order
    confs = np.array(confs)[inds]
    id_tpls = np.array(id_tpls)[inds]
    """

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
    spin_flip_mos = list()
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
            spin_flip_mos.extend(trans_list)

    logging.debug("Transitions from MOs: {}".format(from_mos))
    logging.debug("Transitions to MOs: {}".format( to_mos))
    handle_with_spin_flip = list()
    for occ_mo in from_mos:
        occ_spin, occ_ind = occ_mo
        logging.debug("Handling transition from "
                      "MO {} with spin {}.".format(occ_ind, occ_spin)
        )
        virt_inds = [(virt_spin, virt_ind) for virt_spin, virt_ind in to_mos
                     if (occ_spin == virt_spin)]
        if len(virt_inds) == 1:
            mo_pairs.append((occ_ind, virt_inds[0][1]))
            to_mos.remove(virt_inds[0])
        else:
            # Don't return anything because this state can't be described
            # completly without considering all MO transitions.
            logging.warning("Ambiguous transitions. Can't handle this.")
            return []

    return mo_pairs


def set_mo_transitions(sf_states, trans_dict):
    sf_states = sorted(sf_states, key=lambda sfs: sfs.energy)
    # Determine ground state based on state with minimum energy
    ground_state = sf_states[0]
    gs_conf_line = significant_confs(ground_state.confs)[0]
    gs_conf = gs_conf_line[2]
    gs_weight = gs_conf_line[4]
    logging.info("Found {} GS configuration ({:.1%}) in state {}.".format(
        gs_conf, gs_weight, ground_state.state))

    for sfs in sf_states:
        trans_tpl = (ground_state.state, sfs.state)
        try:
            osc = trans_dict[trans_tpl]
        except KeyError:
            osc = 0
            logging.warning("Oscillator strength below threshold for "
                            "transition between states {} -> {}!".format(
                                *trans_tpl)
            )
        sfs.set_ground_state(ground_state, osc)
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
            cd = ConfDiff(mo_pairs, weight)
            sfs.add_confdiff(cd)


def get_img_nums(path):
    png_fns = [f for f in os.listdir(path) if f.endswith(".png")]
    mo_mobjs = [re.match("mo_(\d+).+\.png", png) for png in png_fns]
    mo_nums = [int(mobj.groups()[0]) for mobj in mo_mobjs]
    mo_nums = sorted(list(set(mo_nums)))
    return mo_nums


def load_json(fn):
    verbose_fn = os.path.splitext(fn)[0] + ".json"
    try:
        with open(verbose_fn) as handle:
            json_data = json.load(handle, object_pairs_hook=OrderedDict)
    except IOError:
        active_spaces = None
        # When no .json file is defined try to look for pictures in the
        # current directory and try to use them.
        imgs = get_img_nums(".")
        # Create a fake 'irreps' dict were all jobiphs belong to
        # the totalsymmetric irrep 'A'
        irreps = {i: "A" for i in range(8)}
        logging.warning("Irrep-Hack used!")
        return active_spaces, imgs, irreps

    active_spaces = {int(key): json_data[key]["as"] for key in json_data}
    irreps = {int(key): json_data[key]["irrep"] for key in json_data}
    imgs = dict()
    for key in json_data:
        for_jobiph = json_data[key]
        if len(irreps) > 1:
            img_fns = ["mo_{}.irrep{}.png".format(img, key)
                       for img in for_jobiph["imgs"]]
        else:
            img_fns = ["mo_{}.png".format(img) for img in for_jobiph["imgs"]]
        as_ = for_jobiph["as"]
        img_dict = {mo: img_fn for mo, img_fn in zip(as_, img_fns)}
        imgs[int(key)] = img_dict

    return active_spaces, imgs, irreps


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


def load_mo_images(path):
    png_fns = [png for png in os.listdir(path)
               if png.endswith(".png")]
    # Filter for MO imgs with matching names
    mo_re = "mo_(\d+)(?:\.job)(\d+)?\.png"
    img_dict = dict()
    for png in png_fns:
        mobj = re.match(mo_re, png)
        if mobj:
            mo = int(mobj.groups()[0])
            job = int(mobj.groups()[1])
            img_dict.setdefault(job, list()).append((mo, png))
    for job in img_dict:
        img_dict[job] = sorted(img_dict[job], key=lambda tpl: tpl[0])

    return img_dict


def set_images(sf_states, image_dict):
    for sfs in sf_states:
        images_for_jobiph = image_dict[sfs.jobiph]
        sfs.set_images(images_for_jobiph)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                        description="Parse a &rassi-output from MOLCAS.")
    parser.add_argument("fn", help="Filename of the RASSI-output.")
    parser.add_argument("--sym", nargs="+", default=None,
            help="Symmetries of the WF in the JobIph-files.")
    parser.add_argument("--docx", action="store_true",
            help="Export data to a .docx-table.")
    parser.add_argument("--html", action="store_true",
            help="Export data to a .html-file.")
    parser.add_argument("--gsroot", nargs="+")
    args = parser.parse_args()
    fn = args.fn

    # Try to load the MO images
    image_dict = load_mo_images(".")

    with open(fn) as handle:
        text = handle.read()
    fn_base = os.path.splitext(fn)[0]

    sf_states, trans_dict = parse_rassi(text)
    
    print_table_by_attr(sf_states, attrs=("state", "root", "mult",
                                          "dE_global_eV"))

    grouped_by_mult = group_sf_states_by(sf_states, "mult")
    for mult in grouped_by_mult:
        by_mult = grouped_by_mult[mult]
        set_mo_transitions(by_mult, trans_dict)
        by_mult = sorted(by_mult, key=lambda sfs: sfs.dE_gs_eV)
        # Add a new numbering according to the energy of the states.
        # The ground state is 0.
        [setattr(sfs, "state_rel", i) for i, sfs in enumerate(by_mult)]
        try:
            set_images(by_mult, image_dict)
        except KeyError:
            jobiphs = set([sfs.jobiph for sfs in by_mult])
            jobiph_strings = ["JOB{:0>3}".format(j) for j in jobiphs]
            logging.warning("Couldn't find MO images for "
                            "states {}".format(jobiph_strings))
        print_table_by_attr(by_mult, attrs=("state_rel", "state", "root",
                                            "mult", "sym", "dE_gs_eV",
                                            "dE_gs_nm", "osc", "confdiffsw",)
        )
        fn_base_mult = "{}.mult{}".format(fn_base, mult)
        if args.html:
            make_html(by_mult, fn_base_mult)
        if args.docx:
            docx_attrs = ("state_rel", "sym", "dE_gs_nm",
                          "dE_gs_eV", "osc", "confdiffs",
                          "weights")
            make_docx(by_mult, docx_attrs, fn_base_mult)
