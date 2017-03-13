#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import namedtuple, OrderedDict
import logging
logging.basicConfig(level=logging.INFO)
import os.path
import re
import sys

from docx import Document
from jinja2 import Environment, FileSystemLoader
import numpy as np
import simplejson as json
from tabulate import tabulate

from helper_funcs import chunks, swap_chars_in_str, print_bt_table
import rex
import td
from SpinFreeState import SpinFreeState
from Transition import Transition


THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def make_docx(output, verbose_confs_dict, irreps, fn_base):
    """Export the supplied excited states into a .docx-document."""
    docx_fn = fn_base + ".docx"
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
    trans_dict = make_trans_dict(trans)

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
            energy_rel=en_rel
        )
        sf_states.append(sfs)

    return sf_states, trans_dict


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


def flip_spins(diffs):
    """Replace ['- 2', '+ u'] with ['- 2', '+ d'] to fake closed shells."""
    flipped_diffs = diffs
    for i in range(len(diffs) - 1):
        if (diffs[i] == "- 2") and (diffs[i+1] == "+ u"):
            flipped_diffs[i+1] = "+ d"
    return flipped_diffs


def conf_diff(c1, c2):
    """Compute transitions between two configurations that are outputted
    by MOLCAS &rasscf calculations.
    E.g. 2222000 2200 and 222u000 220d would give [(3, 10)]
    """

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
    mo_pairs = ()

    # Determine occupied orbitals
    occ_indices = [i for i, mo in enumerate(c1) if (mo is "2")]

    diffs = one_by_one_diff(c1, c2)
    # If we have an open shell configuration then we
    # fake a closed shell configuration by flipping
    # spins from 'u' to 'd' in orbitals that are doubly
    # occupied in the closed shell configuration.
    #
    # Basically we replace sublists that look like this
    # ['- 2', '+ u'] with ['- 2', '+ d']
    if ("u" in c2) and not ("d" in c2):
        logging.info("Open shell configuration found!")
        diffs = flip_spins(diffs)
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

    if (len(from_mos) >= 2) or (len(to_mos) >= 2):
        logging.warning("Double excitations not supported!")
        return list()
    # Correlate information in both lists
    for from_tpl in from_mos:
        i, from_spin = from_tpl
        for to_tpl in to_mos:
            j, to_spin = to_tpl
            if to_spin is spin_pairs[from_spin]:
                #mo_pairs.append((i, j))
                mo_pairs = (i, j)
                from_mos.remove(from_tpl)
                to_mos.remove(to_tpl)
                continue

    return mo_pairs


def new_conf_diff(c1, c2):
    """Compute transitions between two configurations that are outputted
    by MOLCAS &rasscf calculations.
    E.g. 2222000 2200 and 222u000 220d would give (3, 10)
    """

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
    mo_pairs = ()

    # Determine occupied orbitals
    occ_indices = [i for i, mo in enumerate(c1) if (mo in ("u", "d", "2"))]

    print("c1", c1)
    print("c2", c2)
    diffs = np.array(one_by_one_diff(c1, c2))
    print(diffs)
    # If we have an open shell configuration then we
    # fake a closed shell configuration by flipping
    # spins from 'u' to 'd' in orbitals that are doubly
    # occupied in the closed shell configuration.
    #
    # Basically we replace sublists that look like this
    # ['- 2', '+ u'] with ['- 2', '+ d']
    if ("u" in c1) and not ("d" in c1):
        logging.info("Open shell configuration found!")
    offset = 0
    offsets = [1 if i.startswith("+") else 0 for i in diffs]
    indices = np.arange(len(diffs))
    cum_offsets = np.array(offsets, dtype=int).cumsum()
    # We have to correct the indices in diffs because for every
    # difference we introduce a new item for every difference.
    corrected_indices = indices - cum_offsets
    combined = zip(diffs, corrected_indices)
    # Only keep differences and drop items that are the same in both configs.
    diffs = [d for d in zip(diffs, corrected_indices) if d[0][0] in ("+", "-")]
    unhandled = list()
    for d_pair in chunks(diffs, 2):
        (first_occ, first_ind), (sec_occ, sec_ind) = d_pair
        pair_tpl = (first_occ[-1], first_ind, sec_occ[-1], sec_ind)
        first_tpl = (first_occ[-1], first_ind)
        sec_tpl = (sec_occ[-1], sec_ind)
        # This is always an excitation FROM this MO
        if first_occ == "- 2":
            print("excitation from")
            from_mos.append(pair_tpl)
        # This is always an excitation INTO this MO
        elif sec_occ == "+ 2":
            print("excitation to")
            to_mos.append(pair_tpl)
        # This is an excitation FROM this MO
        elif first_occ in ("- u", "- d") and sec_occ == "+ 0":
            print("excitation from")
            from_mos.append(pair_tpl)
        # This is an excitation INTO this MO
        elif first_occ == "- 0":
            print("excitation to")
            to_mos.append(pair_tpl)
        else:
            unhandled.append(pair_tpl)
            print("unhandled", pair_tpl)

    new_occ_dict = {"u": "d",
                    "d": "u"
    }
    if len(from_mos) is 1 and len(to_mos) is 1:
        print(from_mos)
        ffirst_occ, ffirst_ind, fsec_occ, fsec_ind = from_mos[0]
        tfirst_occ, tfirst_ind, tsec_occ, tsec_ind = to_mos[0]
        #return (from_mos[0][1], to_mos[0][1])
        return (fsec_ind, tsec_ind) 
    """
    # Try to match occupations
    mo_pairs = list()
    print("From_mos", from_mos)
    print("to_mos", to_mos)
    for d_pair in from_mos:
        ffirst_occ, ffirst_ind, fsec_occ, fsec_ind = d_pair
        print("checking {} -> {}".format(ffirst_occ, fsec_occ))
        if fsec_occ == "0":
            new_occ = ffirst_occ
        else:
            print("transiton from {}".format(ffirst_occ))
            new_occ = new_occ_dict[ffirst_occ]
        print("new_occ", new_occ)
        new_occs = (new_occ, "2")
        print("new_occs", new_occs)
        # If from_occ is 'u' ('d') then the newly occupied
        # orbital must have 'd' ('u') or '2' occupation. 
        # Find the target MO by chosing the MO with the right (new) occupation
        new_inds = [tind for focc, find, tocc, tind in to_mos
                    if tocc in new_occs]
        print("new_inds", new_inds)
        assert(len(new_inds) is 1)
        to_ind = new_inds[0]
        print("to_ind", to_ind)
        mo_pairs.append(from_ind, to_ind)
    return mo_pairs
    """


    """
    for i, d in enumerate(diffs):
        # Modify index i with offset because ndiff also outputs
        # chars that are in c1 and not in c2 but in the end we
        # just want whats in c2 and not in c1.
        index_in_c2 = i - offset
        previous_item = diffs[i-1]
        # Look for differences in the second configuration compared
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
            # Only increment offset once here, because + and - always come in 
            # pairs and only result in one additional list item.
            offset += 1

    spin_pairs = {
        "u": "d",
        "d": "u"
    }
    print(from_mos)
    print(to_mos)

    if (len(from_mos) is 1) and (len(to_mos) is 2):
        # Get remaining occupation of from_mos[0]
        remaining_occ = from_mos[0][-1]
        print("remaining occ", remaining_occ)
        # If the remaining occupation is 'u' ('d') then the newly occupied
        # orbital must have 'd' ('u') or '2' occupation. 
        new_occ_dict = {"u": ("d", "2"),
                        "d": ("u", "2")
        }
        # Find the target MO by chosing the MO with the right (new) occupation
        to_mos = [mo for mo in to_mos if mo[1] in new_occ_dict[remaining_occ]]
        assert(len(to_mos) is 1)
    if (len(from_mos) >= 2) or (len(to_mos) >= 2):
        logging.warning("Double excitations not supported!")
        return list()
    # Correlate information in both lists
    for from_tpl in from_mos:
        i, from_spin = from_tpl
        for to_tpl in to_mos:
            j, to_spin = to_tpl
            if to_spin is spin_pairs[from_spin]:
                #mo_pairs.append((i, j))
                mo_pairs = (i, j)
                from_mos.remove(from_tpl)
                to_mos.remove(to_tpl)
                continue
    """

    return mo_pairs


def set_mo_transitions(sf_states, trans_dict):
    # Determine ground state based on state with minimum energy
    energies_rel = np.array([sfs.energy_rel for sfs in sf_states])
    gs_index_arr = np.where(energies_rel==energies_rel.min())[0]
    gs_index = gs_index_arr[0]
    ground_state = sf_states[gs_index]
    gs_conf_line = significant_confs(sf_states[gs_index].confs)[0]
    gs_conf = gs_conf_line[2]
    gs_weight = gs_conf_line[4]
    logging.info("Found {} GS configuration ({:.1%}) in root {}.".format(
        gs_conf, gs_weight, gs_index+1))

    # Create a dict to hold verbose information about the
    # configurations for later printing
    mo_transitions = dict()
    for i, sfs in enumerate(sf_states):
        print("STATE {}".format(sfs.state))
        # Filter for configurations with CI coefficients >= 0.2 (8%)
        conf_tpls = significant_confs(sfs.confs)
        for conf_tpl in conf_tpls:
            conf_id, whoot, conf, ci, weight = conf_tpl
            mo_pair = new_conf_diff(gs_conf, conf)
            if mo_pair:
                transition = Transition(ground_state, sfs, mo_pair, -1, weight)
                sfs.add_transition(transition)


def blub():
    """
    if active_spaces:
        mo_pairs = conf_diff(gs_conf, conf)
        print(mo_pairs)
        if not mo_pairs:
            continue
        from_index, to_index = mo_pairs
        # Convert the key to a string because our dict
        # we loaded from the .json-file has string-keys
        jobiph_mos = active_spaces[sfs.jobiph]
        verbose_from = jobiph_mos[from_index]
        verbose_to = jobiph_mos[to_index]
        # Save information in verbose_confs_dict for later
        # printing.
        # Use a tuple holding the jobiph number and the root
        # as key.
        id_tpl = (jobiph, root_num)
        verbose_tpl = (verbose_from, verbose_to, weight)
        verbose_confs_dict.setdefault(sfs.state,
                                      list()).append(verbose_tpl)
    """


def create_output_list(sf_states, trans_dict, verbose_confs_dict):
    output = list()
    trs_dct = lambda root: trans_dict["1,{}".format(root)]
    for i, sfs in enumerate(sf_states):
        try:
            energy = sfs.energy
            output.append((
                sfs.state,
                sfs.jobiph,
                sfs.root,
                energy,
                hartree2eV(energy),
                hartree2nm(energy),
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

    return output


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


def make_states_list(output, irreps):
    # Sort by energy
    by_energy = sorted(output, key=lambda state: state[4])
    states = list()
    for i, state in enumerate(by_energy):
        id, jobiph, root, E, EeV, Enm, f = state
        key_tpl = (jobiph, root)
        states.append({"id": i,
                       "sym": irreps[jobiph],
                       "Enm": Enm,
                       "EeV": EeV,
                       "f": f,
                       "key": key_tpl
        })
    return states


def make_html(output, verbose_confs, irreps, fn_base, imgs):
    j2_env = Environment(loader=FileSystemLoader(THIS_DIR,
                                                 followlinks=True))
    tpl = j2_env.get_template("templates/html.tpl")
    states = make_states_list(output, irreps)
    rendered = tpl.render(states=states,
                          verbose_confs=verbose_confs,
                          imgs=imgs)
    out_fn = os.path.join(
                os.getcwd(), fn_base + ".html")
    with open(out_fn, "w") as handle:
        handle.write(rendered)


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
    object_attrs = [obj.as_list(attrs) for obj in objects]
    headers = objects[0].get_headers(attrs)
    #object_attrs = [[getattr(obj, attr) for attr in attrs] for obj in objects]
    print(tabulate(object_attrs, headers=headers,
                   tablefmt="fancy_grid", floatfmt=floatfmt))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse a &rassi-output" \
                                     " from MOLCAS.")
    parser.add_argument("fn", help="Filename of the RASSI-file.")
    parser.add_argument("--spectrum", nargs=2, type=float,
                        help="Output an spectrum.")
    parser.add_argument("--booktabs", action="store_true",
            help="Format output so that it's usable in latex-tables.")
    parser.add_argument("--sym", nargs="+", default=None,
            help="Symmetries of the WF in the JobIph-files.")
    parser.add_argument("--docx", action="store_true",
            help="Export data to a .docx-table.")
    parser.add_argument("--html", action="store_true",
            help="Export data to a .html-file.")
    args = parser.parse_args()
    fn = args.fn

    with open(fn) as handle:
        text = handle.read()

    active_spaces, imgs, irreps = load_json(fn)
    sf_states, trans_dict= parse_rassi(text)

    print_table_by_attr(sf_states, attrs=("state", "root", "mult", "energy_rel"))

    grouped_by_mult = group_sf_states_by(sf_states, "mult")
    for mult in grouped_by_mult:
        by_mult = grouped_by_mult[mult]
        set_mo_transitions(by_mult, trans_dict)
        print_table_by_attr(by_mult, attrs=("state", "root", "mult",
                                            "energy", "transitions")
        )

    fn_base = os.path.splitext(fn)[0]
    if args.spectrum:
        from_l, to_l = args.spectrum
        make_spectrum(output, from_l, to_l)
        sys.exit()
    if args.booktabs:
        print_booktabs(output, args.sym)
        sys.exit()
    if args.docx:
        make_docx(output, verbose_confs_dict, irreps, fn_base)
    if args.html:
        make_html(output, verbose_confs_dict, irreps, fn_base, imgs)

