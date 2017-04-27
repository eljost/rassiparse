#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
logging.basicConfig(level=logging.INFO)
import os
import re
import sys

from jinja2 import Environment, FileSystemLoader
from matplotlib import rc
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

from rassiparse import (significant_confs, conf_diff,
                        parse_rassi, handle_spin_free_states,
                        load_mo_names)
from SpinOrbitState import SpinOrbitState


TEX = False

def get_block_lines(text, regex):
    mobj = re.search(regex, text, re.DOTALL)
    block = mobj.groups()[0].strip()
    as_lists = [line.strip().split() for line in block.split("\n")]
    return as_lists


def normalize_energies(energies):
    energies = np.array(energies, dtype=np.float)
    energies -= energies.min()
    return energies


def parse_sorassi(text):
    # Only consinder the spin-orbit section
    so_section_re = "Spin-orbit section(.+)"
    so_section = re.search(so_section_re, text, re.DOTALL).groups()[0]
    # Coupling part
    coupling_re = "Real part\s*Imag part\s*Absolute(.+?)--"
    coupling_lists = get_block_lines(so_section, coupling_re)
    couplings = [(int(c[0]), int(c[3]), float(c[8]))
                  for c in coupling_lists]
    # Weights of original states in so-states
    weight_re = "Spin-free states, spin, and weights[\s\-]+(.+?)--"
    weight_lists = get_block_lines(so_section, weight_re)
    # Convert 'weight_lists' to an array and pull out the energies
    # to normalize them, that is setting the ground state energy to
    # 0.0 au
    weight_array = np.array(weight_lists, dtype=np.float)
    so_energies = weight_array[:, 1]
    weight_array[:, 1] = normalize_energies(so_energies)
    logging.warning("Energies get normalized too early. Fixme plz :)")
    # Transition moments
    trans_re = "Total A \(sec-1\)[\s\-]+(.+?)--"
    trans_lists = get_block_lines(so_section, trans_re)
    from_ids, to_ids, oscs, *_ = zip(*trans_lists)
    oscs = [float(osc) for osc in oscs]
    # Check out if MOLCAS left out any transitions.
    # Only consider excitations out of the GS
    to_ids = [int(to_id) for from_id, to_id
              in zip(from_ids, to_ids)
              if from_id == "1"]
    # If we got 40 SO states we expect ids from 2..40
    expected_ids = range(2, len(weight_lists)+1)
    missing_ids = set(expected_ids) - set(to_ids)
    for mid in missing_ids:
        logging.warning("No transition for SO states 1->{} found!".format(mid))
        # Insert f=-1 to symbolize the missing transition
        oscs.insert(mid-2, -1)
    # We expect len(expected_ids) oscillator strengths, because we
    # got no transition between the ground state and itself.
    oscs = oscs[:len(expected_ids)]

    so_states = list()
    for row, osc in zip(weight_array, oscs):
        state, E_global, *weight_line = row
        so_state = SpinOrbitState(state, E_global, weight_line, osc)
        so_states.append(so_state)

    return so_states, couplings


def split_states(states):
    all_ens = np.array([state.dE_global_eV for state in states])

    sing_inds = [i for i, state in enumerate(states)
                 if state.spin == 0.0]
    spins = [state.spin for state in states]
    trip_inds = [i for i, state in enumerate(states)
                 if state.spin == 1.0]
    sing_ens = all_ens[sing_inds]
    trip_ens = all_ens[trip_inds]
    assert((len(sing_ens) + len(trip_ens) == len(states)))

    return sing_ens, sing_inds, trip_ens, trip_inds, all_ens


def fuse_labels(ax, x, ys, labels):
    # y-values are in eV
    ys_diffs = np.array([ys[i+1]-ys[i] for i in range(len(ys)-1)])
    rel_diffs = ys_diffs / ys[1:]
    new_ys = [ys[0], ]
    new_labels = [labels[0], ]
    for i, yd in enumerate(rel_diffs, 1):
        if abs(yd) < 0.01:
            new_labels[-1] += " ," + labels[i]
        else:
            new_ys.append(ys[i])
            new_labels.append(labels[i])

    for y, l in zip(new_ys, new_labels):
        ax.text(x, y, l, va="center")


def to_texmath(sstr):
    sstr = sstr.replace("{}", "{{{}}}")
    sstr = "$" + sstr + "$"
    return sstr



def plot_states(sf_states, so_states, couplings=None, usetex=False):
    label_bases = ("S_{}", "SO_{}", "SO_{}", "T_{}")
    cpl_base = "{}-{} ({:.0f} cm⁻¹)"

    if usetex:
        rc("text", usetex=True)
        label_bases = [to_texmath(lb) for lb in label_bases]
        cpl_base = "${{{}}}-{{{}}} ~ ({:.0f} ~ \mathrm{{cm}}^{{-1}})$"
    # Spin-free states
    (sf_sing_ens, sf_sing_inds,
     sf_trip_ens, sf_trip_inds, sf_ens) = split_states(sf_states)
    # Spin-orbit states
    (so_sing_ens, so_sing_inds,
     so_trip_ens, so_trip_inds, so_ens) = split_states(so_states)

    fig, ax = plt.subplots()
    ax.set_ylabel("deltaE / eV")
    ax.set_xlabel("States")
    xlabels = ["SF singlet", "SO singlet", "SO triplet", "SF triplet"]
    ax.set_xticks(range(4))
    ax.set_xticklabels(xlabels)
    kwargs = dict(color="k", ls=" ", marker="_", ms=40)

    inds = (sf_sing_inds, so_sing_inds, so_trip_inds, sf_trip_inds)
    for i, ens in enumerate((sf_sing_ens, so_sing_ens,
                             so_trip_ens, sf_trip_ens)):
        xs = np.full_like(ens, i)
        ax.plot(xs, ens, **kwargs)
        # Add label
        labels = [label_bases[i].format(ind+1) for ind in inds[i]]
        fuse_labels(ax, i+.1, ens, labels)
        """
        for lx, ly, ind in zip(label_xs, ens, inds[i]):
            label_base = label_bases[i]
            ax.text(lx, ly, label_base.format(ind))
        """

    # Horizontal lines at 405 and 365 nm
    ax.axhline(y=3.4, color="k", linestyle="--")
    ax.axhline(y=3.06, color="k", linestyle="--")

    # Add couplings
    for from_id, to_id, abs_cpl in couplings:
        from_x = 1 if from_id-1 in so_sing_inds else 2
        to_x = 1 if to_id-1 in so_sing_inds else 2
        try:
            from_y = so_ens[from_id-1]
            to_y = so_ens[to_id-1]
            ax.add_line(Line2D((from_x, to_x),
                               (from_y, to_y)))
            x_text = -0.1 + from_x + (to_x - from_x) / 2
            y_text = from_y + (to_y - from_y) / 2
            cpl_str = cpl_base.format(from_id, to_id, abs_cpl)
            ax.text(x_text, y_text, cpl_str)
        except IndexError:
            logging.warning("IndexError")

    # Drop the first singlets energies
    all_ens = np.concatenate((sf_sing_ens[1:],
                              so_sing_ens[1:],
                              sf_trip_ens,
                              so_trip_ens))
    # The first two singlet states are not shown
    #ax.set_ylim((all_ens.min()*.9, all_ens.max()*1.1))
    ax.set_xlim(-0.5, 3.5)
    plt.tight_layout()
    plt.show()


def get_ground_state_conf(sf_states):
    sig_confs = significant_confs(sf_states[0].confs)
    sig_confs = sorted(sig_confs, key=lambda cf: -cf[-1])
    return sig_confs[0][0]


def make_html(so_states):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    j2_env = Environment(loader=FileSystemLoader(this_dir,
                                                 followlinks=True))
    tpl = j2_env.get_template("templates/sohtml.tpl")

    rendered = tpl.render(so_states=so_states)
    out_fn = os.path.join(
                os.getcwd(), "sorassi" + ".html")
    with open(out_fn, "w") as handle:
        handle.write(rendered)


def make_img_dict(sf_states, imgs, gs_conf):
    img_dict = dict()
    mo_fn_base = "mo_{}.irrep{}.png"
    for sf_state in sf_states:
        for conf, coef, weight in significant_confs(sf_state.confs):
            mo_pair = conf_diff(gs_conf, conf)
            if not mo_pair:
                continue
            from_mo, to_mo = mo_pair
            from_img = imgs[from_mo]
            to_img = imgs[to_mo]
            # Construct filenames
            from_fn = mo_fn_base.format(from_img, sf_state.jobiph)
            to_fn = mo_fn_base.format(to_img, sf_state.jobiph)
            img_dict.setdefault(sf_state.id, list()).append(
                (from_fn, to_fn, weight)
            )
    return img_dict


def set_sf_states(so_states, sf_states):
    for sos in so_states:
        for sfs_state in sos.sfs_states:
            sos.sf_states.append(sf_states[sfs_state-1])
            assert(sos.sf_states[-1].state == sfs_state)


def parse_args(args):
    parser = argparse.ArgumentParser("Parse SO-RASSI-calculations.")
    parser.add_argument("fn", help="SO-RASSI output to parse.")
    parser.add_argument("--cthresh", type=int, default=200,
                        help="Coupling threshold.")
    parser.add_argument("--cbelow", type=float, help="Only draw couplings "
                        "below this energy (in eV), e.g. around an "
                        "excitation energy. Helps to avoid cluttered plots "
                        "when many couplings are present.")
    parser.add_argument("--tex", action="store_true", default=False,
                        help="Use tex to render text.")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    fn_base = os.path.splitext(args.fn)[0]
    mo_names_dict = load_mo_names(fn_base)

    with open(args.fn) as handle:
        text = handle.read()
    # Spin-free section
    sf_states, trans_dict = parse_rassi(text)

    # Load images and mo_names if available and parse
    # the spin free states
    handle_spin_free_states(sf_states, trans_dict, mo_names_dict)

    so_states, couplings = parse_sorassi(text)
    set_sf_states(so_states, sf_states)
    # Only keep couplings above or equal to a threshold
    couplings = [c for c in couplings if c[2] >= args.cthresh]
    if args.cbelow:
        so_states_below = [sos.sostate for sos in so_states
                           if sos.dE_global_eV <= args.cbelow]
        logging.info("Drawing couplings between states below {} eV. "
                     "These are the SO states {}.".format(args.cbelow,
                     so_states_below)
        )

        couplings = [c for c in couplings
                     if (c[0] in so_states_below) and
                     (c[1] in so_states_below)]

    plot_states(sf_states, so_states, couplings, usetex=args.tex)
    make_html(so_states)


if __name__ == "__main__":
    run()
