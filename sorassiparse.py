#!/usr/bin/env python
# -*- coding: utf-8 -*-


import re

from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np


from SFState import SFState
from SOState import SOState


def get_block_lines(text, regex):
    mobj = re.search(regex, text, re.DOTALL)
    block = mobj.groups()[0].strip()
    as_lists = [line.strip().split() for line in block.split("\n")]
    return as_lists


def normalize_energies(energies):
    energies = np.array(energies, dtype=np.float)
    energies -= energies.min()
    energies *= 27.2114
    return energies


def parse_sorassi(text):
    # Spin-free section
    sf_state_re = "state\s*(\d+).+?symmetry\s*=\s*(\d+).+?Spin multiplic=\s*(\d+)"
    sf_states_lists = re.findall(sf_state_re, text, re.DOTALL)
    # Transform sf_states into a dict
    # State: (symmetry, multiplicity)
    sf_states_dict = {int(id): (int(sym), int(mult))
                      for id, sym, mult in sf_states_lists}

    sf_energies_re = "SF State.+?\n(.+?)\+\+"
    sf_energies_lists = get_block_lines(text, sf_energies_re)
    sf_states, sf_energies, *_ = zip(*sf_energies_lists)
    sf_states = np.array(sf_states, dtype=np.int)
    #sf_energies = normalize_energies(sf_energies)
    sf_states = list()
    for key in sf_states_dict:
        sym, mult = sf_states_dict[key]
        sf_state = SFState(key, sym, mult, sf_energies[key-1])
        sf_states.append(sf_state)

    # Spin-orbit section
    so_section_re = "Spin-orbit section(.+)"
    so_section = re.search(so_section_re, text, re.DOTALL).groups()[0]
    # Coupling part
    coupling_re = "Real part\s*Imag part\s*Absolute(.+?)--"
    coupling_lists = get_block_lines(so_section, coupling_re)
    couplings = [(int(c[0]), int(c[3]), float(c[8]))
                  for c in coupling_lists]
    print(couplings)
    # Weights of original states in so-states
    weight_re = "Spin-free states, spin, and weights[\s\-]+(.+?)--"
    weight_lists = get_block_lines(so_section, weight_re)
    weight_array = np.array(weight_lists, dtype=np.float)
    # Transition moments
    trans_re = "Total A \(sec-1\)[\s\-]+(.+?)--"
    trans_lists = get_block_lines(so_section, trans_re)

    so_states = list()
    for row in weight_array:
        id, so_energy, *weight_line = row
        so_state = SOState(id, so_energy, weight_line)
        so_states.append(so_state)

    return sf_states, so_states, couplings


def split_states(states):
    all_ens = np.array([state.energy for state in states])
    all_ens = normalize_energies(all_ens)

    sing_inds = [i for i, state in enumerate(states)
                 if state.spin == 0.0]
    trip_inds = [i for i, state in enumerate(states)
                 if state.spin == 1.0]
    sing_ens = all_ens[sing_inds]
    trip_ens = all_ens[trip_inds]
    assert((len(sing_ens) + len(trip_ens) == len(states)))

    return sing_ens, sing_inds, trip_ens, trip_inds, all_ens


def plot_states(sf_states, so_states, couplings=None):
    # Spin-free states
    (sf_sing_ens, sf_sing_inds,
     sf_trip_ens, sf_trip_inds, sf_ens) = split_states(sf_states)
    # Spin-orbit states
    (so_sing_ens, so_sing_inds,
     so_trip_ens, so_trip_inds, so_ens) = split_states(so_states)

    fig, ax = plt.subplots()
    kwargs = dict(color="k", ls=" ", marker="_", ms=40)
    for i, ens in enumerate((sf_sing_ens, sf_trip_ens,
                            so_sing_ens, so_trip_ens)):
        xs = np.full_like(ens, i)
        ax.plot(xs, ens, **kwargs)
    ax.axhline(y=3.4, color="k", linestyle="--")
    for from_id, to_id, abs_cpl in couplings:
        from_x = 2 if from_id-1 in so_sing_inds else 3
        to_x = 2 if to_id-1 in so_sing_inds else 3
        from_y = so_ens[from_id-1]
        to_y = so_ens[to_id-1]
        """
        ax.annotate(str(coupling_abs),
                    xy=(to_x, to_y), xycoords="data",
                    #textcoords="offset points",
                    xytext=(from_x, from_y),
                    arrowprops=dict(arrowstyle="<->"))
        """
        ax.add_line(Line2D((from_x, to_x),
                           (from_y, to_y)))
        x_text = -0.1 + from_x + (to_x - from_x) / 2
        y_text = from_y + (to_y - from_y) / 2
        cpl_str = "{}<->{} ({:.1f})".format(from_id, to_id, abs_cpl)
        ax.text(x_text, y_text, cpl_str)

    # Drop the first singlets energies
    all_ens = np.concatenate((sf_sing_ens[1:],
                              so_sing_ens[1:],
                              sf_trip_ens,
                              so_trip_ens))
    # The first two singlet states are not shown
    #ax.set_ylim((all_ens.min()*.9, all_ens.max()*1.1))
    ax.set_xlim(-0.5, 3.5)
    plt.show()

if __name__ == "__main__":
    fn = "/home/carpx/Code/sorassiparse/logs/12_sorassi.out"
    with open(fn) as handle:
        text = handle.read()
    sf_states, so_states, couplings = parse_sorassi(text)
    print(couplings)
    plot_states(sf_states, so_states, couplings)
