#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Parse caspt2 output and search for intruder states

import re
import sys

import numpy as np
from tabulate import tabulate

import rex

class CASPT2Root:

    def __init__(self, total_e, ref_weight, problems):
        self.total_e = total_e
        self.ref_weight = ref_weight
        self.problems = problems

    def filter_for(self, cond):
        return tuple([item for item in self.problems if cond(item)])

    def small_shifted_denoms(self, shift, thresh=0.1):
        # Check if added level shift accidentally creates near zero
        # denominators
        return self.filter_for(
            lambda item: abs(item[4] + shift) <= thresh)

    def high_energy_contribs(self, thresh=0.001):
        return self.filter_for(
            lambda item: abs(item[7]) >= thresh)

    def high_coefficients(self, thresh=0.05):
        return self.filter_for(
            lambda item: abs(item[6]) >= thresh)

    def unsafe(self):
        return self.filter_for(
            # Denomiator below 1 and the absolute value of the energy-
            # contribution greater than 1e-4
            lambda item: (item[4] < 1) and (abs(item[7]) >= 0.001))

    def check(self, shift):
        self.problem_set = set(self.small_shifted_denoms(shift))
        self.problem_set.update(self.high_energy_contribs())
        self.problem_set.update(self.unsafe())
        self.problem_set.update(self.high_coefficients())

def parse_file(fn, headers):
    with open(fn) as handle:
        text = handle.read()

    # Search all reference weights
    ref_weights = rex.find_floats(text, "Reference weight:")

    # Determine number of calculated roots from the number of found reference
    # weights.
    number_of_roots = len(ref_weights)

    # Limit number of found 'Total energies' to the number of the actual roots,
    # because this string appears multiple times in the output
    total_energies = rex.find_floats(text, "Total energy:", number_of_roots)

    # Prepare regular expressions to parse listed problems (small denom., etc.)
    start_match = "\s+".join(headers)
    end_match = "--"
    findall_re = start_match + "(.+?)" + end_match
    roots_raw = re.findall(findall_re, text, re.DOTALL)

    # Extract reports on small energy denominators, large coefficients etc.
    prob_re_list = ["[A-Z]+", int, "[\w\.]+", "[\w\. ]+",
                   float, float, float, float]
    prob_re, prob_conv = rex.join_re(prob_re_list)

    root_problems = [rex.match_lines(root.split("\n"), prob_re, prob_conv)
        for root in roots_raw]
    zipped = zip(total_energies, ref_weights, root_problems)
    roots = [CASPT2Root(tot_e, ref_w, root_prob)
        for tot_e, ref_w, root_prob in zipped]

    return roots

def run(fn, shift):
    print("Checking file '{0}'.\n".format(fn))

    headers = ("CASE", "SYMM", "ACTIVE-MIX", "NON-ACTIVE INDICES",
        "DENOMINATOR", "RHS VALUE", "COEFFICIENT", "CONTRIBUTION")

    roots = parse_file(fn, headers)

    for i, root in enumerate(roots):
        root.check(shift)
        if root.problem_set:
            print("#" * 10 + "Root {0}".format(i + 1) + "#"*10)
            print(tabulate(root.problem_set, headers=headers, floatfmt=".5f"))
            print()

    tot_energies, ref_weights = zip(*[(r.total_e, r.ref_weight) for r in roots])
    # Normalized to the energy of root 1
    tot_e_norm = np.array(tot_energies) - tot_energies[0]
    tot_e_ev = tot_e_norm * 27.21138
    zipped = zip(range(1, len(roots) + 1), tot_energies, tot_e_norm,
                tot_e_ev, ref_weights)
    print(tabulate(zipped,
        headers=("Root", "E in a.u.", "E-E(GS) in a.u", "E-E(GS) in eV",
                "Ref. weight"),
        floatfmt=".5f"))
    print()
    print("Ref. weights:\t{:.1%}+-{:.1%}".format(
            np.average([r.ref_weight for r in roots]),
            np.std([r.ref_weight for r in roots])))

if __name__ == "__main__":
    fn = sys.argv[1]
    shift = float(sys.argv[2])

    run(fn, shift)
