#!/usr/bin/env python
# -*- coding: utf-8 -*-

from WithHeader import WithHeader

class ConfDiff(WithHeader):
    headers = {
        "delta_E": "ΔE / au",
        "delta_E_eV": "ΔE / eV",
        "delta_E_nm": "ΔE / nm",
        "osc": "f",
    }
    
    def __init__(self, mo_pairs, weight):
        super(ConfDiff, self).__init__()

        if mo_pairs == []:
            mo_pairs = None
        self.mo_pairs = mo_pairs
        self.weight = weight
        self.mo_nums = None

    def set_mo_nums(self, mo_num_list):
        if not self.mo_pairs:
            return
        self.mo_nums = [(mo_num_list[from_mo], mo_num_list[to_mo])
                        for from_mo, to_mo
                        in self.mo_pairs]

    def __str__(self):
        if not self.mo_pairs:
            return "({:.1%})".format(self.weight)
        if self.mo_nums:
            iterate_over = self.mo_nums
        else:
            iterate_over = self.mo_pairs

        mo_pair_str = ", ".join(["{} -> {}".format(from_, to)
                                     for from_, to in iterate_over])
        return "{} ({:.1%})".format(mo_pair_str, self.weight)
