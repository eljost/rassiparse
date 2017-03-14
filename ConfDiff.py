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
            mo_pairs = [("?", "?")]
        self.mo_pairs = mo_pairs
        self.weight = weight

        self.mo_pair_str = ", ".join(["{} -> {}".format(from_, to)
                                     for from_, to in self.mo_pairs])

    def __str__(self):
        return "{} ({:.1%})".format(self.mo_pair_str, self.weight)
