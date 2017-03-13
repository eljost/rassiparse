#!/usr/bin/env python
# -*- coding: utf-8 -*-

from WithHeader import WithHeader

class Transition(WithHeader):
    headers = {
        "delta_E": "ΔE / au",
        "delta_E_eV": "ΔE / eV",
        "delta_E_nm": "ΔE / nm",
        "osc": "f",
    }
    
    def __init__(self, from_state, to_state, mo_pair, osc, weight):
        super(Transition, self).__init__()

        self.from_state = from_state
        self.to_state = to_state
        self.mo_pair = mo_pair
        self.osc = osc
        self.weight = weight

        self.from_energy = from_state.energy
        self.to_energy = to_state.energy

        self.delta_E = self.to_energy - self.from_energy
        self.delta_E_eV = self.delta_E * 27.211386
        self.delta_E_nm = 45.5640 / self.delta_E 

    def __str__(self):
        return "{} -> {} ({:.1%})".format(*self.mo_pair, self.weight)
