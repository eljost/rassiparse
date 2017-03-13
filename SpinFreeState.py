#!/usr/bin/env python
# -*- coding: utf-8 -*-

from WithHeader import WithHeader

class SpinFreeState(WithHeader):
    headers = {
        "state": "State",
        "jobiph": "JobIph",
        "root": "Root",
        "sym": "Sym.",
        "mult" : "2S+1",
        "energy": "E / au",
        "energy_rel": "ΔE / au",
        "transitions": "Transitions",
    }

    def __init__(self, state, jobiph, root, sym, mult,
                 confs, energy, energy_rel):
        super(SpinFreeState, self).__init__()

        self.state = state
        self.jobiph = jobiph
        self.root = root
        self.sym = sym
        self.mult = mult
        self.confs = confs
        self.energy = energy
        self.energy_rel = energy_rel

        self._transitions = list()

    @property
    def transitions(self):
        return ", ".join([str(t) for t in self._transitions])

    def add_transition(self, transition):
        self._transitions.append(transition)
