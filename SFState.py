#!/usr/bin/env python
# -*- coding: utf-8 -*-


class SFState:

    def __init__(self, id, symmetry, multiplicity, energy, confs):
        self.id = id
        self.symmetry = symmetry
        self.energy = energy
        self.multiplicity = multiplicity
        self.confs = confs

        self._spin = (self.multiplicity - 1) / 2

    @property
    def spin(self):
        return self._spin
