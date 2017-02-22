#!/usr/bin/env python
# -*- coding: utf-8 -*-


from State import State

class SFState(State):

    def __init__(self, id, symmetry, multiplicity, Eau, confs):
        super(SFState, self).__init__(id, Eau)
        self.symmetry = symmetry
        self.multiplicity = multiplicity
        self.confs = confs

        self._spin = (self.multiplicity - 1) / 2

    @property
    def spin(self):
        return self._spin
