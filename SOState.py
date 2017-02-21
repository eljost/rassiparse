#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


class SOState:

    def __init__(self, id, energy, weight_line):
        self.id = id
        self.energy = energy
        self.weights = np.array(weight_line).reshape((5, 3))

        self._spin = np.rint(np.sum(self.weights[:, 1] * self.weights[:, 2]))

    @property
    def spin(self):
        return self._spin

    def add_coupling(self, coupling_line):
        pass
