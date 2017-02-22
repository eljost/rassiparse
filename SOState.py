#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from State import State


class SOState(State):

    def __init__(self, id, Eau, weight_line, osc):
        super(SOState, self).__init__(id, Eau)
        self.weights = np.array(weight_line).reshape((5, 3))
        self.osc = osc

        self._spin = np.rint(np.sum(self.weights[:, 1] * self.weights[:, 2]))

    @property
    def spin(self):
        return self._spin
