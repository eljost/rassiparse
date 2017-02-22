#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging

import numpy as np

from State import State


class SOState(State):

    def __init__(self, id, Eau, weight_line, osc):
        super(SOState, self).__init__(id, Eau)
        self.weights = np.array(weight_line).reshape((5, 3))
        self.osc = osc

        self._spin = np.rint(np.sum(self.weights[:, 1] * self.weights[:, 2]))
        self._sf_state = None

    @property
    def spin(self):
        return self._spin

    @property
    def sf_state(self):
        sfs, _, weight = self.weights[0]
        if (weight >= .9):
            return sfs
        else:
            logging.warning("No clear correspondence to one spin free state")
            return None
