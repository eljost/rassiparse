#!/usr/bin/env python
# -*- coding: utf-8 -*-


class State(object):

    _mult_label = {
        0.0: "S",
        0.5: "D",
        1.0: "T",
        1.5: "Q",
    }

    def __init__(self, id, Eau):
        self.id = int(id)
        self.Eau = Eau
        self._Enm = 45.56335 / self.Eau
        self._EeV = self.Eau * 27.211386

    @property
    def Enm(self):
        return self._Enm

    @property
    def EeV(self):
        return self._EeV

    @property
    def mult_label(self):
        return self._mult_label[self.spin]
