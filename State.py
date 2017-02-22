#!/usr/bin/env python
# -*- coding: utf-8 -*-


class State(object):

    def __init__(self, id, Eau):
        self.id = id
        self.Eau = Eau
        self._Enm = 45.56335 / self.Eau
        self._EeV = self.Eau * 27.211386

    @property
    def Enm(self):
        return self._Enm

    @property
    def EeV(self):
        return self._EeV
