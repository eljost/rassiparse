#!/usr/bin/env python3

import numpy as np

class SpinOrbitState:

    def __init__(self, sostate, E_global, weight_line, osc):
        self.sostate = int(sostate)
        self.E_global = E_global
        self.osc = osc

        self.dE_global_eV = self.E_global * 27.211386
        self.weights = np.array(weight_line).reshape((5, 3))
        self.spin = np.rint(np.sum(self.weights[:, 1] * self.weights[:, 2]))
