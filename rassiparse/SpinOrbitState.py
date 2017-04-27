#!/usr/bin/env python3

import numpy as np

class SpinOrbitState:

    def __init__(self, sostate, E_global, weight_line, osc):
        self.sostate = int(sostate)
        self.E_global = E_global
        self.osc = osc

        self.dE_global_eV = self.E_global * 27.211386
        self.sf_states_arr = np.array(weight_line).reshape((5, 3))
        self.spin = np.rint(np.sum(self.sf_states_arr[:, 1] *
                                   self.sf_states_arr[:, 2]))
        self.sfs_states = list()
        self.sf_states = list()

        # Try to determine corresponding spin free state
        weight_sum = 0
        weight_sum_thresh = .8
        for sfs_state, sfs_spin, sfs_weight in self.sf_states_arr:
            weight_sum += sfs_weight
            self.sfs_states.append(int(sfs_state))

            if weight_sum >= weight_sum_thresh:
                break
