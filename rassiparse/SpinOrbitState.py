#!/usr/bin/env python3

import numpy as np

class SpinOrbitState:

    def __init__(self, sostate, E_global, weight_line, osc):
        self.sostate = int(sostate)
        self.E_global = E_global
        self.osc = osc

        self.dE_global = self.E_global
        self.dE_global_eV = self.dE_global * 27.211386
        self.dE_global_nm = 45.5640 / self.dE_global
        self.sf_states_arr = np.array(weight_line).reshape((5, 3))
        self.spin = np.rint(np.sum(self.sf_states_arr[:, 1] *
                                   self.sf_states_arr[:, 2]))
        self.sfs_states = list()
        self.sf_states = list()
        self.weights = list()

        # Try to determine corresponding spin free state
        weight_sum = 0
        weight_sum_thresh = .8
        for sfs_state, sfs_spin, weight in self.sf_states_arr:
            weight_sum += weight
            self.sfs_states.append(int(sfs_state))
            self.weights.append(weight)

            if weight_sum >= weight_sum_thresh:
                break

    def set_sf_states(self, sf_states):
        for sfs_state, weight in zip(self.sfs_states, self.weights):
            self.sf_states.append((sf_states[sfs_state-1], weight))
            assert(self.sf_states[-1][0].state == sfs_state)
