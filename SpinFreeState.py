#!/usr/bin/env python
# -*- coding: utf-8 -*-

from WithHeader import WithHeader

class SpinFreeState(WithHeader):
    headers = {
        "state": ("State", "{}"),
        "jobiph": ("JobIph", "{}"),
        "root": ("Root", "{}"),
        "sym": ("Sym.", "{}"),
        "mult" : ("2S+1", "{}"),
        "energy": ("E / au", "{:.6f}"),
        "confdiffs": ("Transitions", "{}"),
        "dE_global": ("ΔE / au", "{:.2f}"),
        "dE_global_eV": ("ΔE / eV", "{:.2f}"),

        "state_rel" : ("#", "{}"),
        "dE_gs": ("ΔE / au", "{:.2f}"),
        "dE_gs_eV": ("ΔE / eV", "{:.2f}"),
        "dE_gs_nm": ("λ / nm", "{:.1f}"),
        "osc": ("f", "{:.4f}"),
    }

    def __init__(self, state, jobiph, root, sym, mult,
                 confs, energy, dE_global):
        super(SpinFreeState, self).__init__()

        self.state = state
        self.jobiph = jobiph
        self.root = root
        self.sym = sym
        self.mult = mult
        self.confs = confs
        self.energy = energy
        self.dE_global = dE_global

        self._confdiffs = list()

    @property
    def confdiffs(self):
        return ", ".join([str(t) for t in self._confdiffs])

    def add_confdiff(self, confdiff):
        self._confdiffs.append(confdiff)

    @property
    def dE_global_eV(self):
        return self.dE_global * self.hartree2eV

    def set_ground_state(self, ground_state, osc):
        self.ground_state = ground_state
        self.dE_gs = self.energy - self.ground_state.energy
        self.dE_gs_eV = self.dE_gs * self.hartree2eV
        self.dE_gs_nm = self.hartree2nm / self.dE_gs
        self.osc = osc

    def set_images(self, image_list):
        mo_nums, mo_fns, = zip(*[mo_tpl for mo_tpl in image_list])
        self.mo_nums = mo_nums
        self.mo_fns = mo_fns
        for cd in self._confdiffs:
            cd.set_mo_nums(self.mo_nums)
