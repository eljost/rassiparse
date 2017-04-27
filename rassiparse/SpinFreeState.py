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
        "confdiffsw": ("Transitions", "{}"),
        "confdiff_strs": ("Transitions", "{}"),
        "confdiffocc_strs": ("Configuration", "{}"),
        "confdiffoccw_strs": ("Configuration", "{}"),
        "weights" : ("%", "{}"),
        "dE_global": ("ΔE / au", "{:.2f}"),
        "dE_global_eV": ("ΔE / eV", "{:.2f}"),

        "state_rel": ("#", "{}"),
        "dE_gs": ("ΔE / au", "{:.2f}"),
        "dE_gs_eV": ("ΔE / eV", "{:.2f}"),
        "dE_gs_nm": ("λ / nm", "{:.1f}"),
        "osc": ("f", "{:.4f}"),
    }

    _mult_label = {
        1: "S",
        2: "D",
        3: "T",
        4: "Q",
        5: "Quintet", 
        6: "Sextet",
    }


    def __init__(self, state, jobiph, root, sym, mult,
                 confs, energy, dE_global):
        super(SpinFreeState, self).__init__()

        self.state = state
        self.jobiph = jobiph
        self.root = root
        self.sym = sym
        self.mult = mult
        self.mult_label = self._mult_label[mult]
        self.confs = confs
        self.energy = energy
        self.dE_global = dE_global

        self._confdiffs = list()
        self.spin = (self.mult - 1) / 2

    def _cds_sorted(self):
        return sorted([cd for cd in self._confdiffs],
                      key=lambda cd: -cd.weight)

    @property
    def confdiffsw(self):
        return ", ".join([cd.str_with_weight() for cd in self._cds_sorted()])

    @property
    def confdiff_strs(self):
        return ", ".join([str(cd) for cd in self._cds_sorted()])

    @property
    def confdiffoccw_strs(self):
        return ", ".join(["{} ({:.1%})".format(cd.conf_occ, cd.weight)
                          for cd in self.confdiffs])

    @property
    def confdiffocc_strs(self):
        return ", ".join(["{}".format(cd.conf_occ)
                          for cd in self.confdiffs])

    @property
    def confdiffs(self):
        return self._confdiffs

    @property
    def confdiff_images(self):
        return [(cd.weight, cd.mo_images, cd)
                for cd in self._cds_sorted()
                if (cd.mo_pairs is not None) and (cd.mo_images is not None)]

    def add_confdiff(self, confdiff):
        self._confdiffs.append(confdiff)

    @property
    def dE_global_eV(self):
        return self.dE_global * self.hartree2eV

    @property
    def weights(self):
        return ", ".join(["{:.0f}".format(cd.weight*100)
                          for cd in self._confdiffs])

    def set_ground_state(self, ground_state, osc):
        self.ground_state = ground_state
        self.dE_gs = self.energy - self.ground_state.energy
        self.dE_gs_eV = self.dE_gs * self.hartree2eV
        self.dE_gs_nm = self.hartree2nm / self.dE_gs
        self.osc = osc

    def set_mo_images(self, image_list):
        mo_nums, mo_images, = zip(*[mo_tpl for mo_tpl in image_list])
        self.mo_nums = mo_nums
        self.mo_images = mo_images
        for cd in self._confdiffs:
            cd.set_mo_nums_images(self.mo_nums, self.mo_images)

    def set_mo_names(self, names_list):
        self.mo_names = names_list
        for cd in self._confdiffs:
            cd.set_mo_names(self.mo_names)

    def as_str_list(self, attrs, newlines=False):
        str_list = super(SpinFreeState, self).as_str_list(attrs)
        if not newlines:
            return str_list
        if "weights" in attrs:
            weights_str = "\n".join(["{:.0f}".format(cd.weight*100)
                                     for cd in self._confdiffs])
            weight_index = attrs.index("weights")
            str_list[weight_index] = weights_str
        if "confdiffs" in attrs:
            cdstr = "\n".join([str(cd) for cd in self._confdiffs])
            cd_idx = attrs.index("confdiffs")
            str_list[cd_idx] = cdstr

        return str_list

    def __str__(self):
        return "SpinFreeState {} from JOB{:0>3}".format(self.state,
                                                        self.jobiph)
