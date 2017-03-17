#!/usr/bin/env python
# -*- coding: utf-8 -*-

class ConfDiff:
    
    def __init__(self, mo_pairs, weight):

        if mo_pairs == []:
            mo_pairs = None
        self.mo_pairs = mo_pairs
        self.weight = weight

        self.mo_nums = None

    def set_mo_nums_images(self, mo_num_list, mo_images):
        if not self.mo_pairs:
            return
        self.mo_nums = list()
        self.mo_images = list()
        for from_mo, to_mo in self.mo_pairs:
            self.mo_nums.append(
                        (mo_num_list[from_mo], mo_num_list[to_mo])
            )
            self.mo_images.append(
                        (mo_images[from_mo], mo_images[to_mo])
            )

    def str_with_weight(self):
        if not self.mo_pairs:
            return "({:.1%})".format(self.weight)
        mo_pair_str = self.__str__()
        return "{} ({:.1%})".format(mo_pair_str, self.weight)

    def __str__(self):
        if not self.mo_pairs:
            return ""
            #return None
        if self.mo_nums:
            iterate_over = self.mo_nums
        else:
            iterate_over = self.mo_pairs

        return ", ".join(["{} → {}".format(from_, to)
                         for from_, to in iterate_over])