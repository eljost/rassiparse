#!/usr/bin/env python3

from collections import namedtuple

def conv_coupling_line(coupling_line):
    i1, s1, ms1, i2, s2, ms2, real, imag, abs_ = coupling_line
    return (int(i1), float(s1), float(ms1),
            int(i2), float(s2), float(ms2),
            float(real), float(imag), float(abs_)
    )

Coupling = namedtuple("Coupling",
                      "i1 s1 ms1 i2 s2 ms2 real imag abs"
)
