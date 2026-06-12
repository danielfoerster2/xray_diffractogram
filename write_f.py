#!/usr/bin/python3

import xraydb
import numpy as np
import sys

energy = float(sys.argv[1])
outfile = sys.argv[2]

el = sys.argv[3:]

q_vals = np.loadtxt("q_vals.dat")

with open(outfile, "w") as f:
    for q in q_vals:
        print(*[xraydb.f0(e, q/(10*4*np.pi))[0] for e in el], file=f)

    for e in el:
        print(xraydb.f1_chantler(e, energy), xraydb.f2_chantler(e, energy), file=f)


