#!/usr/bin/python3

import xraydb
import numpy as np
import sys

energy = float(sys.argv[1])

el = []
for i in range(len(sys.argv)-2):
    el.append(sys.argv[i+2])

q_vals = np.loadtxt("q_vals.dat")

# q_vals are momentum-transfer magnitudes |Q| = 4*pi*sin(theta)/lambda in nm^-1
# (the Fortran works in nm). xraydb.f0 expects sin(theta)/lambda in Ang^-1, so we
# divide by 4*pi (Q -> sin(theta)/lambda) and by 10 (nm^-1 -> Ang^-1).
with open("f_vals.dat", "w") as f:
    for q in q_vals:
        print(*[xraydb.f0(e, q/(4*np.pi)/10)[0] for e in el], file=f)

    for e in el:
        print(xraydb.f1_chantler(e, energy), xraydb.f2_chantler(e, energy), file=f)


