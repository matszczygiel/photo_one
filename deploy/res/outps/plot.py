#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

exp_file = "/home/mateusz/Documents/photo/h/analitical/dip.dat"

infiles = sys.argv[1:]
exp_data = np.loadtxt(exp_file)

include = -1

plt.xlim((exp_data[0, 0], exp_data[-1, 0]))
plt.ylim((0.5, 7.))
xs = np.arange(exp_data[0, 0], exp_data[-1, 0], 0.01)
cs = scipy.interpolate.CubicSpline(exp_data[:, 0], exp_data[:, 2])

plt.plot(exp_data[0, 0], exp_data[0, 2], 'o', markersize=3, label="numerical")
plt.plot(xs, cs(xs))


for infile in infiles:
    file = open(infile, 'r')
    inp = file.readlines()
    file.close()

    stride = 8
    phind = range(1, len(inp), stride)

    phot = []
    for i in phind:
        words = inp[i].split()
        phot.append(float(words[-1]))

    csind = range(6, len(inp), stride)
    crss = []
    for i in csind:
        words = inp[i].split()
        crss.append(float(words[include]))

    xs = np.arange(phot[0], phot[-1], 0.01)
    cs = scipy.interpolate.CubicSpline(phot, crss)

    plt.plot(phot, crss, 'o', markersize=3, label=infile)
    plt.plot(xs, cs(xs))


plt.title("H", va='bottom')
plt.legend()
plt.savefig("res.png")