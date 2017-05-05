#!/usr/bin/env python3

from __future__ import print_function

import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams["text.usetex"] = True
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

# font sizes
mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'


def plot_sim(filename, ax, label, lines=False):

    data = np.loadtxt(filename)

    x = data[:,0]
    rho = data[:,1]
    u = data[:,4]
    p = data[:,5]
    e = data[:,6]

    if lines:
        ax.plot(x, rho, label=label, color="C0", zorder=-100)
    else:
        ax.scatter(x, rho, marker="x", s=7, label=label, color="C1")


if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("--sim", help="simulation output file", default=None, type=str)
    p.add_argument("--exact", help="exact solution", default=None, type=str)
    p.add_argument("--label", help="label for the simulation", default="simulation", type=str)
    p.add_argument("-o", help="output png file name", default="plot.png", type=str)
    args = p.parse_args()

    # plot
    fig, axes = plt.subplots(nrows=1, ncols=1, num=1)

    if args.sim is not None:
        plot_sim(args.sim, axes, args.label)

    if args.exact is not None:
        plot_sim(args.exact, axes, "exact", lines=True)


    axes.legend(frameon=False, fontsize="medium")

    fig.set_size_inches(6, 5)
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

    plt.tight_layout()

    plt.savefig(args.o, bbox_inches="tight", dpi=100)
