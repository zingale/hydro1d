#!/usr/bin/env python3

from __future__ import print_function

import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'


def plot_sim(filename, axes, label):

    data = np.loadtxt(filename)

    x = data[:,0]
    rho = data[:,1]
    u = data[:,4]
    p = data[:,5]
    e = data[:,6]

    ax = axes[0]
    ax.scatter(x, rho, marker="x", s=7, label=label, color="C1")

    ax = axes[1]
    ax.scatter(x, u, marker="x", s=7, label=label, color="C1")

    ax = axes[2]
    ax.scatter(x, p, marker="x", s=7, label=label, color="C1")

    ax = axes[3]
    ax.scatter(x, e, marker="x", s=7, label=label, color="C1")


def plot_exact(filename, axes):

    # get the exact solution
    exact = np.loadtxt(filename)

    x_exact   = exact[:,0]
    rho_exact = exact[:,1]
    u_exact   = exact[:,2]
    p_exact   = exact[:,3]
    e_exact   = exact[:,4]


    ax = axes[0]

    ax.plot(x_exact, rho_exact, label="exact")

    ax.set_ylabel(r"$\rho$")
    ax.set_xlim(0,1.0)
    ax.set_ylim(0,1.1)

    ax = axes[1]
    
    ax.plot(x_exact, u_exact, label="exact")

    ax.set_ylabel(r"$u$")
    ax.set_xlim(0,1.0)

    ax = axes[2]

    ax.plot(x_exact, p_exact, label="exact")

    ax.set_ylabel(r"$p$")
    ax.set_xlim(0,1.0)
    
    ax = axes[3]

    ax.plot(x_exact, e_exact, label="exact")

    ax.set_xlabel(r"x")
    ax.set_ylabel(r"$e$")
    ax.set_xlim(0,1.0)


if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("--sim", help="simulation output file", default=None, type=str)
    p.add_argument("--exact", help="exact solution", default=None, type=str)
    p.add_argument("--label", help="label for the simulation", default="simulation", type=str)
    p.add_argument("-o", help="output png file name", default="plot.png", type=str)
    args = p.parse_args()
    
    # plot
    fig, axes = plt.subplots(nrows=4, ncols=1, num=1)

    if args.sim is not None:
        plot_sim(args.sim, axes, args.label)
    
    if args.exact is not None:
        plot_exact(args.exact, axes)


    axes[0].legend(frameon=False, fontsize="medium")

    fig.set_size_inches(7,10.0)
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)                                  

    plt.tight_layout()

    plt.savefig(args.o, bbox_inches="tight", dpi=100)
