from __future__ import print_function
import os,sys
import argparse

import pysam
from collections import defaultdict

import math

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd


import networkx as nx
import argparse, os
import pysam

from collections import defaultdict


def plot_V_all_tools(args):
    # sns.plt.clf()
    with sns.plotting_context("paper", font_scale=2.4):
        rc={'axes.labelsize': 16.0, 'font.size': 16.0, 'legend.fontsize': 10.0, 'axes.titlesize': 14, 'xtick.labelsize': 16.0, 'ytick.labelsize' : 16.0}
        sns.set(rc=rc)

        # sns.set()
        indata = pd.read_csv(args.infile)
        indata = indata.dropna()

        g = sns.pointplot(x="dataset", y="V", hue="tool", data=indata)
        g.set(ylim=(0.5, 1))
        plt.xticks(rotation=90)
        # plt.legend(loc='lower right')
        plt.tight_layout()
        # g.set_title('Clustering accuracy per class size')
        # g.set_ylabel("Performance")
        # g.set_xlabel("Class size")
        plt.savefig(os.path.join(args.outfolder, "v.pdf" ))
        plt.clf()

        g = sns.pointplot(x="dataset", y="C", hue="tool", data=indata)
        g.set(ylim=(0.5, 1))
        plt.xticks(rotation=90)
        # plt.legend(loc='lower right')
        plt.tight_layout()
        # g.set_title('Clustering accuracy per class size')
        # g.set_ylabel("Performance")
        # g.set_xlabel("Class size")
        plt.savefig(os.path.join(args.outfolder, "C.pdf" ))
        plt.clf()


        g = sns.pointplot(x="dataset", y="H", hue="tool", data=indata)
        g.set(ylim=(0.5, 1))
        plt.xticks(rotation=90)
        # plt.legend(loc='lower right')
        plt.tight_layout()
        # g.set_title('Clustering accuracy per class size')
        # g.set_ylabel("Performance")
        # g.set_xlabel("Class size")
        plt.savefig(os.path.join(args.outfolder, "H.pdf" ))
        plt.clf()


        g = sns.pointplot(x="dataset", y="percent_nontrivially_clustered", hue="tool", data=indata)
        g.set(ylim=(0, 100))
        plt.xticks(rotation=90)
        # plt.legend(loc='lower right')
        plt.tight_layout()
        # g.set_title('Clustering accuracy per class size')
        # g.set_ylabel("Performance")
        # g.set_xlabel("Class size")
        plt.savefig(os.path.join(args.outfolder, "percent_nontrivially_clustered.pdf" ))
        plt.clf()

        plt.close()


def plot_V_all_tools2(args):
    with sns.plotting_context("paper", font_scale=2.4):
        # sns.set(style="whitegrid")
        rc={'axes.labelsize': 20.0, 'font.size': 20.0, 'legend.fontsize': 20.0, 'axes.titlesize': 20, 'xtick.labelsize': 24.0, 'ytick.labelsize' : 20.0}
        sns.set(rc=rc, style="whitegrid")

        # Load the dataset
        # crashes = sns.load_dataset("car_crashes")
        # print(crashes)
        indata = pd.read_csv(args.infile)
        # indata = indata.dropna()   

        # Make the PairGrid
        
        g = sns.PairGrid(indata.sort_values("dataset", na_position="first"), hue="tool", hue_order=["linclust", "isoseq3", "carnac", "qubric"],
                         x_vars=indata.columns[2:], y_vars=["dataset"], 
                         height=10, aspect=.5)

        # Draw a dot plot using the stripplot function
        g.map(sns.stripplot, size=20, orient="h", jitter=0.0, #palette="ch:s=1,r=-.1,h=1_r",
               linewidth=1, edgecolor="w", alpha=.8)

        # Use the same x axis limits on all columns and add better labels
        # g.set(xlim=[(0, 1), (0, 1),(0, 1), (0, 100) ], xlabel=["", "", "", "", ""], ylabel=["", "", "", "", ""])
        g.set(xlabel="", ylabel="")
        g.axes[0,0].set_xlim((0.5,1))
        g.axes[0,1].set_xlim((0.5,1))
        g.axes[0,2].set_xlim((0.5,1))
        g.axes[0,3].set_xlim((0,100))
        g.set(xticks=[0.5,0.75,1.0])
        g.axes[0,3].set_xticks((0,50,100))

        # Use semantically meaningful titles for the columns
        titles = ["V-measure", "Completeness", "Homogeneity",
                  "%-nontrivially clustered"]

        for ax, title in zip(g.axes.flat, titles):

            # Set a different title for each axes
            ax.set(title=title)

            # Make the grid horizontal instead of vertical
            ax.xaxis.grid(False)
            ax.yaxis.grid(True)

        sns.despine(left=True, bottom=True)
        plt.tight_layout()

        g = g.add_legend()
        # handles = ["linclust", "isoseq3", "carnac", "qubric"]
        # plt.legend()
        # plt.yticks(rotation=90)

        plt.savefig(os.path.join(args.outfolder, "test.pdf" ))
        plt.clf()

        plt.close()



def main(args):
    plot_V_all_tools2(args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--infile', type=str, help='csv infile')
    parser.add_argument('--outfolder', type=str, help='Output file with results')
    args = parser.parse_args()

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)



