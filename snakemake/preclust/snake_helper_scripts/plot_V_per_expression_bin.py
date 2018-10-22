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


def plot_V_per_expression_bin(args):
    # sns.plt.clf()
    with sns.plotting_context("paper", font_scale=2.4):
        rc={'axes.labelsize': 16.0, 'font.size': 16.0, 'legend.fontsize': 10.0, 'axes.titlesize': 14, 'xtick.labelsize': 16.0, 'ytick.labelsize' : 16.0}
        sns.set(rc=rc)

        # sns.set()
        indata = pd.read_csv(args.infile)
        # g = sns.pointplot(x="class_size", y="measure", data=indata, hue="measure_type")

        # g = sns.factorplot(x="w", y="P_minimizer_shared", col="error_rate", row="k", col_order = [0.1, 0.05, 0.02], row_order = [10, 15, 20],
        #                     hue="hash", data=indata,
        #                     size=3, aspect=1.6, palette="Set2",
        #                     dodge=True, cut=0, bw=.2) #kind="violin",
        
        g = sns.FacetGrid(indata, row="measure_type", size=3, aspect=2.2, row_order=["V-measure", "Completeness", "Homogeneity"], legend_out=True)
        # sns.set(style="whitegrid", palette="muted")
        # (g.map(sns.pointplot, "class_size", "measure", "Dataset", hue_order=["ALZ_PB", "BHAM_ONT", "RC0_PB", "HUM_PB", "ZEB_PB", "ENS_PB_100k", "ENS_PB_500k", "ENS_PB_1M"], palette=sns.color_palette("Set1", 8)).despine(left=True).add_legend(title="DATASET", label_order=["ALZ_PB", "BHAM_ONT", "RC0_PB", "HUM_PB", "ZEB_PB", "ENS_PB_100k", "ENS_PB_500k", "ENS_PB_1M"]))
        g.map(sns.pointplot, "class_size", "measure", "Dataset", hue_order=["ALZ_PB", "BHAM_ONT", "RC0_PB", "HUM_PB", "ZEB_PB", "ENS_PB_100k", "ENS_PB_500k", "ENS_PB_1M"], palette=sns.color_palette("Set1", 8))

        # g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=16)
        # g.set_yticklabels(["",0,0.2,0.4,0.6,0.8,1.0])
        # g.set(yscale="log")
        # g.set(xlim=(0, 6))
        # g.set(ylim=(0.75, 1))

        # resize figure box to -> put the legend out of the figure
        # box = g.ax.get_position() # get position of figure
        # g.ax.set_position([box.x0, box.y0, box.width * 0.85, box.height]) # resize position

        # # Put a legend to the right side
        # g.ax.legend(loc='lower center', bbox_to_anchor=(1.25, 0.5), ncol=1)
        plt.legend(loc='upper left', bbox_to_anchor=(1,0.5))
        # plt.legend(loc='lower right')
        plt.tight_layout()
        # g.set_title('Clustering accuracy per class size')
        # g.set_ylabel("Performance")
        # g.set_xlabel("Class size")
        plt.savefig(args.outfile)
        plt.close()



def main(args):
    plot_V_per_expression_bin(args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--infile', type=str, help='csv infile')
    parser.add_argument('--outfile', type=str, help='Output file with results')
    args = parser.parse_args()

    # if not os.path.exists(args.outfolder):
    #     os.makedirs(args.outfolder)

    main(args)



