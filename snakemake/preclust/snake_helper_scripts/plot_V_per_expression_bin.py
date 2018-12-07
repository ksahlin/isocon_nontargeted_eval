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
        sns.set_style("whitegrid")

        # sns.set()
        indata = pd.read_csv(args.infile)
        # g = sns.pointplot(x="class_size", y="measure", data=indata, hue="measure_type")

        # g = sns.factorplot(x="w", y="P_minimizer_shared", col="error_rate", row="k", col_order = [0.1, 0.05, 0.02], row_order = [10, 15, 20],
        #                     hue="hash", data=indata,
        #                     size=3, aspect=1.6, palette="Set2",
        #                     dodge=True, cut=0, bw=.2) #kind="violin",
        
        g = sns.FacetGrid(indata, col="measure_type", size=4, aspect=1.6, col_order=["Completeness", "Homogeneity"], sharey=False, legend_out=True)
        # sns.set(style="whitegrid", palette="muted")
        (g.map(sns.pointplot, "class_size", "measure", "dataset", hue_order=["ALZ", "RC0", "HUM", "ZEB", "SIM-100k", "SIM-500k", "SIM-1000k", "ONT"], palette=sns.color_palette("colorblind", 8)).despine(left=True).add_legend(title="DATASET")) # .despine(left=True).add_legend(title="DATASET", label_order=["ALZ_PB", "BHAM_ONT", "RC0_PB", "HUM_PB", "ZEB_PB", "ENS_PB_100k", "ENS_PB_500k", "ENS_PB_1M"]))
        # g.map(sns.pointplot, "class_size", "measure", "Dataset", hue_order=["ALZ_PB", "BHAM_ONT", "RC0_PB", "HUM_PB", "ZEB_PB", "ENS_PB_100k", "ENS_PB_500k", "ENS_PB_1M"], palette=sns.color_palette("Set1", 8))
        axes = g.axes.flatten()
        axes[0].set_title("")
        axes[1].set_title("")
        axes[0].set_ylabel("Completeness")
        axes[1].set_ylabel("Homogeneity")
        for ax in axes:
            ax.set_xlabel("Class size")
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
        # plt.legend(loc='upper left', bbox_to_anchor=(1,0.5))
        # plt.legend(loc='lower right')
        # plt.tight_layout()
        # g.set_title('Clustering accuracy per class size')
        # g.set_ylabel("Performance")
        # g.set_xlabel("Class size")
        plt.savefig(args.outfile)
        plt.close()

def plot_V_per_expression_bin2(args):
    with sns.plotting_context("paper", font_scale=2.4):
        # sns.set(style="whitegrid")
        rc={'axes.labelsize': 24.0, 'legend.fontsize': 28.0, 'axes.titlesize': 28, 'xtick.labelsize': 24.0, 'ytick.labelsize' : 24.0}
        sns.set(rc=rc, style="whitegrid")

        indata = pd.read_csv(args.infile)
 
        # print(indata.sort_values(by="sort_by_column", axis=0, na_position="first"))
        # Make the PairGrid
        
        g = sns.PairGrid(indata, hue="class_size", hue_order=["1-5", "6-10", "11-20", "21-50", ">50"],
                         x_vars=indata.columns[2:], y_vars=["dataset"],
                         height=10, aspect=.5)
        # Draw a dot plot using the stripplot function
        g.map(sns.stripplot, size=20, orient="h", jitter=0.1, #palette="ch:s=1,r=-.1,h=1_r",
               linewidth=1, edgecolor="w", alpha=.8)
        # y_vars = ["ALZ_PB", "RC0_PB", "HUM_PB", "ZEB_PB", "BHAM_ONT", "ENS_100", "ENS_500", "ENS_1M"],
        # Use the same x axis limits on all columns and add better labels
        # g.set(xlim=[(0, 1), (0, 1),(0, 1), (0, 100) ], xlabel=["", "", "", "", ""], ylabel=["", "", "", "", ""])
        g.set(xlabel="", ylabel="")
        g.axes[0,0].set_xlim((0,5000))
        # g.axes[0,1].set_xlim((0.75,1))
        # g.axes[0,2].set_xlim((0.75,1))
        # g.axes[0,3].set_xlim((0.9,1))
        g.set(xticks=[0.8,0.9,1.0])
        g.axes[0,0].set_xticks((0, 5000, 10000))
        g.axes[0,1].set_xticks((0.8, 0.9, 1.0))
        g.axes[0,2].set_xticks((0.8, 0.9, 1.0))
        g.axes[0,3].set_xticks((0.8,0.9,1.0))

        # Use semantically meaningful titles for the columns
        titles = ["#classes", "V-measure", "Completeness", "Homogeneity"]

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

        plt.savefig(args.outfile )
        plt.clf()

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



