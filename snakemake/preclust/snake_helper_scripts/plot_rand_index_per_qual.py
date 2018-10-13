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

import argparse, os


def plot_rand_index_per_error_batch(csv_file, outfile):
    # sns.plt.clf()
    with sns.plotting_context("paper", font_scale=2.4):
        sns.set()
        indata = pd.read_csv(csv_file)
        # ARI,error_rate,nr_samples
        g = sns.pointplot(x="error_rate", y="RI", data=indata)

        # g = sns.factorplot(x="w", y="P_minimizer_shared", col="error_rate", row="k", col_order = [0.1, 0.05, 0.02], row_order = [10, 15, 20],
        #                     hue="hash", data=indata,
        #                     size=3, aspect=1.6, palette="Set2",
        #                     dodge=True, cut=0, bw=.2) #kind="violin",
        
        # g = sns.FacetGrid(data, row="Family", col="mutation_rate", size=3, aspect=1.6, row_order=["TSPY13P", "HSFY2", "DAZ2"], col_order=[0.01, 0.001, 0.0001], legend_out=True)
        # sns.set(style="whitegrid", palette="muted")
        # (g.map(sns.violinplot, "read_count", args.y_axis, "TOOL", cut=0, hue_order=["ISOCON", "ICE"], palette=sns.color_palette("muted", 2)).despine(left=True).add_legend(title="TOOL", label_order=["ISOCON", "ICE"]))
        # g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=16)
        # g.set_yticklabels(["",0,0.2,0.4,0.6,0.8,1.0])
        # g.set(yscale="log")
        # g.set(xlim=(0, 6))
        g.set(ylim=(0, 1))
        plt.legend(loc='lower right')
        plt.tight_layout()
        g.set_title('Rand index per error rate')
        g.set_ylabel("Rand index")
        g.set_xlabel("Read error rate")
        g.set_xticklabels(g.get_xticklabels(), rotation=90)
        plt.savefig(outfile)
        plt.close()




def main(args):
    plot_rand_index_per_error_batch(args.infile, args.outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--infile', type=str, help='Csv infile with results')
    parser.add_argument('--outfile', type=str, help='Output file with results')
    args = parser.parse_args()

    # if not os.path.exists(args.outfolder):
    #     os.makedirs(args.outfolder)

    main(args)




