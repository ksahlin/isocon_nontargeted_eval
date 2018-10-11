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

from sklearn.metrics.cluster import v_measure_score, completeness_score, homogeneity_score

def parse_inferred_clusters_tsv(tsv_file, args):
    infile = open(tsv_file , "r")
    clusters = {}
    for line in infile:
        cluster_id, read_acc = line.strip().split("\t")
        if args.simulated:
            read_acc = "_".join([item for item in read_acc.split("_")[:-1] ])
        elif args.ont:
            read_acc = read_acc.split("_runid")[0]
        else:
            read_acc = read_acc.split("_strand")[0] 

        clusters[read_acc] = int(cluster_id)
    return clusters


def parse_true_clusters(ref_file):
    classes = defaultdict(dict)
    class_index = 1
    class_ranges = {}
    ref_id_to_chrom = {}
    alignment_counter = defaultdict(int)
    
    prev_chrom = -1
    curr_class_id = -1
    prev_class_start = -1
    prev_class_stop = -1
    prev_read_id = ""
    unique_reads = set()
    unclassified = 0
    for read in ref_file.fetch(until_eof=True):
        unique_reads.add(read.query_name)
        if read.is_unmapped:
            unclassified += 1
            continue
        if read.is_secondary or read.is_supplementary: # deal with supplementary alignments!!
            continue
        # print(read.query_name, read.flag)
        assert prev_read_id != read.query_name

        chrom = read.reference_name
        if chrom != prev_chrom:
            curr_class_id += 1
            classes[read.query_name] = curr_class_id
            prev_chrom = chrom
            prev_class_start = read.reference_start
            prev_class_stop = read.reference_end

        else:
            read_ref_start = read.reference_start
            read_ref_end = read.reference_end
            if read_ref_start > prev_class_stop:
                curr_class_id += 1
                classes[read.query_name] = curr_class_id
                prev_class_start = read.reference_start
                prev_class_stop = read.reference_end
            else:
                classes[read.query_name] = curr_class_id
            

        prev_class_stop = max(read.reference_end, prev_class_stop)
        prev_read_id = read.query_name



        # classes[read.query_name] = int(read.reference_id)  # chrom
        ref_id_to_chrom[int(read.reference_id)] = chrom
        alignment_counter[int(read.reference_id)] += 1
        # if chrom not in class_ranges:
        #     class_ranges[chrom] = {}

        # print(chrom, read_ref_start, read_ref_end)
        # for start, stop in class_ranges[chrom]:
        #     if start <= read_ref_start and  read_ref_end <= stop:
        #         # entirly within
        #     elif

        # classes[read.query_name] =  #read.reference_name.split("|")[0]  #"|".join(read.reference_name.split("|")[:2])
    print(ref_id_to_chrom)
    print()
    print(alignment_counter)
    print()
    return classes, len(unique_reads), unclassified


def parse_true_clusters_simulated(ref_file):
    classes = defaultdict(dict)
    for read in ref_file.fetch(until_eof=True):
        classes[read.query_name] = read.reference_name.split("|")[0] # by gene id 
        # classes[read.query_name] = read.reference_name.split("|")[1] # by transcript id 
    return classes


def compute_V_measure_per_expression_bin(clusters, classes, tmp_file):
    """
        split classes dict into dicts of classes with 0-5, 6-10, 11-50 , >50 sizes. 
        Compute homog, completeness, V measure on each of these subdata sets.
    """


    class_expression_levels = {}
    inv_map = {}
    for read_acc, cl_id in classes.items():
        inv_map.setdefault(cl_id, []).append(read_acc)

    for read_acc, cl_id in classes.items():
        class_expression_levels[read_acc] = len(inv_map[cl_id])
        # print(len(inv_map[cl_id]))


    class_sorted_by_expression, cluster_sorted_by_expression = {}, {}
    print(len(clusters), len(classes))
    # not_found_id = 1000000
    clustered_but_unaligned = 0
    for read_acc in classes:
        if read_acc in clusters:
            # print(class_expression_levels[read_acc])
            if  class_expression_levels[read_acc] in class_sorted_by_expression:
                class_sorted_by_expression[ class_expression_levels[read_acc] ].append( classes[read_acc] )
                cluster_sorted_by_expression[ class_expression_levels[read_acc] ].append( clusters[read_acc] )
            else:
                class_sorted_by_expression[ class_expression_levels[read_acc] ] = [classes[read_acc] ]
                cluster_sorted_by_expression[ class_expression_levels[read_acc] ] = [ clusters[read_acc] ]
        else:
            pass

    # print(class_sorted_by_expression)
    # bin together into
    bins= [(0,5), (6,10), (11,20), (21,50), (50, 10000000000)]
    # final_out = []
    tmp_file_ = open(tmp_file, "w")
    tmp_file_.write("class_size,measure,nr_samples,measure_type\n")
    for  (l, u) in bins:
        class_list = []
        cluster_list = []
        unique_classes = set()
        for expr_level in  class_sorted_by_expression:
            print(expr_level, l, u)
            if l <= expr_level <= u:
                for class_id, cluster_id in zip(class_sorted_by_expression[expr_level], cluster_sorted_by_expression[expr_level]):
                    class_list.append(class_id)
                    cluster_list.append(cluster_id)
                    unique_classes.add(class_id)


        v_score = v_measure_score(class_list, cluster_list)
        compl_score = completeness_score(class_list, cluster_list)
        homog_score = homogeneity_score(class_list, cluster_list)
        print("Bin size: {0}-{1}:".format(l, u))
        print("Nr samples:", len(unique_classes))
        print("V:", v_score, "Completeness:", compl_score, "Homogeneity:", homog_score)
        if l != 50:
            tmp_file_.write("{0},{1},{2},{3}\n".format(str(l)+ "-" + str(u), v_score, len(unique_classes), "V-measure" ))
            tmp_file_.write("{0},{1},{2},{3}\n".format(str(l)+ "-" + str(u), compl_score,  len(unique_classes), "Completeness" ))
            tmp_file_.write("{0},{1},{2},{3}\n".format(str(l)+ "-" + str(u), homog_score, len(unique_classes), "Homogeneity" ))
        else:
            tmp_file_.write("{0},{1},{2},{3}\n".format(">" + str(l), v_score, len(unique_classes), "V-measure" ))
            tmp_file_.write("{0},{1},{2},{3}\n".format(">" + str(l), compl_score,  len(unique_classes), "Completeness" ))
            tmp_file_.write("{0},{1},{2},{3}\n".format(">" + str(l), homog_score, len(unique_classes), "Homogeneity" ))

    tmp_file_.close
        # final_out.append( ((l,u), v_score, compl_score, homog_score ) )

    return tmp_file_.name



def plot_V_per_expression_bin(csv_file, outfile):
    # sns.plt.clf()
    with sns.plotting_context("paper", font_scale=2.4):
        sns.set()
        indata = pd.read_csv(csv_file)
        g = sns.pointplot(x="class_size", y="measure", data=indata, hue="measure_type")

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
        g.set_title('Clustering accuracy per class size')
        g.set_ylabel("Performance")
        g.set_xlabel("Class size")
        plt.savefig(outfile)
        plt.close()



def main(args):

    clusters = parse_inferred_clusters_tsv(args.clusters, args)
    if args.simulated:
        ref_file = pysam.AlignmentFile(args.classes, "r", check_sq=False)
        classes  = parse_true_clusters_simulated(ref_file)
        tot_nr_reads = len(classes) # by simulation we know classes of all reads, they are therefore the same number.
    else:
        ref_file = pysam.AlignmentFile(args.classes, "rb", check_sq=False)
        classes, tot_nr_reads, unclassified  = parse_true_clusters(ref_file)

    tmp_file_ = compute_V_measure_per_expression_bin(clusters, classes, args.outfile+"_tmp.csv")
    plot_V_per_expression_bin(tmp_file_, args.outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--clusters', type=str, help='Inferred clusters (tsv file)')
    parser.add_argument('--classes', type=str, help='A sorted and indexed bam file.')
    parser.add_argument('--simulated', action="store_true", help='Simulated data, we can simply read correct classes from the ref field.')
    parser.add_argument('--ont', action="store_true", help='ONT data, parsing accessions differently.')
    parser.add_argument('--outfile', type=str, help='Output file with results')
    args = parser.parse_args()

    # if not os.path.exists(args.outfolder):
    #     os.makedirs(args.outfolder)

    main(args)



