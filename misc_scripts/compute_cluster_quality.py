import networkx as nx
import argparse, os
import pysam

from collections import defaultdict

from sklearn.metrics.cluster import v_measure_score, completeness_score, homogeneity_score

def parse_inferred_clusters_tsv(tsv_file):
    infile = open(tsv_file , "r")
    clusters = {}
    for line in infile:
        cluster_id, read_acc = line.strip().split("\t")
        read_acc = "_".join([item for item in read_acc.split("_")[:-1] ])
        clusters[read_acc] = int(cluster_id)
    return clusters

# def parse_true_clusters(sam_file):
#     ccs_dict = defaultdict(dict)
#     for read in ccs_file.fetch(until_eof=True):
#         ccs_read = CCS(read.query_name, read.query_alignment_sequence, read.query_qualities, read.get_tag("np"))
#         ccs_dict[read.query_name] = ccs_read
#         # ccs_dict[read.query_name]["seq"] = read.query_alignment_sequence
#         # print(read.query_qualities)
#         # sys.exit()
#         # ccs_dict[read.query_name]["qual"] = read.query_qualities
#         # ccs_dict[read.query_name]["np"] = read.get_tag("np")
#         assert len(read.query_alignment_sequence) == len(read.query_qualities)
        
#         # if ccs_read.np > 10:
#         #     print(ccs_read.np, ccs_read.positions_with_p_error_higher_than(0.01))
#     return ccs_dict


def parse_true_clusters_simulated(ref_file):
    classes = defaultdict(dict)
    for read in ref_file.fetch(until_eof=True):
        classes[read.query_name] = read.reference_name.split("|")[0]
    return classes

def compute_V_measure(clusters, classes):
    class_list, cluster_list = [], []
    print(len(clusters), len(classes))
    not_found_id = 1000000
    for read in classes:
        class_list.append( classes[read] )
        if read not in clusters:
            cluster_list.append(not_found_id)
            not_found_id += 1
        else:
            cluster_list.append( clusters[read] )


    v_score = v_measure_score(class_list, cluster_list)
    compl_score = completeness_score(class_list, cluster_list)
    homog_score = homogeneity_score(class_list, cluster_list)
    print(v_score, compl_score, homog_score)


def main(args):

    clusters = parse_inferred_clusters_tsv(args.clusters)
    ref_file = pysam.AlignmentFile(args.classes, "r", check_sq=False)
    if args.simulated:
        classes  = parse_true_clusters_simulated(ref_file)
    else:
        truth  = parse_true_clusters(args.classes)

    value = compute_V_measure(clusters, classes)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--clusters', type=str, help='Inferred clusters (tsv file)')
    parser.add_argument('--classes', type=str, help='A sam file.')
    parser.add_argument('--simulated', action="store_true", help='Simulated data, we can simply read correct classes from the ref field.')
    parser.add_argument('--outfolder', type=str, help='Output path of results')
    args = parser.parse_args()

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)