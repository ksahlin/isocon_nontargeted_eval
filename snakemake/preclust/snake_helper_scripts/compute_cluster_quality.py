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

# def compute_V_measure(clusters, classes):
#     class_list, cluster_list = [], []
#     print(len(clusters), len(classes))
#     not_found_id = 1000000
#     for read in classes:
#         class_list.append( classes[read] )
#         if read not in clusters:
#             cluster_list.append(not_found_id)
#             not_found_id += 1
#         else:
#             cluster_list.append( clusters[read] )


#     v_score = v_measure_score(class_list, cluster_list)
#     compl_score = completeness_score(class_list, cluster_list)
#     homog_score = homogeneity_score(class_list, cluster_list)
#     print(v_score, compl_score, homog_score)


def compute_V_measure(clusters, classes):
    class_list, cluster_list = [], []
    print(len(clusters), len(classes))
    # not_found_id = 1000000
    clustered_but_unaligned = 0
    for read in clusters:
        if read in classes:
            class_list.append( classes[read] )
            cluster_list.append( clusters[read] )
        else:
            # print("Read was clustered but unaligned:", read)
            clustered_but_unaligned +=1

    # added the unprocessed reads to the measure
    not_clustered = set(classes.keys()) - set(clusters.keys())
    highest_cluster_id = max(clusters.values())
    highest_cluster_id += 1
    for read in not_clustered:
        class_list.append( classes[read] )
        cluster_list.append( highest_cluster_id )
        highest_cluster_id += 1

    v_score = v_measure_score(class_list, cluster_list)
    compl_score = completeness_score(class_list, cluster_list)
    homog_score = homogeneity_score(class_list, cluster_list)
    print("Not inluded in clustering but aligned:", len(not_clustered))
    print("V:", v_score, "Completeness:", compl_score, "Homogeneity:", homog_score)
    print("Nr reads clustered but unaligned (i.e., no class and excluded from V-veasure): ", clustered_but_unaligned)
    return v_score, compl_score, homog_score, clustered_but_unaligned



def compute_V_measure_non_singletons(clusters, classes):
    cluster_dict = {}
    for read_acc, cl_id in clusters.items():
        if cl_id not in cluster_dict:
            cluster_dict[cl_id] = [read_acc]
        else:
            cluster_dict[cl_id].append(read_acc)

    nontrivial_clustered_reads = []
    for cl_id in cluster_dict:
        if len(cluster_dict[cl_id]) <= 1:
            continue
        else:
            for read in cluster_dict[cl_id]:
                nontrivial_clustered_reads.append(read)

    class_list, cluster_list = [], []
    print(len(nontrivial_clustered_reads), len(clusters), len(classes))
    # not_found_id = 1000000
    clustered_but_unaligned = 0
    for read in nontrivial_clustered_reads:
        if read in classes:
            class_list.append( classes[read] )
            cluster_list.append( clusters[read] )
        else:
            # print("Read was clustered but unaligned:", read)
            clustered_but_unaligned +=1


    v_score = v_measure_score(class_list, cluster_list)
    compl_score = completeness_score(class_list, cluster_list)
    homog_score = homogeneity_score(class_list, cluster_list)
    print("NONTRIVIAL CLUSTERS: V:", v_score, "Completeness:", compl_score, "Homogeneity:", homog_score)
    print("NONTRIVIAL CLUSTERS: Nr reads clustered but unaligned (i.e., no class and excluded from V-veasure): ", clustered_but_unaligned)
    return v_score, compl_score, homog_score, clustered_but_unaligned

# def singleton_stats():

def get_cluster_information(clusters, classes):


    # class distribution
    class_dict = {}
    for read_acc, class_id in classes.items():
        if class_id not in class_dict:
            class_dict[class_id] = [read_acc]
        else:
            class_dict[class_id].append(read_acc)

    total_nr_classes = len(class_dict)
    class_distribution = sorted([len(cl) for cl in class_dict.values()])
    singleton_classes = set([acc_list[0] for cl_id, acc_list in class_dict.items() if len(acc_list) == 1 ])
    min_class_size = min(class_distribution)
    max_class_size = max(class_distribution)
    mean_class_size = sum(class_distribution) / float(len(class_distribution))
    median_class_size = class_distribution[ int(len(class_distribution)/2) ] if len(class_distribution) % 2 == 1 else ( class_distribution[ int(len(class_distribution)/2) ] + class_distribution[ int(len(class_distribution)/2) -1] ) / 2.0

    # cluster distribution
    cluster_dict = {}
    for read_acc, cl_id in clusters.items():
        if cl_id not in cluster_dict:
            cluster_dict[cl_id] = [read_acc]
        else:
            cluster_dict[cl_id].append(read_acc)

    total_nr_clusters = len(cluster_dict) 
    cluster_distribution = sorted([len(cl) for cl in cluster_dict.values()])
    
    # in case unclustered reads are missing from output (as for isoseq3)
    omitted_from_output_singletons = set(classes.keys()) - set(clusters.keys())
    cluster_distribution = [1 for i in range(len(omitted_from_output_singletons))] + cluster_distribution

    singleton_clusters = set([acc_list[0] for cl_id, acc_list in cluster_dict.items() if len(acc_list) == 1 ])
    min_cluster_size = min(cluster_distribution)
    max_cluster_size = max(cluster_distribution)
    mean_cluster_size = sum(cluster_distribution) / float(len(cluster_distribution))
    median_cluster_size = cluster_distribution[ int(len(cluster_distribution)/2) ] if len(cluster_distribution) % 2 == 1 else ( cluster_distribution[ int(len(cluster_distribution)/2) ] + cluster_distribution[ int(len(cluster_distribution)/2) -1] ) / 2.0

    unaligned_but_nontrivially_clustered = set(clusters.keys()) - singleton_clusters - set(classes.keys())

    # not_considered = set([read for read in classes if read not in clusters ])

    not_clustered_classes = defaultdict(int)
    clustered_classes = defaultdict(int)
    reads_not_clustered = defaultdict(list)
    for read in clusters:
        if read in classes:
            class_id = classes[read]
        else:
            class_id = "unaligned"
        
        if read in singleton_clusters:
            not_clustered_classes[class_id] += 1
            reads_not_clustered[class_id].append(read)
        else:
            clustered_classes[class_id] += 1

    print()
    print("Not clustered:")
    for cl_id in reads_not_clustered:
        print("{0}\t{1}".format(cl_id, reads_not_clustered[cl_id] ))

    print()

    print("UNCLUSTERED:", "Tot classes:", len(not_clustered_classes), sorted(not_clustered_classes.items(), key=lambda x: x[1], reverse=True))
    print("CLUSTERED:", "Tot classes:", len(clustered_classes), sorted(clustered_classes.items(), key=lambda x: x[1], reverse=True))
    print("MIXED:", "Tot classes containing both:", len( set(clustered_classes.keys()) & set(not_clustered_classes.keys())))
    print("Total number of classes (unique gene ID):", total_nr_classes)
    return total_nr_classes, len(singleton_classes), min_class_size, max_class_size, mean_class_size, median_class_size, total_nr_clusters, len(singleton_clusters), min_cluster_size, max_cluster_size, mean_cluster_size, median_cluster_size, len(unaligned_but_nontrivially_clustered) 



def main(args):

    clusters = parse_inferred_clusters_tsv(args.clusters, args)
    if args.simulated:
        ref_file = pysam.AlignmentFile(args.classes, "r", check_sq=False)
        classes  = parse_true_clusters_simulated(ref_file)
        tot_nr_reads = len(classes) # by simulation we know classes of all reads, they are therefore the same number.
    else:
        ref_file = pysam.AlignmentFile(args.classes, "rb", check_sq=False)
        classes, tot_nr_reads, unclassified  = parse_true_clusters(ref_file)

    v_score, compl_score, homog_score, clustered_but_unaligned = compute_V_measure(clusters, classes)
    NT_v_score, NT_compl_score, NT_homog_score, _ = compute_V_measure_non_singletons(clusters, classes)
    total_nr_classes, singleton_classes, min_class_size, max_class_size, mean_class_size, median_class_size, total_nr_clusters, singleton_clusters, min_cluster_size, max_cluster_size, mean_cluster_size, median_cluster_size, unaligned_but_nontrivially_clustered  = get_cluster_information(clusters, classes)
    tot_nr_reads_included_inclustering = len(clusters)

    outfile = open(args.outfile, "w")

    outfile.write("CLASSES\n")

    #reads, unaligned, classes, singleton, min,max, mean,median

    outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format("tot_nr_reads", "unclassified", "total_nr_classes", "singleton_classes", "min_class_size", "max_class_size", "mean_class_size", "median_class_size"))
    outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(tot_nr_reads, unclassified, total_nr_classes, singleton_classes, min_class_size, max_class_size, mean_class_size, median_class_size))


    # Reads_nontrivially_clustered_(%), Singletons_(%), Reads_Nontrivially_clustered_but_unaligned, V, c,h ,V_nt, c_nt,h_nt,  non_singleton_clusters, min, max, median, mean

    Reads_nontrivially_clustered_percent = round(1.0 - float(singleton_clusters)/tot_nr_reads, 2)
    Singletons_percent = round(1.0 - Reads_nontrivially_clustered_percent)
    Reads_Nontrivially_clustered_but_unaligned = unaligned_but_nontrivially_clustered
    V, c,h = v_score, compl_score, homog_score
    V_nt, c_nt, h_nt = NT_v_score, NT_compl_score, NT_homog_score
    non_singleton_clusters = total_nr_clusters - singleton_clusters
    min_, max_, median, mean = min_cluster_size, max_cluster_size, mean_cluster_size, median_cluster_size


    outfile.write("CLUSTERS\n")
    outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n".format("Reads_nontrivially_clustered_percent", "Singletons_percent", "Reads_Nontrivially_clustered_but_unaligned", \
                                                                                    "V", "c","h", "V_nt", "c_nt", "h_nt", "non_singleton_clusters", "min_", "max_", "median", "mean" ))
    outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n".format(Reads_nontrivially_clustered_percent, Singletons_percent, Reads_Nontrivially_clustered_but_unaligned, \
                                                                                    V, c,h, V_nt, c_nt,h_nt, non_singleton_clusters, min_, max_, median, mean))


    # outfile.write("ALL CLUSTERS\n")
    # outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format("v_score", "compl_score", "homog_score", "clustered_but_unaligned", \
    #                                                              "total_nr_clusters", "singleton_clusters", "min_cluster_size", "max_cluster_size", "mean_cluster_size", "median_cluster_size", "unaligned_but_nontrivially_clustered"))
    # outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(v_score, compl_score, homog_score, clustered_but_unaligned, \
    #                                                             total_nr_clusters, singleton_clusters, min_cluster_size, max_cluster_size, mean_cluster_size, median_cluster_size, unaligned_but_nontrivially_clustered))

    # outfile.write("NONTRIVIAL CLUSTERS\n")
    # outfile.write("{0}\t{1}\t{2}\n".format("v_score", "compl_score", "homog_score"))
    # outfile.write("{0}\t{1}\t{2}\n".format(NT_v_score, NT_compl_score, NT_homog_score))
    outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--clusters', type=str, help='Inferred clusters (tsv file)')
    parser.add_argument('--classes', type=str, help='A sorted and indexed bam file.')
    parser.add_argument('--simulated', action="store_true", help='Simulated data, we can simply read correct classes from the ref field.')
    parser.add_argument('--outfile', type=str, help='Output file with results')
    args = parser.parse_args()

    # if not os.path.exists(args.outfolder):
    #     os.makedirs(args.outfolder)

    main(args)