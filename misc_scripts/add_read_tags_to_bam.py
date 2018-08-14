import argparse, os
import pysam

from collections import defaultdict
import errno

def mkdir_p(path):
    print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def add_tags(args):
    infile = pysam.AlignmentFile(args.inbam, "rb")
    outfile = pysam.AlignmentFile(args.outbam, "wb", template=infile)
    for read in infile.fetch(until_eof=True):
        read.set_tag("zm",4194413)
        read.set_tag("RG","fd5bcde2")
        outfile.write(read)

    infile.close()
    outfile.close()

def main(args):
    add_tags(args)
 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--inbam', type=str, help='Inferred clusters (tsv file)')
    parser.add_argument('--outbam', type=str, help='A sorted and indexed bam file.')
    args = parser.parse_args()
    
    path_, file_prefix = os.path.split(args.outbam)
    mkdir_p(path_)
    main(args)


