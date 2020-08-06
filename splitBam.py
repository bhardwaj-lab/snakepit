#!/usr/bin/env python

import pysam
import csv
import argparse
import os
import sys
import textwrap

def parse_args(defaults=None):
    """
    Parse arguments from the command line.
    """
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False
    )

    # Workflow options
    parser.add_argument("-i",
                          dest="input_bam",
                          metavar="STR",
                          help="Input BAM",
                          type=str,
                          required=True)

    parser.add_argument("-c",
                          dest="clusters_tsv",
                          metavar="STR",
                          help="TSV file with barcodes '\\t' cluster_name",
                          type=str,
                          required=True)
    parser.add_argument("-t",
                          dest="barcodeTag",
                          metavar="STR",
                          help="tag for barcodes in the BAM file",
                          type=str,
                          default="BC",
                          required=True)

    parser.add_argument("-o",
                          dest="outprefix",
                          metavar="STR",
                          help="output prefix",
                          type=str,
                          default=".")

    parser.add_argument("-h", "--help",
                         action="help",
                         help="show this help message and exit")

    return parser


def main():

    # get command line arguments
    parser = parse_args()
    args = parser.parse_args()
    tag = args.barcodeTag

    cluster_dict = {}
    with open(args.clusters_tsv) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        #skip header
        #header = next(csv_reader)
        for row in csv_reader:
            cluster_dict[row[0]] = row[1]

    clusters = set(x for x in cluster_dict.values())


    fin = pysam.AlignmentFile(args.input_bam, "rb")

    # open the number of bam files as the same number of clusters, and map the out file handler to the cluster id, write to a bam with wb
    fouts_dict = {}
    for cluster in clusters:
        fout = pysam.AlignmentFile(args.outprefix + cluster + ".bam", "wb", template = fin)
        fouts_dict[cluster] = fout

    for read in fin:
        tags = read.tags
        CB_list = [ x for x in tags if x[0] == tag]
        if CB_list:
            cell_barcode = CB_list[0][1]
        # the bam files may contain reads not in the final clustered barcodes
        # will be None if the barcode is not in the clusters.csv file
        else:
            continue
        cluster_id = cluster_dict.get(cell_barcode)
        if cluster_id:
            fouts_dict[cluster_id].write(read)

    ## do not forget to close the files
    fin.close()
    for fout in fouts_dict.values():
        fout.close()


if __name__ == "__main__":
    main()
