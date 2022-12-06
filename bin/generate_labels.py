#!/usr/bin/env python

import sys
import pysam
import argparse
import pandas as pd

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Generate labels from bam and bed")
    parser.add_argument("--bed", required=True,
                        help="Path to bed containing labels (listed in name column)")
    parser.add_argument("--bam", required=True,
                        help="Path to metagenome bam file")
    args = parser.parse_args()

    amr_genes = pd.read_csv(args.bed, sep='\t', names=['chrom', 'chromStart',
                                                       'chromEnd', 'name',
                                                       'score', 'strand'])
    labels = {'Read Name': [], 'ARO': []}
    unique_amr_genes = amr_genes['name'].unique()
    for amr_gene in unique_amr_genes:
        bed = amr_genes[amr_genes['name'] == amr_gene]

        bed.to_csv("temp.bed", header=False, sep='\t', index=False)
        amr_reads = pysam.view(args.bam, "-ML", "temp.bed")

        for read in amr_reads.strip().split('\n'):
            read_data = read.split('\t')
            labels['Read Name'].append(read_data[0])
            labels['ARO'].append(str(amr_gene))

    labels = pd.DataFrame(labels)
    labels = labels.drop_duplicates()
    labels['ARO'] = labels.groupby('Read Name')['ARO'].transform(lambda x: ";".join(x))
    labels = labels.drop_duplicates()
    labels.to_csv(sys.stdout, sep='\t', index=False)
