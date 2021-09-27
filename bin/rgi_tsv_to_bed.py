#!/usr/bin/env python

import sys
import csv
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Convert RGI tsv to a .bed file")

    parser.add_argument("-i", "--input_rgi_tsv", required=True,
                        help='Path to RGI tsv output file ($prefix.txt)')
    parser.add_argument("-o", "--output_bed", required=True,
                        help='Path to output bed file')

    args = parser.parse_args()


    bedfile = open(args.output_bed, 'w')
    bed_writer = csv.DictWriter(bedfile, delimiter='\t',
                                fieldnames=['chrom', 'chromStart', 'chromEnd',
                                            'name', 'score', 'strand'])

    with open(args.input_rgi_tsv) as rgi_csv:
        reader = csv.DictReader(rgi_csv, delimiter='\t')
        for row in reader:
            contig_name = "_".join(row['Contig'].split('_')[:-1])

            if row['Cut_Off'] == 'Perfect':
                score = 1000
            elif row['Cut_Off'] == 'Strict':
                score = 500
            bed_writer.writerow({'chrom': contig_name,
                                 'chromStart': row['Start'],
                                 'chromEnd': row['Stop'],
                                 'name': row['ARO'],
                                 'score': score,
                                 'strand': row['Orientation']})
