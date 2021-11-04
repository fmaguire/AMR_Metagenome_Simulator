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
            if row['Model_type'] == 'rRNA gene variant model':
                contig_name = row['Contig']
            else:
                contig_name = "_".join(row['Contig'].split('_')[:-1])

            # append snp information to ARO as name is only field that
            # can include this in bed format
            if row['Model_type'] in ['rRNA gene variant model',
                                     'protein variant model']:
                name = f"{row['ARO']}:{row['SNPs_in_Best_Hit_ARO'].replace(' ', '')}"
            else:
                name = row['ARO']

            if row['Cut_Off'] == 'Perfect':
                score = 1000
            elif row['Cut_Off'] == 'Strict':
                score = 500

            bed_writer.writerow({'chrom': contig_name,
                                 'chromStart': row['Start'],
                                 'chromEnd': row['Stop'],
                                 'name': name,
                                 'score': score,
                                 'strand': row['Orientation']})
