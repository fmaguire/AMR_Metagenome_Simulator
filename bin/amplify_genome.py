#!/usr/bin/env python

import argparse
import csv
import copy
from Bio import SeqIO

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Amplify genomes")
    parser.add_argument("-i", "--input_genome", required=True,
                        help="Genome to amplify")
    parser.add_argument("-a", "--input_annotation", required=True,
                        help="BED file of annotations for genome")
    parser.add_argument("-c", "--copy_number", required=True, type=int,
                        help="Number of copies to make for genome")
    parser.add_argument("-o", "--output_prefix", required=True, type=str,
                        help="Prefix string for output")


    args = parser.parse_args()

    with open(args.input_genome) as fh:
        genome = [record for record in SeqIO.parse(fh, 'fasta')]


    with open(args.input_annotation) as fh:
        reader = csv.reader(fh, delimiter='\t')
        annotations = [record for record in reader]

    output_genome = open(f"{args.output_prefix}_amplified.fna", 'w')
    output_bed_fh = open(f"{args.output_prefix}_amr_amplified.bed", 'w', newline='')
    output_bed = csv.writer(output_bed_fh, delimiter='\t')

    for copy_ix in range(args.copy_number):
        for record in genome:
            new_record = copy.deepcopy(record)
            new_record.id = new_record.id + f"_{copy_ix}"
            SeqIO.write(new_record, output_genome, 'fasta')
        for annotation in annotations:
            new_annotation = copy.deepcopy(annotation)
            new_annotation[0] = annotation[0] + f"_{copy_ix}"
            output_bed.writerow(new_annotation)

    output_genome.close()
    output_bed_fh.close()


