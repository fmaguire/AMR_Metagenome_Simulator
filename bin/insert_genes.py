#!/usr/bin/env python

from Bio import SeqIO
import random
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Randomly insert genes into genome")

    parser.add_argument("-i", "--input_genome", required=True,
                        help='Path to fasta containing genome')
    parser.add_argument("-g", "--insert_genes", required=True,
                        help='Path to fasta containing genes to insert')
    parser.add_argument("-o", "--output", required=True,
                        help='Path to output modified genome')

    args = parser.parse_args()

    genes_to_insert = []
    with open(args.insert_genes) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            genes_to_insert.append(record)


    genome_contigs = []
    with open(args.input_genome) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            genome_contigs.append(record)

    for gene in genes_to_insert:
        contig_ix = random.randint(0, len(genome_contigs) - 1)
        contig_seq = genome_contigs[contig_ix].seq
        insert_pos = random.randint(0, len(contig_seq) - 1)
        contig_seq = contig_seq[:insert_pos] + gene.seq + contig_seq[insert_pos: ]
        genome_contigs[contig_ix].seq = contig_seq


    with open(args.output, 'w') as fh:
        SeqIO.write(genome_contigs, fh, 'fasta-2line')
