#!/usr/bin/env python

import argparse
import random
from Bio import SeqIO
import numpy as np


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Randomly split a fasta file into N chunks")
    parser.add_argument("-i", "--input_fasta", required=True,
                        help='Path to fasta file to split')
    parser.add_argument("-n", "--number_of_chunks", required=True, type=int,
                         help="Number of pieces to randomly split fasta into")
    args = parser.parse_args()

    insertion_records = []
    with open(args.input_fasta) as in_fh:
        for record in SeqIO.parse(in_fh, 'fasta'):
            insertion_records.append(record)

    # shuffle list
    random.shuffle(insertion_records)

    split_points = np.random.choice(len(insertion_records) - 2,
                                     args.number_of_chunks - 1,
                                     replace=False) + 1
    split_points.sort()
    split_records = np.split(np.array(insertion_records,
                                      dtype='object'),
                             split_points)

    for chunk_ix, chunk in enumerate(split_records):
        with open(f"insert_chunk_{chunk_ix}.fna", 'w') as out_fh:
            SeqIO.write(chunk, out_fh, 'fasta')
