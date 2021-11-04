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

    # number of empty chunks
    if args.number_of_chunks > len(insertion_records) + 1:
        empty_samples = args.number_of_chunks - len(insertion_records)
        number_of_chunks = len(insertion_records) - 1
    else:
        number_of_chunks = args.number_of_chunks
        empty_samples = 0


    print(number_of_chunks)
    split_points = np.random.choice(len(insertion_records) - 2,
                                     number_of_chunks - 1,
                                     replace=False) + 1
    split_points.sort()
    split_records = np.split(np.array(insertion_records,
                                      dtype='object'),
                             split_points)

    # add empty chunks
    for i in range(empty_samples):
        split_records.append([])

    random.shuffle(split_records)
    for chunk_ix, chunk in enumerate(split_records):
        with open(f"insert_chunk_{chunk_ix}.fna", 'w') as out_fh:
            SeqIO.write(chunk, out_fh, 'fasta-2line')
