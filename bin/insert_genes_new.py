#!/usr/bin/env python
"""
#refactored insert_gene.py script
#instead of relying on RGI for the BED file, just make the BED file as genes are inserted
#this makes this workflow broadly applicable to non-amr inserts
#also cuz the BED file made from RGI contains false positives so you cannot map with 100% accuracy afterwards.
"""
from Bio import SeqIO
import random
import itertools
import argparse

if __name__ == "__main__":
    class Inserts():
        def __init__(self, name, contig, start, length) -> None:
            self.name = name
            self.contig = contig
            self.origInsertIndex = start
            self.currentInsertIndex = self.origInsertIndex
            self.length = length
            self.end = self.origInsertIndex + self.length
            self.positions = [x for x in range(self.origInsertIndex, self.end)]
        
        #updates the position array of this insert using the new insert.
        def Update(self, contigName, insertStartPos, insertLength):
            if (contigName == self.contig): #same contig, may need to update
                if (insertStartPos < (self.end)): #inserted into a position that already had an existing insert
                    if (insertStartPos >= self.currentInsertIndex):
                        #inserted into existing insert
                        for i in range(len(self.positions)):
                            if self.positions[i] >= insertStartPos:
                                self.positions[i] = self.positions[i] + insertLength
                        self.currentInsertIndex = self.currentInsertIndex
                        self.end = self.end + insertLength 
                    else:
                        #inserted into a location before this existing insert
                        self.positions = [i + insertLength for i in self.positions]
                        self.currentInsertIndex = self.currentInsertIndex + insertLength
                        self.end = self.end + insertLength 
                else:
                    #inserted into a location after this insert. do not give a shit
                    self.positions = self.positions
                    self.currentInsertIndex = self.currentInsertIndex
                    self.end = self.end

            else:
                self.positions = self.positions
                self.currentInsertIndex = self.currentInsertIndex
                self.end = self.end

        #consecutive number to ranges
        def ranges(self, i):
            for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
                b = list(b)
                yield b[0][1], b[-1][1]

        def toTSV(self):
            range = list(self.ranges(self.positions))
            return "{}\t{}\t{}\t{}\t{}".format(self.name, self.contig, self.length, self.origInsertIndex, range)
        
        def toBED(self):
            range = list(self.ranges(self.positions))
            out = []
            for r in range:
                contig = self.contig.replace(":","|").replace("-","|")
                start = r[0]
                stop = r[1] + 1
                aro = self.name.split("|")[0].replace("ARO:","")
                out.append("{}\t{}\t{}\t{}\t{}\t{}".format(contig, start, stop, aro, "9999", "+"))
                
            return out
            #NZ_CP007772.1|1403213|1411695_1	61336	62487	3003796	1000	+



    parser = argparse.ArgumentParser("Randomly insert genes into genome")

    parser.add_argument("-i", "--input_genome", required=True,
                        help='Path to fasta containing genome')
    parser.add_argument("-g", "--insert_genes", required=True,
                        help='Path to fasta containing genes to insert')
    parser.add_argument("-o", "--output", required=True,
                        help='Path to output modified genome')

    args = parser.parse_args()


    #source = "source.fasta" 
    #insert = "insert.fasta" 
    #output = "out.tsv"
    source = args.input_genome
    insert = args.insert_genes
    output = args.output

    genome_contigs = {}
    count = 0
    with open(source) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            genome_contigs[count] = [record.id, str(record.seq)]
            count = count + 1

    inserts = []
    with open(insert) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            contig_ix = random.randint(0, len(genome_contigs) - 1) #get a contig number
            contig_seq = genome_contigs[contig_ix][1] #get the contig sequence

            insert_name = record.id
            insert_seq = str(record.seq)
            insert_pos = random.randint(0, len(contig_seq) - 1) #get the insert position
            insert_len = len(insert_seq)
            insertData = Inserts(insert_name, genome_contigs[contig_ix][0], insert_pos, insert_len)

            contig_seq = contig_seq[:insert_pos] + insert_seq + contig_seq[insert_pos: ]
            genome_contigs[contig_ix][1] = contig_seq

            for insert in inserts:
                insert.Update(genome_contigs[contig_ix][0], insert_pos, insert_len)
            
            inserts.append(insertData)

    with open(output, 'w') as fh:
        for contig in list(genome_contigs.values()):
            id = ">" + contig[0]
            seq = contig[1]
            fh.write(id + "\n")
            fh.write(seq + "\n")

    with open (output + ".inserts.tsv", 'w') as fh:
        for i in inserts:
            for b in i.toBED():
                fh.write(b + "\n")
