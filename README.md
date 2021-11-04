# AMR Metagenome Simulator

Simple AMR metagenome simulator implemented in nextflow. 

This takes a spreadsheet of input fasta files (corresponding to genomes/plasmids etc) and specified relative copy number, as well as a fasta file containing genes to randomly insert into those genomes.

1. Splits fasta containing AMR genes to insert. Results in as many chunks as there are genomes each with a random number of AMR genes from the insertion file (or no genes) so that every insertion gene is included 1 time.
2. Randomly insert these AMR genes from chunks into the input genomes (caveat: could destroy/split an existing AMR gene)
3. Annotate the input genomes with inserted AMR genes using RGI and convert RGI output to a BED file.
4. Copy each input genome the specified number of times (from input CSV) to get relative abundance (note: this could be automated/randomised). Update the BED file of AMR gene locations to include all the copies.
5. Simulate the metagenome from these AMR-inserted/copied input genomes using ART. This stage also fixes the BAM file the simulator generates to include the actual contig names (note: currently hardcoded params for simulation)
6. Use samtools and the BED files to get all the reads which overlap with an AMR gene in the metagenome. A final table of labels is extracted from this. (note: check overlap requirement - could be total so missing partial AMR reads)

## Installation

There are relatively few dependencies: nextflow, ART, RGI, samtools, biopython, and pysam

[Nextflow](https://www.nextflow.io/) can be installed following the instructions on their site.

The other dependencies are just stored in a conda env yaml (`env.yaml`) and can be installed using `conda env create -f env.yaml`

There are better ways to do this but its so simple that a single env works

## Example Execution

    ./nextflow run amr_simulate.nf --input_data test/data/test_input.csv --insertion_genes test/data/insert_seqs.fna  -resume


This will generate an output folder called `simulated_metagenome`

    metagenome_unsorted.bed : BED file containing all the AMR gene locations in the metagenome
    simulated_metagenome.fna : Fasta file containing all the modified genomes in the metagenome
    simulated_metagenome_1.fq.gz : Simulated metagenomic R1
    simulated_metagenome_2.fq.gz : Simulated metagenomic R2
    simulated_metagenome_error_free.bam : Error free ground-truth of all simulated read origin locations
    AMR_metagenome_labels.tsv : TSV with read names and the AROs that read overlaps


