nextflow.enable.dsl=2

// this has the disadvantage of potentially 
// inserting new AMR genes into existing AMR genes but 
// every tool will be equally disadvantages by this
process SPLIT_INSERT_GENES {
    input:
        path(insertion_fasta)
        val(number_of_chunks)
    output:
        path "insert_chunk_*.fna", emit: insert_chunks
    script:
        """
        split_insertion_genes.py --input_fasta $insertion_fasta --number_of_chunks $number_of_chunks
        """
}

process INSERT_GENES {
    input:
        tuple val(index), path(genome), val(copy_number), path(insert_fasta)
    output:
        tuple path("genome_with_insertions.fna"), val(copy_number)
    script:
        """
        insert_genes.py --input_genome $genome --insert_genes $insert_fasta --output genome_with_insertions.fna
        """
}

process ANNOTATE_AMR {
    input:
        tuple path(genome), val(copy_number) 
    output:
        tuple path(genome), val(copy_number), path("rgi_labels.bed")
    script:
        """
        rgi main --input_sequence $genome --output_file rgi_annotation --input_type contig --alignment_tool BLAST
        rgi_tsv_to_bed.py --input_rgi_tsv rgi_annotation.txt --output_bed rgi_labels.bed
        """
}

process COPY_NUMBER {
    publishDir 'simulated_metagenome', mode: 'copy' 
    input:
        tuple path(genome), val(copy_number), path(genome_bed)
    output:
        path("genome_copied.fna"), emit: genome
        path("rgi_labels_copied.bed"), emit: bed
    script:
        """
        amplify_genome.py --input_genome $genome --input_annotation $genome_bed --copy_number $copy_number 
        """
}

process SIMULATE_READS {
    publishDir 'simulated_metagenome', mode: 'copy' 
    input:
        path(metagenome_fasta)
    output:
        path("*.fq.gz"), emit: reads
        path("metagenome_errFree.bam"), emit: bam
        path(metagenome_fasta), emit: metagenome_fasta
    script:
        """
        art_illumina --in $metagenome_fasta --paired --len 250 --fcov 2 --mflen 600 --sdev 50 --seqSys MSv3 --errfree --noALN --out metagenome 

        gzip metagenome1.fq metagenome2.fq
        samtools view -S -b metagenome_errFree.sam > metagenome_errFree.bam 
        """
}


process GET_LABELS {
// extract reads which overlap with amr annotation bed file
// samtools view -L all_AMR.bed errfree.sam 
    publishDir 'simulated_metagenome', mode: 'copy'
    input:
        path metagenome_bed, 
        path metagenome_bam
    output:
        path metagenome_bed,
        path "reads_with_amr.tsv"
    script:
        """
        echo 'test' > reads_with_amr.tsv 
        """
}



workflow {

    // inputs metadata file with copy number of each genome and path to genome
    // Fasta containing genes to randomly insert into genomes 
    params.random_seed = 42
    int genome_index = 0
    input_data = Channel.fromPath(params.input_data)
                        .splitCsv(header:true)
                        .map{ row-> [ genome_index++, file(row.fasta_path), row.copy_number ] }

    // split amr fasta into the same number of genomes 
    number_of_genomes = input_data.count()
    insert_genes = Channel.fromPath(params.insertion_genes)

    int insert_index = 0
    gene_chunks = SPLIT_INSERT_GENES(insert_genes, number_of_genomes)
    gene_chunks = gene_chunks.flatten().map{ [insert_index++, it] }

    input_data = input_data.join( gene_chunks )
    input_data_with_insertions = INSERT_GENES(input_data)

    annotated_genomes = ANNOTATE_AMR(input_data_with_insertions)

    amplified_input = COPY_NUMBER(annotated_genomes)
    
    metagenome_fasta = amplified_input.genome.collectFile(name: "metagenome.fna", 
                                                          newLine: true)
    metagenome_bed = amplified_input.bed.collectFile(name: "metagenome.bed",
                                                      newLine: true)
    
    simulation = SIMULATE_READS(metagenome_fasta)

    labels = GET_LABELS(metagenome_bed, simulation.bam)
    
}
