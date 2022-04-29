nextflow.enable.dsl=2


// this has the disadvantage of potentially 
// inserting new AMR genes into existing AMR genes but 
// every tool will be equally disadvantages by this

process SPLIT_INSERT_GENES {
    //publishDir 'results/0.split', mode: 'copy' 
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
    publishDir 'simulated_metagenome', mode: 'copy' 
    input:
        tuple val(index), val(genome_name), path(genome), val(copy_number), path(insert_fasta)
    output:
        tuple val(genome_name), path("*_with_insertions.fna"), val(copy_number), path("*_with_insertions.fna.insertsOnly.fna")
    script:
        """
        insert_genes.py --input_genome $genome --insert_genes $insert_fasta --output ${genome_name}_with_insertions.fna
        """
}

process ANNOTATE_AMR {
    //publishDir 'results/2.annotate', mode: 'copy' , pattern: "*.bed"
    input:
        tuple val(genome_name), path(genome), val(copy_number), path(inserts)
    output:
        tuple val(genome_name), path(genome), val(copy_number), path("*_rgi_labels.bed")
    script:
        """
        rgi main --input_sequence $genome --output_file ${genome_name}_rgi_annotation --input_type contig --alignment_tool BLAST
        rgi_tsv_to_bed.py --input_rgi_tsv ${genome_name}_rgi_annotation.txt --output_bed ${genome_name}_rgi_labels.bed
        """
}

process COPY_NUMBER {
    //publishDir 'results/3.copy', mode: 'copy' 
    input:
        tuple val(genome_name), path(genome), val(copy_number), path(genome_bed)
    output:
        path("*_amplified.fna"), emit: genome
        path("*_amr_amplified.bed"), emit: bed
    script:
        """
        amplify_genome.py --input_genome $genome --input_annotation $genome_bed --copy_number $copy_number --output_prefix $genome_name
        """
}

process SIMULATE_READS {
    //publishDir 'results/4.simulated_metagenome', mode: 'copy' 
    publishDir 'simulated_metagenome', mode: 'copy', pattern: "simulated_metagenome*"
    input:
        path(metagenome_fasta)
    output:
        path("*.fq.gz"), emit: reads
        path("simulated_metagenome_error_free.bam"), emit: bam
        path(metagenome_fasta), emit: metagenome_fasta
    script:
        """
		art_illumina --in $metagenome_fasta --paired --len 250 --fcov 2 --mflen 600 --sdev 50 --seqSys MSv3 --errfree --noALN --out metagenome 
        gzip -cvf metagenome1.fq > simulated_metagenome_1.fq.gz
        rm metagenome1.fq
        gzip -cvf metagenome2.fq > simulated_metagenome_2.fq.gz
        rm metagenome2.fq

        samtools view -S -b metagenome_errFree.sam > metagenome_errFree.bam 
        samtools index metagenome_errFree.bam
        add_contig_names_to_art_bam.py --input_art_bam metagenome_errFree.bam --output_bam unsorted_bam_with_contig_names.bam

        samtools sort unsorted_bam_with_contig_names.bam > simulated_metagenome_error_free.bam
        """
}

process GET_LABELS {
    //publishDir 'results/5.read_labels', mode: 'copy'
    publishDir 'simulated_metagenome', mode: 'copy'
    input:
        path metagenome_bed
        path metagenome_bam
    output:
        path metagenome_bed //
        path "AMR_metagenome_labels.tsv"
    script:
        """
        generate_labels.py --bed $metagenome_bed --bam $metagenome_bam > AMR_metagenome_labels.tsv
        """
}


workflow {
    // inputs metadata file with copy number of each genome and path to genome
    // Fasta containing genes to randomly insert into genomes 
    int genome_index = 0
    input_data = Channel.fromPath(params.input_data)
                        .splitCsv(header:true)
                        .map{ row-> [ genome_index++, row.genome_name, file(row.fasta_path), row.copy_number ] }

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
    
    metagenome_fasta = amplified_input.genome.collectFile(name: "simulated_metagenome.fna")
    metagenome_bed = amplified_input.bed.collectFile(name: "metagenome_unsorted.bed")
    
    simulation = SIMULATE_READS(metagenome_fasta)

    labels = GET_LABELS(metagenome_bed, 
                        simulation.bam)
    
}
