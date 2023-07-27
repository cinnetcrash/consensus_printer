#!/usr/bin/env nextflow

params.reads = "$baseDir/*.fastq" // Input reads
params.ref = "ref_gen/reference.fasta" // Reference genome
params.outdir = "consensus_output_batch_1" // Output directory

Channel
    .fromPath(params.reads)
    .map { file -> [file.baseName, file] }
    .set { ch_fastq }

process TrimAdapters {
    container 'https://hub.docker.com/r/genomicpariscentre/porechop' // replace with the actual image name

    input:
    tuple val(sample_name), path(fastq) from ch_fastq
    
    output:
    tuple val(sample_name), path("${sample_name}_trimmed.fastq") into ch_trimmed

    """
    porechop -i $fastq -o ${sample_name}_trimmed.fastq
    """
}

process AlignSort {
    container 'https://github.com/Niema-Docker/minimap2_samtools' // replace with the actual image name

    input:
    tuple val(sample_name), path(fastq) from ch_trimmed
    
    output:
    tuple val(sample_name), path("${sample_name}_sorted.bam"), path("${sample_name}_sorted.bam.bai") into ch_sorted
    
    """
    minimap2 -ax map-ont ${params.ref} $fastq | samtools sort -o ${sample_name}_sorted.bam
    samtools index ${sample_name}_sorted.bam
    """
}

process CallVariants {
    container 'https://hub.docker.com/r/biocontainers/bcftools' // replace with the actual image name

    input:
    tuple val(sample_name), path(bam), path(bai) from ch_sorted
    
    output:
    tuple val(sample_name), path("${sample_name}_variants.vcf.gz"), path("${sample_name}_variants.vcf.gz.csi") into ch_variants

    """
    bcftools mpileup -Ou -f ${params.ref} $bam | bcftools call -mv -Oz -o ${sample_name}_variants.vcf.gz
    bcftools index ${sample_name}_variants.vcf.gz
    """
}

process AlignmentStats {
    container 'https://hub.docker.com/r/biocontainers/samtools/' // replace with the actual image name

    input:
    tuple val(sample_name), path(bam), path(bai) from ch_sorted

    output:
    path("${sample_name}_alignment_stats.txt") into ch_stats

    """
    samtools flagstat $bam > ${sample_name}_alignment_stats.txt
    """
}

process GenerateConsensus {
    container 'https://hub.docker.com/r/biocontainers/samtools/' // replace with the actual image name

    input:
    tuple val(sample_name), path(bam), path(bai) from ch_sorted

    output:
    path("${sample_name}_consensus.fasta") into ch_consensus

    """
    samtools mpileup -aa -A -d 0 -Q 20 $bam | ivar consensus -p ${sample_name}_consensus.fasta -t 0.01 -q 8
    """
}
