// sort bam or sam file by coordinate
// accepts a tuple sampleID, bam_file, disc_bam, split_bam

params.outdir = 'BAM'

process SORT_BAM {
    label 'highcores'
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.bam*"

    input:
    tuple val(sample), path(bamfile), path(discbam), path(splitbam)
    
    output:
    tuple val(sample), file("${sample}.sorted.bam"), file("${sample}.disc.bam"), file("${sample}.split.bam"), emit: all_bams
    tuple val(sample), file("${sample}.sorted.bam"), file("${sample}.sorted.bam.bai"), emit: main_bam
    path("*.sorted.bam.bai"), emit: idx_files

    script:
    """
    samtools sort -@ 10 -m12G -T bam_tmp -o ${sample}.sorted.bam $bamfile
    samtools sort -@ 10 -m12G -T disc_tmp -o ${sample}.disc.bam $discbam
    samtools sort -@ 10 -m12G -T split_tmp -o ${sample}.split.bam $splitbam
    samtools index ${sample}.sorted.bam
    samtools index ${sample}.disc.bam
    samtools index ${sample}.split.bam
    """
}