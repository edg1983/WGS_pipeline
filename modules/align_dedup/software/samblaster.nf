// mark duplicates and output marked bam and split/disc bams
// only the log is published since resulting bam/sam usually get passed to sort

params.outdir = 'BAM'

process SAMBLASTER {
    label 'lowcores'

    publishDir "${params.outdir}", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample), path('merged.bam')

    output:
    tuple val(sample), path('dedup.bam'), path('disc.sam'), path('split.sam'), emit: bam_files
    path "${sample}.dedup.log" , emit: logfiles

    script:
    """
    samtools view -h merged.bam \
    | samblaster --excludeDups --addMateTags -d disc.sam -s split.sam 2> ${sample}.dedup.log \
    | samtools view -Sb - > dedup.bam
    """
}