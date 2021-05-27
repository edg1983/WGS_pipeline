params.outdir = 'reports'

process BCFTOOLS_STATS {
    label 'singlecore'
    
    publishDir "$params.outdir", mode: 'copy', pattern: "*.stats"

    input:
        tuple file(vcf_file), file(vcf_index)

    output:
        path "*.stats", emit: stats

    script:
    """
    bcftools stats -s- $vcf_file > ${vcf_file}.stats
    bcftools index --stats $vcf_file > ${vcf_file}_varsPerChrom.stats
    """
}