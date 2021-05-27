// This is used in svtools cohort calling to regenotype single samples

process SVTOOLS_REGENOTYPE {
    label 'lowcores'

    input:
        tuple val(sampleID), file(bam_file), file(bai_file)
        tuple file(merged_vcf), file(coordinates)

    output:
        tuple val(sampleID), file("${params.cohort_id}_${sampleID}.gt.vcf")
        
    script:
    """
    (grep "^#" $merged_vcf && grep -v "^#" $merged_vcf | awk '{OFS="\\t"};{ \$6="."; print }') \
    | svtools genotype -B $bam_file -l ${sampleID}.bam.json \
    | sed 's/PR...=[0-9\\.e,-]*\\(;\\)\\{0,1\\}\\(\\t\\)\\{0,1\\}/\\2/g' - \
    > ${params.cohort_id}_${sampleID}.gt.vcf
    """
}