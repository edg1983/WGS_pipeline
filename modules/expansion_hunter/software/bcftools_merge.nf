params.outdir = 'VCF'

process BCFTOOLS_MERGE {
    label 'singlecore'
    publishDir "$params.outdir", mode: 'copy'

    input:
    file(vcf_files)

    output:
    file("Merged_vcf.vcf.gz")

    script:
    def input_vcfs = vcf_files.collect{ "$it" }.findAll{ !it.contains('tbi') }.join(" ")
    """
    bcftools merge -m all -Oz -o Merged_vcf.vcf.gz $input_vcfs
    """
}