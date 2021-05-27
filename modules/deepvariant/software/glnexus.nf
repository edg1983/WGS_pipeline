nextflow.enable.dsl=2

/*
Use GLNexus to merge multiple g.vcf generate by deepvariant into a single cohort
Output file will be cohort_id.vcf.gz
*/

params.cohort_id = 'Cohort'
params.mode = "WGS"

params.outdir = "VCF"
output_vcf = file("${params.outdir}")
output_vcf.mkdirs()

process GLNEXUS_MERGE_GVCF {
    label 'highcores'
    publishDir "$output_vcf", mode: 'copy', pattern: "${params.cohort_id}.vcf.gz*"

    input:
        file(gvcf_files)

    output:
        tuple val("${params.cohort_id}"), file("${params.cohort_id}.vcf.gz"), file("${params.cohort_id}.vcf.gz.csi")

    script:
    def input_vcfs = gvcf_files.collect{ "$it" }.join(" ")

    """
    glnexus_cli -m 64 -t 10 \
    -c DeepVariant${params.mode} \
    --dir GLnexusDB \
    $input_vcfs \
    | bcftools view - \
    | bgzip -@ 10 -c > ${params.cohort_id}.vcf.gz

    tabix -p vcf -m 12 --csi ${params.cohort_id}.vcf.gz
    """
}