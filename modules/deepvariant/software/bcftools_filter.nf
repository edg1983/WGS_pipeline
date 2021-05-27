params.outdir = 'processed_vcfs'

process BCFTOOLS_FILTER {
    label 'singlecore'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(prefix), file(vcf_file), file(vcf_index)
        val(filter_string) 
        // a string like 'QUAL < 20 || N_PASS(GQ >= 20) == 0'

    output:
        tuple file("${prefix}.PASS.vcf.gz"), file("${prefix}.PASS.vcf.gz.csi"), emit: pass_vcf
        tuple file("${prefix}.processed.vcf.gz"), file("${prefix}.processed.vcf.gz.csi"), emit: processed_vcf

    script:
    """
    bcftools filter -m + -s LowQual \
    -e '$filter_string' \
    -Oz -o ${prefix}.processed.vcf.gz \
    $vcf_file
    tabix -p vcf -m 12 --csi ${prefix}.processed.vcf.gz

    zgrep "#\\|PASS" ${prefix}.processed.vcf.gz \
    | bgzip -c > ${prefix}.PASS.vcf.gz
    tabix -p vcf -m 12 --csi ${prefix}.PASS.vcf.gz
    """
}