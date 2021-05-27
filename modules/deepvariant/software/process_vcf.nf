// These processes require bcftools and VCF_toolbox
// https://github.com/edg1983/VCF_toolbox
process BCFTOOLS_NORM {
    label 'singlecore'

    input:
        tuple val(prefix), file(vcf_file), file(vcf_index)
        tuple file(genome), file(genome_fai)

    output:
        tuple val(prefix), file("${prefix}.norm.vcf.gz"), file("${prefix}.norm.vcf.gz.csi"), emit: norm_vcf
    
    script:
    """
    VCF_toolbox annotate_multiallele -v $vcf_file -o temp.vcf.gz

    bcftools norm -m- -f $genome \
    -Oz -o ${prefix}.norm.vcf.gz temp.vcf.gz

    tabix -p vcf -m 12 --csi ${prefix}.norm.vcf.gz
    """
}

process FIX_DEEPVAR_VCF {
    label 'singlecore'

    input:
        tuple val(prefix), file(vcf_file), file(vcf_index)

    output:
        tuple val(prefix), file("${prefix}.fix.vcf.gz"), file("${prefix}.fix.vcf.gz.csi"), emit: fixed_vcf

    script:
    """
    VCF_toolbox deepvar_fix \
    -v $vcf_file -o ${prefix}.fix.vcf.gz \
    --fixGQ --setMissing --fixhalfGT --computeAB
    tabix -p vcf -m 12 --csi ${prefix}.fix.vcf.gz
    """
}