// use svtools to add CN annotation for SVs based on CNVnator root files

params.outdir = 'single_cnVCF'
output_dir = file(params.outdir)
output_dir.mkdirs()

process SVTOOLS_CN_COHORT {
    label 'lowcores'
    publishDir "$params.outdir", mode: 'copy', pattern: '*.cn.vcf'

    input:
        tuple val(sampleID), file(lumpy_gt_file), file(root_file)
        tuple file(merged_vcf), file(coordinates)

    output:
        file("${params.cohort_id}_${sampleID}.cn.vcf")

    script:
    """
    svtools copynumber --cnvnator cnvnator -w 100 \
    -s $sampleID \
    -r $root_file \
    -c $coordinates \
    -i $lumpy_gt_file \
    > ${params.cohort_id}_${sampleID}.cn.vcf
    """
}

process SVTOOLS_CN_SINGLE {
    label 'lowcores'
    publishDir "$params.outdir", mode: 'copy', pattern: '*.cn.vcf'

    input:
        tuple val(sampleID), file(lumpy_gt_file), file(root_file)
        // lumpy_gt_file is compressed .vcf.gz

    output:
        file("${sampleID}.cn.vcf.gz")

    script:
    """
    zgrep -v "KN\\|JTFH\\|HLA" $lumpy_gt_file > clean.vcf
    create_coordinates -i clean.vcf -o coordinates
    
    svtools copynumber --cnvnator cnvnator -w 100 \
    -s $sampleID \
    -r $root_file \
    -c coordinates \
    -i $lumpy_gt_file \
    > ${sampleID}.cn.vcf
    bgzip ${sampleID}.cn.vcf
    """
}