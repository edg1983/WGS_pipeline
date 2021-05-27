nextflow.enable.dsl=2

params.cohort_id = 'Cohort'

workflow SVTOOLS_COHORT_VCF {
    take:
        cn_files // cn.vcf files (from svtools copynumber)
        merged_file // merged.vcf, coordinates
    
    main:
        make_cohort_vcf(cn_files.collect(), merged_file)
        prune_cohort_vcf(make_cohort_vcf.out)

    emit:
        cohort_vcf = prune_cohort_vcf.out
}

process make_cohort_vcf {
    label 'singlecore'

    input: 
        file(cn_vcf)
        tuple file(merged_vcf), file(coordinates)
    
    output:
        file 'merged.sv.gt.cn.vcf'

    script:

    """
    ls *.cn.vcf > cnlist.txt
    svtools vcfpaste -m ${merged_vcf} -f cnlist.txt -q > merged.sv.gt.cn.vcf
    """
}

process prune_cohort_vcf {
    label 'singlecore'

    input:
        file(family_vcf)
    
    output:
        file '*.pruned.vcf.gz'

    script:

    """
    cat $family_vcf | svtools afreq | svtools vcftobedpe | svtools bedpesort \
    | svtools prune -s -d 100 -e \"AF\" | svtools bedpetovcf \
    | bgzip -c > ${params.cohort_id}.sv.pruned.vcf.gz
    """  
}