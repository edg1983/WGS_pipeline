params.outdir = 'SV'
params.outprefix = 'Cohort'

process SVTOOLS_CLASSIFY_LARGESAMPLE {
    label 'singlecore'
    publishDir "$params.outdir", mode: 'copy', pattern: '*.classify.vcf.gz'

    input:
        file(pruned_vcf)
        file(sex_file) // tab-separated file: sampleID, sex (1=male,2=female)
        tuple file(mei_regions), file(mei_index)
    
    output:
        file "*.classify.vcf.gz"

    script:

    """
    zcat $pruned_vcf \
    | svtools classify \
    -g $sex_file \
    -a $mei_regions \
    -m large_sample \
    | bgzip -c > ${params.outprefix}.classify.vcf.gz
    """     
}

process SVTOOLS_CLASSIFY_SMALLSAMPLE {
    label 'singlecore'
    publishDir "$params.outdir", mode: 'copy', pattern: '*.classify.vcf.gz'

    input:
        file(pruned_vcf)
        file(sex_file)
        file(training_bedpe)
        tuple file(mei_regions), file(mei_index)
    
    output:
        file "*.classify.vcf.gz"

    script:
    def outprefix = pruned_vcf.getName().replaceFirst(/\.vcf\.gz/, "")
    """
    zcat $pruned_vcf \
    | svtools vcftobedpe  \
    | svtools varlookup -a stdin -b $training_bedpe -c HQ -d 50 \
    | svtools bedpetovcf \
    | svtools vcfsort \
    | bcftools view -i "INFO/HQ_AF > 0" \
    | bgzip -c > training.vars.vcf.gz

    zcat $pruned_vcf \
    | svtools classify \
    -g $sex_file \
    -a $mei_regions \
    -m naive_bayes \
    -t training.vars.vcf.gz \
    | bgzip -c > ${outprefix}.classify.vcf.gz
    """     
}