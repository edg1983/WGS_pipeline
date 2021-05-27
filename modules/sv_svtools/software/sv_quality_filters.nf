process QUALITY_FILTER_LARGESAMPLE {
    label 'singlecore'

    input:
        file(pruned_vcf)
    
    output:
        file "*.filters.vcf.gz" 

    script:
    def outprefix = pruned_vcf.getName().replaceFirst(/\.vcf\.gz/, "").replaceFirst(/\.raw/, "")
    """
    zcat $pruned_vcf > temp_0.vcf
    
    bcftools filter -m + -s LowConfSmallDEL \
    -e 'INFO/SVTYPE == "DEL" && INFO/SVLEN >= -1000 && N_PASS(AS > 0) == 0' \
    temp_0.vcf > temp_1.vcf

    bcftools filter -m + -s INV_LowMSQ \
    -e 'INFO/SVTYPE == "INV" && INFO/MSQ < 150' \
    temp_1.vcf > temp_2.vcf

    bcftools filter -m + -s INV_LowLumpyEvidence \
    -e 'INFO/SVTYPE == "INV" && ((INFO/SR / INFO/SU) < 0.1 || (INFO/PE / INFO/SU) < 0.1)' \
    temp_2.vcf > temp_3.vcf

    bcftools filter -m + -s BND_LowMSQ \
    -e 'INFO/SVTYPE == "BND" && INFO/MSQ < 250' \
    temp_3.vcf > temp_4.vcf

    bcftools filter -m + -s DELDUP_LowMSQ \
    -e '(INFO/SVTYPE == "DEL" || INFO/SVTYPE == "DUP") && INFO/MSQ < 40' \
    temp_4.vcf > temp_4b.vcf

    bcftools filter -m + -s VerySmallSV \
    -e '(INFO/SVTYPE == "DEL" && INFO/SVLEN > -50) || (INFO/SVTYPE == "DUP" && INFO/SVLEN < 50)' \
    temp_4b.vcf > ${outprefix}.filters.vcf

    bgzip ${outprefix}.filters.vcf
    """    
}

process QUALITY_FILTER_SINGLESAMPLE {
    label 'singlecore'

    input:
        file(pruned_vcf)
    
    output:
        file "${outprefix}.filters.vcf.gz" 

    script:
    def outprefix = pruned_vcf.getName().replaceFirst(/\.vcf\.gz/, "")
    """
    zcat $pruned_vcf > temp_0.vcf
    
    bcftools filter -m + -s LowConfSmallDEL \
    -e 'INFO/SVTYPE == "DEL" && INFO/SVLEN >= -1000 && N_PASS(AS > 0) == 0' \
    temp_0.vcf > temp_1.vcf

    bcftools filter -m + -s INV_LowSQ \
    -e 'INFO/SVTYPE == "INV" && N_PASS(SQ >= 30) == 0' \
    temp_1.vcf > temp_2.vcf

    bcftools filter -m + -s INV_LowLumpyEvidence \
    -e 'INFO/SVTYPE == "INV" && ((INFO/SR / INFO/SU) < 0.1 || (INFO/PE / INFO/SU) < 0.1)' \
    temp_2.vcf > temp_3.vcf

    bcftools filter -m + -s BND_LowSQ \
    -e 'INFO/SVTYPE == "BND" && N_PASS(SQ >= 50) == 0' \
    temp_3.vcf > temp_4.vcf

    bcftools filter -m + -s DELDUP_LowSQ \
    -e '(INFO/SVTYPE == "DEL" || INFO/SVTYPE == "DUP") && N_PASS(SQ >= 10) == 0' \
    temp_4.vcf > temp_4b.vcf

    bcftools filter -m + -s VerySmallSV \
    -e '(INFO/SVTYPE == "DEL" && INFO/SVLEN > -50) || (INFO/SVTYPE == "DUP" && INFO/SVLEN < 50)' \
    temp_4b.vcf > ${outprefix}.filters.vcf

    bgzip ${outprefix}.filters.vcf
    """    
}