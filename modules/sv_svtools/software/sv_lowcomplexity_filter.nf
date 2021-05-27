params.outdir = 'SV'

workflow SV_LOWCOMPLEXITY_FILTER {
    take:
        vcf_file
        lowcomplexity_regions // file: bed file of LCRs

    main:
        convert_to_bed(vcf_file)
        lowcomplexity_filter(vcf_file,convert_to_bed.out,lowcomplexity_regions)
        
    emit:
        lowcomplexity_filter.out
}

process convert_to_bed {
    label 'singlecore'

    input:
        file(temp_vcf)
    
    output:
        file 'temp_5.noBND.forOverlap.bed'

    script:
    """
    svtools vcftobedpe -i $temp_vcf -o temp_5.bedpe
    svtools bedpetobed12 -i temp_5.bedpe -o temp_5.bed
    grep -v "BND" temp_5.bed | tail -n+2 | cut -f1-4 | tr ";" "\t" \
    | cut -f1-5 | sed 's/ID=//g' > temp_5.noBND.forOverlap.bed
    """
}

process lowcomplexity_filter {
    label 'singlecore'
    publishDir "$params.outdir", mode: 'copy'

    input:
        file(vcf_file)
        file(bed_file)
        tuple file(lowcomplexity_regions), file(lowcomplexity_index)
    
    output:
        file "*.LCR.vcf.gz"

    script:
    def outprefix = vcf_file.getName().replaceFirst(/\.vcf\.gz/, "")
    """
    bedtools coverage -wa \
    -a $bed_file \
    -b $lowcomplexity_regions \
    > overlapExcludeRegions.txt
    
    FilterSVByValue.py \
    -v $vcf_file \
    -b overlapExcludeRegions.txt \
    -c 9 -i 5 -t LowComplexity -x 0.7 -l 1000 -o temp_LCR.vcf
    sed 's/:\\./:0/g' temp_LCR.vcf | bgzip -c > ${outprefix}.LCR.vcf.gz
    """
}