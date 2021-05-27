nextflow.enable.dsl=2

workflow SVTOOLS_MERGE_CALLS {
    take:
        vcfgt_files 
        // lumpy vcf files genotyped
    
    main:
        concat_lumpy_calls(vcfgt_files)
        merge_lumpy_calls(concat_lumpy_calls.out)

    emit:
        merged_file = merge_lumpy_calls.out
}

process concat_lumpy_calls {
    label 'singlecore'

    input:
        file(lumpy_vcf)
    
    output:
        file('sorted.vcf')

    script:
    def lumpy_files = lumpy_vcf.collect{ "$it" }.join(" ")
    """
    svtools lsort $lumpy_files > sorted.vcf
    """
}

process merge_lumpy_calls {
    label 'singlecore'

    input:
        file(sorted_vcf)
    
    output:
        tuple file('merged.vcf'), file('coordinates')

    script:
    """
    svtools lmerge -i $sorted_vcf -f 20 | grep -v "KN\\|JTFH\\|HLA" > merged.vcf
    create_coordinates -i merged.vcf -o coordinates
    """
}