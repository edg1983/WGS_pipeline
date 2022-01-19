nextflow.enable.dsl=2

params.outdir = "SV"
output_dir = file(params.outdir)
output_dir.mkdirs()

if (params.samplesize >= 30) {
    include { SVTOOLS_CLASSIFY_LARGESAMPLE as SVTOOLS_CLASSIFY } from './software/svtools_classify' addParams( outdir: params.outdir, outprefix: params.cohort_id )
} else if (params.samplesize >= 10) {
    include { SVTOOLS_CLASSIFY_SMALLSAMPLE as SVTOOLS_CLASSIFY } from './software/svtools_classify' addParams( outdir: params.outdir )
}

if (params.samplesize == 1) {
    include { QUALITY_FILTER_SINGLESAMPLE as SV_QUALITY_FILTER } from './software/sv_quality_filters'
} else {
    include { QUALITY_FILTER_LARGESAMPLE as SV_QUALITY_FILTER } from './software/sv_quality_filters'
}

include { SV_LOWCOMPLEXITY_FILTER } from './software/sv_lowcomplexity_filter' addParams( outdir: params.outdir )

workflow SV_FILTER_SMALLSAMPLE {
    take:
        vcf_file
        sex_definitions // a list of values 'sampleID\tsex' where sex=1(M) / 2(F)
        lowcomplexity_regions
        mei_regions
        training_bedpe

    main:
        sex_file=sex_definitions.collectFile(name: 'sex.txt', newLine: true)
        
        if (params.samplesize >= 10) {
            SVTOOLS_CLASSIFY(vcf_file, sex_file, training_bedpe, mei_regions)
            SV_QUALITY_FILTER(SVTOOLS_CLASSIFY.out)
        } else {
            SV_QUALITY_FILTER(vcf_file)
        }
        SV_LOWCOMPLEXITY_FILTER(SV_QUALITY_FILTER.out, lowcomplexity_regions)
      
    emit:
        SV_LOWCOMPLEXITY_FILTER.out
}

workflow SV_FILTER_LARGESAMPLE {
    take:
        vcf_file
        sex_definitions // a list of values 'sampleID\tsex' where sex=1(M) / 2(F)
        lowcomplexity_regions
        mei_regions

    main:
        sex_file=sex_definitions.collectFile(name: 'sex.txt', newLine: true)

        SVTOOLS_CLASSIFY(vcf_file, sex_file, mei_regions)
        SV_QUALITY_FILTER(SVTOOLS_CLASSIFY.out)
        SV_LOWCOMPLEXITY_FILTER(SV_QUALITY_FILTER.out, lowcomplexity_regions)
      
    emit:
        SV_LOWCOMPLEXITY_FILTER.out
}