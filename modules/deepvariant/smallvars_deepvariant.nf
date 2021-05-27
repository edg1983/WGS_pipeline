nextflow.enable.dsl=2

params.outdir = "VCF"
output_rawvcf = file("${params.outdir}/raw_vcf")
if (params.samplesize == 1) {
    output_deepvar = file("${params.outdir}/raw_vcf")
} else {
    output_deepvar = file("${params.outdir}/single_vcf")
}
output_rawvcf.mkdirs()
output_deepvar.mkdirs()

include { DEEPVARIANT } from "./software/deepvariant" addParams(outdir: "$output_deepvar", mode: params.mode)
include { BCFTOOLS_NORM; FIX_DEEPVAR_VCF } from './software/process_vcf'
include { BCFTOOLS_FILTER } from './software/bcftools_filter' addParams( outdir: params.outdir )

if (params.samplesize > 1) {
    include { GLNEXUS_MERGE_GVCF } from './software/glnexus' addParams( outdir: "$output_rawvcf" )
}

workflow SMALLVARS_DEEPVAR {
    take: 
        ref_genome // file(genome.fa), file(genome.fa.fai)
        regular_bam // val(sampleID), file(bam), file(bam.bai)

    main:
        // Set filter string according to sample size
        // string to be passed to bcftools filter -e
        def filter_string = ''
        if (params.samplesize == 1) {
            filter_string = params.filters.smallvars.single
        }
        if (params.samplesize <= 10) {
            filter_string = params.filters.smallvars.small
        }
        if (params.samplesize > 10) {
            filter_string = params.filters.smallvars.large
        }

        filter_exp = Channel.value(filter_string)

        log.info """\
            --- SMALL VARIANTS ---

            Based on sample size, the following filter string is applied for smallvars
            $filter_string
            =============================================
            """
            .stripIndent()
            
        DEEPVARIANT(ref_genome, regular_bam)
        if (params.samplesize > 1) {
            GLNEXUS_MERGE_GVCF(DEEPVARIANT.out.gvcf.collect())
            BCFTOOLS_NORM(GLNEXUS_MERGE_GVCF.out, ref_genome)
        } else {
            BCFTOOLS_NORM(DEEPVARIANT.out.vcf, ref_genome)
        }
        FIX_DEEPVAR_VCF(BCFTOOLS_NORM.out.norm_vcf)
        BCFTOOLS_FILTER(FIX_DEEPVAR_VCF.out.fixed_vcf, filter_exp)

    emit:
        pass_vcf = BCFTOOLS_FILTER.out.pass_vcf
        processed_vcf = BCFTOOLS_FILTER.out.processed_vcf
}