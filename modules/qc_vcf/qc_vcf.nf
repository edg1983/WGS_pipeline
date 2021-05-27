nextflow.enable.dsl=2

params.outdir = 'vcf_QC'
output_reports = file("${params.outdir}")
output_reports.mkdirs()

include { MULTIQC } from './software/multiqc' addParams( outdir: params.outdir)
include { BCFTOOLS_STATS } from './software/bcftools_stats' addParams( outdir: params.outdir)

workflow QC_VCF {
    take: 
        vcf_file // tuple file(vcf_file), file(vcf_index)

    main:
        BCFTOOLS_STATS(vcf_file)
        MULTIQC(BCFTOOLS_STATS.out.collect())

    emit:
        multiqc_report = MULTIQC.out
        stats_file = BCFTOOLS_STATS.out
}


