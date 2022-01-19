nextflow.enable.dsl=2

params.outdir = 'SV'
params.cohort_id = 'Cohort'
params.samplesize = 1

output_lumpy = file("${params.outdir}/raw_calls/lumpy_VCF")
output_cnvnator = file("${params.outdir}/raw_calls/cnvnator_root")

if (params.samplesize == 1) {
    output_cn_files = file("${params.outdir}/cnVCF")
} else {
    output_cn_files = file("${params.outdir}/single_cnVCF")
}

output_lumpy.mkdirs()
output_cnvnator.mkdirs()
output_cn_files.mkdirs()

include { LUMPY_CALLS } from './software/lumpy' addParams( outdir: "$output_lumpy")
include { SVTYPER } from './software/svtyper' addParams( outdir: "$output_lumpy")
include { CNVNATOR } from './software/cnvnator' addParams( outdir: "$output_cnvnator")
include { FIX_SVTOOLS_VCF } from './software/svtools_fix_vcf' addParams( outdir: params.outdir)

if (params.samplesize == 1) {
    include { SVTOOLS_CN_SINGLE as SVTOOLS_CN } from './software/svtools_copynumber' addParams( outdir: "$output_cn_files")
} else {
    include { SVTOOLS_MERGE_CALLS } from './software/svtools_make_merged'
    include { SVTOOLS_REGENOTYPE } from './software/svtools_regenotype' 
    include { SVTOOLS_CN_COHORT as SVTOOLS_CN } from './software/svtools_copynumber' addParams( outdir: "$output_cn_files")
    include { SVTOOLS_COHORT_VCF } from './software/svtools_make_cohort_vcf' addParams( cohort_id: params.cohort_id)
}


workflow SV_CALL_SVTOOLS {
    take:
        discsplit_bam
        // tuple sampleID, file(bam_file), file(bai_file), file(split_bam), file(disc_bam)
        regular_bam
        // tuple sampleID, file(bam_file), file(bai_file)
        ref_genome
        // tuple file(genome_fasta), file(genome_fai)
        exclude_regions
        // file(bed_file) regions to exclude when performing lumpy
        chrs_folder
        // file(chrs_folder) a folder containing chrs.fa sequences (like chr1.fa ...)

    main:
        //Determine structural variants and CNV
        LUMPY_CALLS(discsplit_bam, exclude_regions)
        SVTYPER(LUMPY_CALLS.out.vcf_files, ref_genome)
        CNVNATOR(regular_bam, chrs_folder)
        
        // When single sample combine CNVnator and lumpy output
        if (params.samplesize == 1) {
            lumpy_root_ch = SVTYPER.out.join(CNVNATOR.out)
            SVTOOLS_CN(lumpy_root_ch)
            FIX_SVTOOLS_VCF(SVTOOLS_CN.out)
        } else {
            // When multiple samples merge calls from all samples using svtools pipeline
            vcf_to_merge_ch = SVTYPER.out.collect({ it[1] })
            SVTOOLS_MERGE_CALLS(vcf_to_merge_ch)
            SVTOOLS_REGENOTYPE(regular_bam, SVTOOLS_MERGE_CALLS.out.merged_file)
            lumpy_root_ch = SVTOOLS_REGENOTYPE.out.join(CNVNATOR.out)
            SVTOOLS_CN(lumpy_root_ch, SVTOOLS_MERGE_CALLS.out.merged_file)
            SVTOOLS_COHORT_VCF(SVTOOLS_CN.out, SVTOOLS_MERGE_CALLS.out.merged_file)
            FIX_SVTOOLS_VCF(SVTOOLS_COHORT_VCF.out.cohort_vcf)
        }
        
    emit:
        sv_vcf = FIX_SVTOOLS_VCF.out
}