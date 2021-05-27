nextflow.enable.dsl=2

/*
This runs deepvariant and output .vcf, .g.vcf and .html report
Default mode is for WGS and output folder is VCF
*/

params.mode = "WGS"
params.outdir = "VCF"

output_vcfs = file("${params.outdir}")
output_htmls = file("${params.outdir}/deepvariant_html")
output_htmls.mkdirs()

process DEEPVARIANT {
    label 'highcores_AVX'
    publishDir "$output_vcfs", mode: 'copy', pattern: "*.vcf.gz*"
    publishDir "$output_htmls", mode: 'copy', pattern: "*.html"

    input:
        tuple file(genome), file(genome_fai)
        tuple val(sampleID), file(bam_file), file(bai_file)

    output:
        path "*.g.vcf.gz", emit: gvcf
        tuple val("${sampleID}"), file("${sampleID}.vcf.gz"), file("${sampleID}.vcf.gz.tbi"), emit: vcf
        path "*.html", emit: html

    script:
    STDCHRS = (1..22).collect { "chr$it"}.join(' ')
    STDCHRS = "$STDCHRS chrX chrY chrM"

    """
    mkdir temp_dir && \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${params.mode} --ref="$genome" \
    --reads="$bam_file" \
    --output_vcf="${sampleID}.vcf.gz" \
    --output_gvcf="${sampleID}.g.vcf.gz" \
    --intermediate_results_dir="temp_dir" \
    --num_shards=10 \
    --regions="$STDCHRS"
    """
}