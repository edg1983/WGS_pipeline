nextflow.enable.dsl=2
/* Additional params provided in configuration file:
    params.exphunter_catalog = .json file of ExpHunter variant catalog
*/

params.outdir = "ExpHunter" /* output directory containing results */
output_exphunter = file("$params.outdir")
output_exphunter.mkdirs()

include { BCFTOOLS_MERGE } from './software/bcftools_merge' addParams( outdir: params.outdir, outprefix: params.cohort_id )

workflow EXPANSION_HUNTER_CALLS {
    take: 
        project /* a tuple containing [sample, sex, file.bam, file.bam.bai] */ 
        genome_ref /* a tuple with genome fasta and its .fai index */
        variant_catalog /* .json file of ExpHunter variant catalog */

    main:
        expansion_hunter(project, genome_ref, variant_catalog)
        compress_exphunter_vcf(expansion_hunter.out.vcf_files)
        if (params.samplesize > 1) {
            BCFTOOLS_MERGE(compress_exphunter_vcf.out.collect())
        }

    emit:
        json = expansion_hunter.out.json_files
        vcf = compress_exphunter_vcf.out
}

process expansion_hunter {
    label 'lowcores'
    publishDir "$params.outdir", mode: 'copy'

    input:
    tuple val(sample), val(sex), file(bam_file), file(bai_file)
    tuple file(genome_fasta), file(genome_index)
    file(variant_catalog)

    output:
    path "${sample}.json", emit: json_files
    path "${sample}.vcf", emit: vcf_files
    path "${sample}_realigned.bam", emit: bam_files

    script:
    """
    ExpansionHunter --reads $bam_file \
    --reference $genome_fasta \
    --variant-catalog $variant_catalog \
    --sex $sex \
    --threads 5 \
    --output-prefix $sample 2> ${sample}.log
    """
}

process compress_exphunter_vcf {
    label 'singlecore'

    input:
    file(vcf_file)

    output:
    tuple file("${vcf_file}.gz"), file("${vcf_file}.gz.tbi" )

    script:
    """
    bgzip $vcf_file
    tabix -p vcf ${vcf_file}.gz
    """
}