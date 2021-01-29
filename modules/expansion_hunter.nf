nextflow.enable.dsl=2
/* Additional params provided in configuration file:
    params.exphunter_catalog = .json file of ExpHunter variant catalog
*/

params.input = "input.txt" /* tab-separated file: sampleID, sex(male/female), bamfile */
params.out_folder = "ExpHunter" /* output directory containing results */
params.ref="/well/gel/HICF2/ref/genomes/GRCh38/bwakit/hs38DH.fa"

output_exphunter = file("$params.out_folder")
output_exphunter.mkdirs()

workflow {
    log.info """\
    ExpnasionHunter - N F   P I P E L I N E    
    =======================================
    This pipeline works for hg38 only
    input file  (--input)               : ${params.input}
    outdir      (--out_folder)          : ${params.out_folder}
    var catalog (--exphunter_catalog)   : ${params.exphunter_catalog}
    """
    .stripIndent()

    designFile = file(params.input)
    project = Channel
                .from(designFile)
                .splitCsv(header: false, sep: '\t')
                .map{row -> return tuple(
                    row[0], row[1],
                    file(row[2]), 
                    file(row[2]+".bai"))}
    genome_ref = tuple(file(params.ref), file("${params.ref}.fai"))
    variant_catalog = file(params.exphunter_catalog)
    ExpansionHunter(project, genome_ref, variant_catalog)
}

workflow ExpansionHunter {
    take: 
        project /* a tuple containing [sample, sex, file.bam, file.bam.bai] */ 
        genome_ref /* a tuple with genome fasta and its .fai index */
        variant_catalog /* .json file of ExpHunter variant catalog */

    main:
        exphunter_call(project, genome_ref, variant_catalog)
        compress_vcf(exphunter_call.out.vcf_files)
        merge_vcf(compress_vcf.out.collect())

    emit:
        exphunter_call.out.json_files
}

process exphunter_call {
    label 'singlecore'
    publishDir "$output_exphunter", mode: 'copy'

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
    --output-prefix $sample 2> ${sample}.log
    """
}

process compress_vcf {
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

process merge_vcf {
    label 'singlecore'
    publishDir "$output_exphunter", mode: 'copy'

    input:
    file(vcf_files)

    output:
    file("ExpHunter_merged.vcf.gz")

    script:
    def input_vcfs = vcf_files.collect{ "$it" }.findAll{ !it.contains('tbi') }.join(" ")
    """
    bcftools merge -m all -Oz -o ExpHunter_merged.vcf.gz $input_vcfs
    """
}