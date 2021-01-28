nextflow.enable.dsl=2
/* Additional params provided in configuration file:
    params.roh_gmaps = '/well/gel/HICF2/ref/ROH_resources/Genetic_maps'
    The above folder should contain one file per chromosome named as {chrom}.genetic_map_hg38.txt
    params.roh_AFfile = '/well/gel/HICF2/ref/ROH_resources/gnomAD_v3_AF.tab.gz'
*/

params.vcf_file = "vcf_file.vcf.gz"
params.out_folder = "ROH" /* output directory containing for ROH files */

output_roh = file("$params.out_folder")
output_roh.mkdirs()

workflow {
    log.info """\
    ROH detection - N F   P I P E L I N E    
    =====================================
    This pipeline works for hg38 only
    input VCF       (--vcf_file)    : ${params.vcf_file}
    outdir          (--out_folder)  : ${params.out_folder}
    AF file         (--roh_AFfile)  : ${params.roh_AFfile}
    """
    .stripIndent()
    
    vcf_file = tuple(file("${params.vcf_file}"), file("${params.vcf_file}.csi"))
    AF_file = tuple(file("${params.roh_AFfile}"), file("${params.roh_AFfile}.csi"))
    genetic_maps = file("${params.roh_gmaps}/*.txt")
    ROH_detection(vcf_file, AF_file, genetic_maps)
}

workflow ROH_detection {
    take: 
        vcf_file /* a tuple with vcf_file and its index */
        AF_file /* a tuple with AF tab file and its index */
        genetic_maps /* a series of files {chrom}.genetic_map_hg38.txt */

    main:
        getSampleIDs(vcf_file)
        sampleIDs = getSampleIDs.out.splitText() { it.trim() }
        bcftools_ROH(vcf_file,AF_file,genetic_maps,sampleIDs)

    emit:
        bcftools_ROH.out
}

process getSampleIDs {
    executor 'local'

    input:
    tuple file(vcf_file), file(vcf_index)

    output:
    file("samples.txt")

    script:
    """
    zgrep -m1 "#CHROM" $vcf_file | cut -f10- | tr "\\t" "\\n" > samples.txt
    """
}

process bcftools_ROH {
    label 'singlecore'
    publishDir "$output_roh", mode: 'copy'

    input:
    tuple file(vcf_file), file(vcf_index)
    tuple file(AF_file), file(AF_index)
    file genetic_maps
    val sampleID

    output:
    file("${sampleID}_ROH.txt")

    script:
    """
    bcftools roh -s $sampleID --AF-file $AF_file -m {CHROM}.genetic_map_hg38.txt -Or -o ${sampleID}_ROH.txt $vcf_file
    """
}