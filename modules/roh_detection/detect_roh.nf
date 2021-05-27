nextflow.enable.dsl=2
/* Additional params provided in configuration file:
    params.roh_gmaps = '/well/gel/HICF2/ref/ROH_resources/Genetic_maps'
    The above folder should contain one file per chromosome named as {chrom}.genetic_map_hg38.txt
    params.roh_AFfile = '/well/gel/HICF2/ref/ROH_resources/gnomAD_v3_AF.tab.gz'
*/

params.outdir = "ROH" /* output directory containing for ROH files */

output_roh = file("$params.outdir")
output_roh.mkdirs()

include { BCFTOOLS_ROH } from './software/bcftools_roh' addParams( outdir: params.outdir)
 
workflow ROH_DETECTION {
    take: 
        vcf_file /* a tuple with vcf_file and its index */
        AF_file /* a tuple with AF tab file and its index */
        genetic_maps /* a series of files {chrom}.genetic_map_hg38.txt */

    main:
        getSampleIDs(vcf_file)
        sampleIDs = getSampleIDs.out.splitText() { it.trim() }
        BCFTOOLS_ROH(vcf_file,AF_file,genetic_maps,sampleIDs)

    emit:
        BCFTOOLS_ROH.out
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