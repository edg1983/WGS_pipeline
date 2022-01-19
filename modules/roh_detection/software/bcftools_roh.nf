params.outdir = 'ROH'
output_dir = file(params.outdir)
output_dir.mkdirs()

process BCFTOOLS_ROH {
    label 'singlecore'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple file(vcf_file), file(vcf_index)
    tuple file(AF_file), file(AF_index)
    file genetic_maps
    val sampleID

    output:
    file("${sampleID}_ROH.txt")

    script:
    """
    bcftools roh -s $sampleID --AF-file $AF_file -m {CHROM}.genetic_map.txt -Or -o ${sampleID}_ROH.txt $vcf_file
    """
}