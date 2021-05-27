params.outdir = 'cnvnator_root'

process CNVNATOR {
    label 'lowcores'
    publishDir "$params.outdir", mode: 'copy', pattern: "*.root"

    input:
        tuple val(sampleID), file(bam_file), file(bai_file)
        path chrs_folder

    output:
        tuple val(sampleID), file("${sampleID}.root")

    script:
    STDCHRS = (1..22).collect { "chr$it"}.join(' ')
    STDCHRS = "$STDCHRS chrX chrY chrM"
    
    """
    cnvnator -root ${sampleID}.root -chrom $STDCHRS -tree $bam_file
    
    cnvnator -root ${sampleID}.root -d $chrs_folder -his 100 
	cnvnator -root ${sampleID}.root -d $chrs_folder -his 1000 

	cnvnator -root ${sampleID}.root -stat 100
	cnvnator -root ${sampleID}.root -stat 1000

	cnvnator -root ${sampleID}.root -partition 100
	cnvnator -root ${sampleID}.root -partition 1000
    """
}