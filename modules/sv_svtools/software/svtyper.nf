params.outdir = "svtyper_VCF"

process SVTYPER {
    label 'lowcores'
    publishDir "$params.outdir", mode: 'copy', pattern: "*.gt.vcf.gz"

    input:
        tuple val(sampleID), file(lumpy_vcf), file(bam_file), file(bai_file)
        tuple file(genome_fasta), file(genome_fai)

    output:
        tuple val(sampleID), path("*.gt.vcf.gz")

    script:
    """
    svtyper-sso --cores 5 --batch_size 1000 --max_reads 2500 -B $bam_file -T $genome_fasta -i $lumpy_vcf | bgzip -c > ${sampleID}_lumpy.gt.vcf.gz
    """   
}