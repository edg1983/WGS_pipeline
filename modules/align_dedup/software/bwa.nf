// Run bwa alignment
// output is sample.bam or read_R1.part.bam if params.multipart_fasta is set

process BWA {
    label 'highcores'
    
    input:
    tuple val(sample), path(reads_R1), path(reads_R2)
    tuple val(genome_fasta), file(genome_files)

    output:
    tuple val(sample), file("*.part.bam"), emit: bam_file   

    script:
    def prefix = reads_R1.getSimpleName()
    """
    bwa mem -t10 -R\"@RG\\tID:${sample}\\tSM:${sample}\\tPL:Illumina\" ${genome_fasta} $reads_R1 $reads_R2 2> ${prefix}.log.bwamem \
    | k8 /opt/bwa.kit/bwa-postalt.js -p ${prefix}.hla ${genome_fasta}.alt \
    | samtools view -1 - > ${prefix}.part.bam
    """
}