nextflow.enable.dsl=2

/* Additional params provided in configuration file:
    bwakit = bwakit folder
    samblaster = samblaster executable
*/

params.input = 'input_file.txt'
params.bam_folder = "BAM"
params.ref="/well/gel/HICF2/ref/genomes/GRCh38/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna"
/* note that bwa index file and .alt file must be available together with the
ref genome fasta in the same location */

output_bam = file(params.bam_folder)
output_bam.mkdirs()

workflow {
    log.info """\
    Alignment - N F   P I P E L I N E    
    =================================
    This pipeline works for hg38 only
    input file  (--input)       : ${params.input}
    genome      (--ref)         : ${params.ref}
    outdir      (--bam_folder)  : ${params.bam_folder}
    """
    .stripIndent()

    designFile = file(params.input)
    fastqfiles = Channel
                   .from(designFile)
                   .splitCsv(header: false, sep: '\t')
                   .map{row -> return tuple(row[0], file(row[1]), file(row[2]))}
    genome_file = file("${params.ref}")
    genome_base = genome_file.getBaseName()
    genome_fasta = genome_file.getName()
    genome_folder = genome_file.getParent()
    genome_data = tuple("$genome_fasta", file("${genome_folder}/${genome_base}.*"))
    align_dedup(fastqfiles, genome_data)
}

workflow align_dedup {
    take:
        /* a tuple with [sampleID, R1.fastq, R2.fastq] */ 
        fastqfiles
        genome_data
        
    main:
        singlePartAlign(fastqfiles, genome_data)
        mergeBams(singlePartAlign.out.part_bams.groupTuple())
        deduplicate(mergeBams.out.merged_bams)
        sortBams(deduplicate.out.dedup_files)
    
    emit:
        bamfiles = sortBams.out.bam_files
        allbams = sortBams.out.all_bams
        /* tuple val(sample), file(bam), file(bam.bai), file(split.bam), file(disc.bam)}") */
}

process singlePartAlign {
    container '/well/gel/HICF2/software/singularity/AlignDedup-v1.0.sif'
    label 'highcores'
    
    input:
    tuple val(sample), path(reads_R1), path(reads_R2)
    tuple val(genome_fasta), file(genome_files)

    output:
    tuple val(sample), file("${prefix}.part.bam"), emit: part_bams    

    script:
    prefix = reads_R1.getSimpleName()
    """
    bwa mem -t10 -R\"@RG\\tID:${sample}\\tSM:${sample}\\tPL:Illumina\" ${genome_fasta} $reads_R1 $reads_R2 2> ${prefix}.log.bwamem \
    | k8 /opt/bwa.kit/bwa-postalt.js -p ${prefix}.hla ${genome_fasta}.alt \
    | samtools view -1 - > ${prefix}.part.bam
    """
}

process mergeBams {
    container '/well/gel/HICF2/software/singularity/AlignDedup-v1.0.sif'
    label 'highcores'

    input:
    /*  A tuple of files to be merged per sample
        like: (sample_id, [part1.bam, part2.bam]) */
    tuple val(sample), file(myfiles)

    output:
    tuple val(sample), file('merged.bam'), emit: merged_bams

    script:
    def input_args = myfiles.collect{ "$it" }.join(" ")
    
    """
    samtools merge -n -f -@ 10 merged.bam $input_args
    """
}

process deduplicate {
    container '/well/gel/HICF2/software/singularity/AlignDedup-v1.0.sif'
    label 'lowcores'

    publishDir "${params.bam_folder}", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample), path('merged.bam')

    output:
    tuple val(sample), path('dedup.bam'), path('split.sam'), path('disc.sam'), emit: dedup_files
    path "${sample}.dedup.log" , emit: logfiles

    script:
    """
    samtools view -h merged.bam \
    | samblaster --excludeDups --addMateTags -d disc.sam -s split.sam 2> ${sample}.dedup.log \
    | samtools view -Sb - > dedup.bam
    """
}

process sortBams {
    container '/well/gel/HICF2/software/singularity/AlignDedup-v1.0.sif'
    label 'highcores'
    publishDir "${params.bam_folder}", mode: 'copy', pattern: "*.bam*"

    input:
    tuple val(sample), path('dedup.bam'), path('split.sam'), path('disc.sam')
    
    output:
    tuple val(sample), file("${sample}.bam"), file("${sample}.bam.bai"), file("${sample}.split.bam"), file("${sample}.disc.bam"), emit: all_bams
    tuple val(sample), file("${sample}.bam"), file("${sample}.bam.bai"), emit: bam_files

    script:
    """
    samtools sort -@ 10 -m12G -T disc.tmp -o ${sample}.disc.bam disc.sam
    samtools sort -@ 10 -m12G -T split.tmp -o ${sample}.split.bam split.sam
    samtools sort -@ 10 -m12G -T dedup.tmp -o ${sample}.bam dedup.bam
    samtools index ${sample}.bam
    """
}