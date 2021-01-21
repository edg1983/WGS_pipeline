nextflow.enable.dsl=2

/* Additional params provided in configuration file:
    params.regions = '/well/gel/HICF2/ref/geneanno/gencode/gencode.v31.annotation.exons.merged.bed.gz'
    params.mosdepth = '/well/gel/HICF2/software/bin-linux/mosdepth'
    params.somalier = '/well/gel/HICF2/software/bin-linux/somalier'
    params.fastqc = '/well/gel/HICF2/software/fastqc-0.11.2/fastqc'
    params.somalier_data = '/well/gel/HICF2/ref/somalier_data' 
*/

params.bam_folder = "BAM" /* directory containing bam files */
params.qc_folder = "QC" /* output folder for QC files */
params.ped = "input.ped"
params.ref = "/well/gel/HICF2/ref/genomes/GRCh38/bwakit/hs38DH.fa"

output_somalier = file("${params.qc_folder}/somalier")
output_files = file("${params.qc_folder}/files")
output_reports = file("${params.qc_folder}/reports")
output_reports.mkdirs()
output_files.mkdirs()
output_somalier.mkdirs()

workflow {
    log.info """\
    BAM files QC - N F   P I P E L I N E    
    ====================================
    This pipeline works for hg38 only
    BAMs file dir   (--bam_folder)  : ${params.bam_folder}
    ped file        (--ped)         : ${params.ped}
    outdir          (--qc_folder)   : ${params.qc_folder}
    applied QC                      : fastQC, somalier, flagstats, coverage
    """
    .stripIndent()
    
    bamfiles = Channel
                .fromPath("${params.bam_folder}/*.bam")
                .map { file -> tuple(file.baseName, file, file+".bai") }
    genome_ref = tuple(file(params.ref), file("${params.ref}.fai"))
    QC_bam(bamfiles, file(params.regions), file(params.somalier_data), genome_ref)
}

workflow QC_bam {
    take:
        /* a tuple containing bam files
        [sample, file.bam, file.bam.bai] */ 
        bamfiles
        regions_file
        somalier_data
        genome_ref
        
    main:
        fastqc(bamfiles)
        coverage(bamfiles, regions_file)
        map_flag_stat(bamfiles)
        somalier(bamfiles, somalier_data, genome_ref)
        sex_and_relatedness(somalier.out.collect(),file("$params.ped"))
        ancestry(somalier.out.collect(), somalier_data)
        qc_files = fastqc.out.mix(map_flag_stat.out,coverage.out,sex_and_relatedness.out,ancestry.out)
        multiqc(qc_files.collect())

    emit:
        multiqc.out
}

process fastqc {
    container '/well/gel/HICF2/software/singularity/QCbam-v1.0.sif'
    label 'singlecore'
    publishDir "$output_files", mode: 'copy'

    input:
    tuple val(sample), file(bam_file), file(bai_file)

    output:
    file("*_fastqc.zip")

    script:
    """
    fastqc --noextract $bam_file
    """
}

process coverage {
    container '/well/gel/HICF2/software/singularity/QCbam-v1.0.sif'
    label 'lowcores'
    publishDir "$output_files", mode: 'copy'

    input:
    tuple val(sample), file(bam_file), file(bai_file)
    file(regions)

    output:
    file("*.{dist.txt,bed.gz}")

    script:
    """
    export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
    export MOSDEPTH_Q1=LOW_COVERAGE  # 1..5
    export MOSDEPTH_Q2=CALLABLE_10X  # 6..10
    export MOSDEPTH_Q3=CALLABLE_20X  # 11..20
    export MOSDEPTH_Q4=CALLABLE_100X # 21..100
    export MOSDEPTH_Q5=HIGH_COV # >100

    if [ -f "$regions" ]
    then
        MOSDEPTH_PRECISION=4 mosdepth -n -t 5 --quantize 0:1:6:11:21:101: --by $regions --thresholds 1,5,10,20 ${sample} ${sample}.bam
    else
        MOSDEPTH_PRECISION=4 mosdepth -n -t 5 --quantize 0:1:6:11:21:101: ${sample} ${sample}.bam
    fi
    """
}

process map_flag_stat {
    container '/well/gel/HICF2/software/singularity/QCbam-v1.0.sif'
    label 'lowcores'
    publishDir "$output_files", mode: 'copy'

    input:
    tuple val(sample), file(bam_file), file(bai_file)

    output:
    file("*.{flagstat,mapstat}")
    
    script:
    """
    samtools flagstat -@ 5 ${sample}.bam > ${sample}.flagstat
    samtools idxstats ${sample}.bam > ${sample}.mapstat
    """
}

process somalier {
    container '/well/gel/HICF2/software/singularity/QCbam-v1.0.sif'
    label 'singlecore'
    publishDir "$output_somalier", mode: 'copy', pattern: "*.somalier"

    input:
    tuple val(sample), file(bam_file), file(bai_file)
    path(somalier_data)
    tuple file(genome_fasta), file(genome_index)

    output:
    file("*.somalier")

    script:
    """  
    somalier extract -d ./ -f $genome_fasta -s ${somalier_data}/somalier_varsets/sites.hg38.vcf.gz ${sample}.bam
    """
}

process sex_and_relatedness {
    container '/well/gel/HICF2/software/singularity/QCbam-v1.0.sif'
    label 'singlecore'
    publishDir "$output_files", mode: 'copy', pattern: '*.tsv'
    publishDir "$output_reports", mode: 'copy', pattern: '*.html'
    
    input:
    file("*.somalier")
    file("ped_file.ped")

    output:
    file 'cohort_relate.*'

    script:
    """
    somalier relate -p ped_file.ped -o cohort_relate *.somalier
    """
}

process ancestry {
    container '/well/gel/HICF2/software/singularity/QCbam-v1.0.sif'
    label 'singlecore'
    publishDir "$output_files", mode: 'copy', pattern: '*.tsv'
    publishDir "$output_reports", mode: 'copy', pattern: '*.html'
    
    input:
    file("*.somalier")
    path(somalier_data)

    output:
    file 'cohort.*'

    script:
    """
    somalier ancestry --labels ${somalier_data}/ancestry-labels-1kg.tsv -o cohort.somalier-ancestry ${somalier_data}/1kg-somalier/*.somalier ++ *.somalier
    """
}

process multiqc {
    container '/well/gel/HICF2/software/singularity/QCbam-v1.0.sif'
    label 'singlecore'
    publishDir "$output_reports", mode: 'copy'

    input:
    file('*')

    output:
    file 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}