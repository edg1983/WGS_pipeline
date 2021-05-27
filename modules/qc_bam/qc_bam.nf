nextflow.enable.dsl=2

/* Additional params provided in configuration file:
    params.regions = '/well/gel/HICF2/ref/geneanno/gencode/gencode.v31.annotation.exons.merged.bed.gz'
    params.somalier_data = '/well/gel/HICF2/ref/somalier_data' 
*/

params.outdir = "QC" /* output folder for QC files */

output_somalier = file("${params.outdir}/somalier")
output_fastqc = file("${params.outdir}/fastQC")
output_coverage = file("${params.outdir}/coverage")
output_mapflagstats = file("${params.outdir}/mapflag_stats")
output_reports = file("${params.outdir}/reports")
output_reports.mkdirs()
output_fastqc.mkdirs()
output_coverage.mkdirs()
output_mapflagstats.mkdirs()
output_somalier.mkdirs()

include { MULTIQC } from './software/multiqc' addParams( outdir: "$output_reports")

workflow QC_BAM {
    take:
        /* a tuple containing bam files
        [sample, file.bam, file.bam.bai] */ 
        bamfiles
        regions_file
        somalier_data
        ped_file
        genome_ref
        
    main:
        fastqc(bamfiles)
        coverage(bamfiles, regions_file)
        map_flag_stat(bamfiles)
        somalier(bamfiles, somalier_data, genome_ref)
        sex_and_relatedness(somalier.out.collect(),ped_file)
        ancestry(somalier.out.collect(), somalier_data)
        qc_files = fastqc.out.mix(map_flag_stat.out,coverage.out,sex_and_relatedness.out,ancestry.out)
        MULTIQC(qc_files.collect())

    emit:
        MULTIQC.out
}

process fastqc {
    label 'singlecore'
    publishDir "$output_fastqc", mode: 'copy'

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
    label 'lowcores'
    publishDir "$output_coverage", mode: 'copy'

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
        MOSDEPTH_PRECISION=4 mosdepth -n -t 5 --quantize 0:1:6:11:21:101: --by $regions --thresholds 1,5,10,20 ${sample} $bam_file
    else
        MOSDEPTH_PRECISION=4 mosdepth -n -t 5 --quantize 0:1:6:11:21:101: ${sample} $bam_file
    fi
    """
}

process map_flag_stat {
    label 'lowcores'
    publishDir "$output_mapflagstats", mode: 'copy'

    input:
    tuple val(sample), file(bam_file), file(bai_file)

    output:
    file("*.{flagstat,mapstat}")
    
    script:
    """
    samtools flagstat -@ 5 $bam_file > ${sample}.flagstat
    samtools idxstats $bam_file > ${sample}.mapstat
    """
}

process somalier {
    label 'singlecore'
    publishDir "$output_somalier", mode: 'copy', pattern: "*.somalier"

    input:
    tuple val(sample), file(bam_file), file(bai_file)
    path(somalier_data)
    tuple file(genome_fasta), file(genome_index)

    output:
    file("${sample}.somalier")

    script:
    """  
    somalier extract -d ./ -f $genome_fasta -s ${somalier_data}/somalier_varsets/sites.hg38.vcf.gz $bam_file
    """
}

process sex_and_relatedness {
    label 'singlecore'
    publishDir "$output_reports", mode: 'copy'
    
    input:
    file("*")
    file("ped_file.ped")

    output:
    file 'cohort_relate.*'

    script:
    """
    somalier relate -p ped_file.ped -o cohort_relate *.somalier
    """
}

process ancestry {
    label 'singlecore'
    publishDir "$output_reports", mode: 'copy'
    
    input:
    file("*")
    path(somalier_data)

    output:
    file 'cohort.*'

    script:
    """
    somalier ancestry --labels ${somalier_data}/ancestry-labels-1kg.tsv -o cohort.somalier-ancestry ${somalier_data}/1kg-somalier/*.somalier ++ *.somalier
    """
}

