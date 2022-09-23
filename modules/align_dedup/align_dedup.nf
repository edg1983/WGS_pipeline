nextflow.enable.dsl=2

params.outdir = 'BAM'
output_bam = file(params.outdir)
output_bam.mkdirs()

if (params.build == "GRCh37") {
    include { BWA_GRCh37 as BWA   } from './software/bwa'
} else if (params.build == "GRCh38") {
    include { BWA_GRCh38 as BWA   } from './software/bwa'
}
include { MERGE_BAMS    } from './software/samtools_merge'
include { SAMBLASTER    } from './software/samblaster' addParams( outdir: params.outdir)
include { SORT_BAM      } from './software/samtools_sort'   addParams( outdir: params.outdir)

workflow ALIGN_DEDUP {
    take:
        /* a tuple with [sampleID, R1.fastq, R2.fastq] */ 
        fastqfiles
        genome_data
        
    main:
        BWA(fastqfiles, genome_data)
        MERGE_BAMS(BWA.out.bam_file.groupTuple())
        SAMBLASTER(MERGE_BAMS.out.merged_bam)
        SORT_BAM(SAMBLASTER.out.bam_files)

        bam_publish_paths = SORT_BAM.out.all_bams
            .map { it -> return it[0]+\
            "\t"+output_bam+"/"+it[1].getName()+\
            "\t"+output_bam+"/"+it[2].getName()+\
            "\t"+output_bam+"/"+it[3].getName() }
  
    emit:
        // tuple val(sample), file(bam), file(discbam), file(splitbam)
        all_bams = SORT_BAM.out.all_bams
        // tuple val(sample), file(bam), file(bai)
        main_bam = SORT_BAM.out.main_bam
        bam_paths = bam_publish_paths
}
