// merge multiple BAMs in a single BAM named sample.bam
// this do not publish files since one usually pass output to sort or dedup

process MERGE_BAMS {
    label 'lowcores'

    input:
    /*  A tuple of files to be merged per sample
        like: (sample_id, [part1.bam, part2.bam]) */
    tuple val(sample), file(myfiles)

    output:
    tuple val(sample), file('merged.bam'), emit: merged_bam

    script:
    def input_args = myfiles.collect{ "$it" }.join(" ")
    
    """
    samtools merge -n -f -@ 5 merged.bam $input_args
    """
}