params.outdir = "lumpy_VCF"

workflow LUMPY_CALLS {
    take:
        input_data
        // tuple sampleID, file(bam_file), file(bai_file), file(split_bam), file(disc_bam)]
        exclude_regions
        // file  bed file of regions to exclude in Lumpy calls

    main:        
        estimate_readlength(input_data)
        estimate_libvalues(estimate_readlength.out) 
        LUMPY(estimate_libvalues.out, exclude_regions)

    emit:
        vcf_files = LUMPY.out.vcf_files
}

process estimate_readlength {
    label 'singlecore'

    input:
        tuple val(sampleID), file(bam_file), file(bai_file), file(split_bam), file(disc_bam)

    output:
        tuple val(sampleID), file(bam_file), file(bai_file), file(split_bam), file(disc_bam), stdout

    script:
    """
    samtools view $bam_file | tail -n+500000 | head -n200000 \
    | awk -F"\\t" '{sum += length(\$10)}; END {print sum / NR}' \
    | xargs printf '%.*f' 0
    """
}

process estimate_libvalues {
    label 'singlecore'
    publishDir "$params.outdir", mode: 'copy', pattern: "*.histo"

    input:
        tuple val(sampleID), file(bam_file), file(bai_file), file(split_bam), file(disc_bam), val(read_length)

    output:
        tuple val(sampleID), file(bam_file), file(bai_file), file(split_bam), file(disc_bam), file("${sampleID}.histo"), val(read_length), stdout

    script:
    """
    samtools view $bam_file | tail -n+100000 \
    | pairend_distro.py -r $read_length -X 4 -N 10000 -o ${sampleID}.histo \
    | tr "\t" "," | xargs printf '%s'
    """
}

process LUMPY {
    label 'lowcores'
    publishDir "$params.outdir", mode: 'copy', pattern: "*.vcf"

    input:
        tuple val(sampleID), file(bam_file), file(bai_file), file(split_bam), file(disc_bam), file(histo_file), val(read_length), val(lib_values)
        file(exclude_regions)

    output:
        tuple val(sampleID), file("${sampleID}_lumpy.vcf"), file(bam_file), file(bai_file), emit: vcf_files
        tuple val(sampleID), file(bam_file), file(bai_file), emit: bam_files

    script:
    """
    mkdir tmp_folder && \
    lumpy -t tmp_folder/lumpy_tmp -P -mw 4 -tt 0 \
    -x $exclude_regions \
    -pe id:${sampleID},bam_file:${disc_bam},histo_file:${histo_file},${lib_values},read_length:${read_length},min_non_overlap:${read_length},discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
    -sr id:${sampleID},bam_file:${split_bam},back_distance:10,weight:1,min_mapping_threshold:20 \
    > ${sampleID}_lumpy.vcf
    """
}