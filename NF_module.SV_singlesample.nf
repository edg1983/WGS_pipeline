nextflow.enable.dsl=2

/*
Additional params from config
chrs_folder: a folder containing chrs.fa sequences (like chr1.fa ...)

In case running as separate workflow you need an input file describing
samples included in the family / cohort.
Input file is a tab-separated text file without header and columns
1. SAMPLEID
2. SEX (1=male, 2=female) 
3. BAM_FILES_PREFIX
Based on BAM_FILES_PREFIX the pipeline search for .bam, .disc.bam, .split.bam 
*/
params.input = 'input_file.txt'
params.merge_vcfs = false

params.out_folder = "SV"
params.ref="/well/gel/HICF2/ref/genomes/GRCh38/bwakit/hs38DH.fa"

output_lumpy = file("${params.out_folder}/lumpy_VCF")
output_cnvnator = file("${params.out_folder}/cnvnator_root")

if (params.merge_vcfs) {
    include { SV_mergesamples } from "./NF_module.SV_mergesamples.nf" addParams(out_folder: "${params.out_folder}/merged_vcf")
}

output_lumpy.mkdirs()
output_cnvnator.mkdirs()

workflow {
    log.info """\
    SV calls (lumpy + cnvnator) - N F   P I P E L I N E    
    ===================================================
    This pipeline works for hg38 only

    input file  (--input)      : ${params.input}
    genome      (--ref)        : ${params.ref}
    merge vcfs  (--merge_vcfs) : ${params.merge_vcfs}
    outdir      (--out_folder) : ${params.out_folder}
    ====================================================
    """
    .stripIndent()
    
    designFile = file(params.input)
    project = Channel
                .from(designFile)
                .splitCsv(header: false, sep: '\t')
                .map{row -> return tuple(
                    row[0], 
                    file(row[2]+".bam"), 
                    file(row[2]+".bam.bai"),
                    file(row[2]+".split.bam"),
                    file(row[2]+".disc.bam"))}
    sex_definitions = Channel
                .from(designFile)
                .splitCsv(header: false, sep: '\t')
                .map{row -> return row[0]+"\t"+row[1]}

    ref_genome = tuple(file("${params.ref}"), file("${params.ref}.fai"))
    chrs_folder = file("${params.chrs_folder}")
    
    SV_singlesample(project, ref_genome, chrs_folder)

    if (params.merge_vcfs) {
        SV_mergesamples(SV_singlesample.out.root_files, SV_singlesample.out.vcfgt_files, sex_definitions)
    }
}

workflow SV_singlesample {
    take:
        /* a tuple with [sampleID, file(bam_file), file(bai_file), file(split_bam), file(disc_bam)] */ 
        project
        ref_genome
        chrs_folder
        
    main:        
        estimate_readlength(project)
        estimate_libvalues(estimate_readlength.out) 
        lumpy(estimate_libvalues.out, ref_genome)
        cnvnator(project,chrs_folder)

    emit:
        root_files = lumpy.out.bam_files.join(cnvnator.out) /* sampleID, bam, bai, root */
        vcfgt_files = lumpy.out.vcfgt_files
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
    publishDir "$output_lumpy", mode: 'copy', pattern: "*.histo"

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

process lumpy {
    label 'lowcores'
    publishDir "$output_lumpy", mode: 'copy', pattern: "*.gt.vcf.gz"

    input:
        tuple val(sampleID), file(bam_file), file(bai_file), file(split_bam), file(disc_bam), file(histo_file), val(read_length), val(lib_values)
        tuple file(genome_fasta), file(genome_fai)

    output:
        path "${sampleID}_lumpy.gt.vcf.gz", emit: vcfgt_files
        tuple val(sampleID), file(bam_file), file(bai_file), emit: bam_files

    script:
    """
    mkdir tmp_folder && \
    lumpy -t tmp_folder/lumpy_tmp -P -mw 4 -tt 0 \
    -x /resources/hg38_CNV_ExcludeRegions.bed \
    -pe id:${sampleID},bam_file:${disc_bam},histo_file:${histo_file},${lib_values},read_length:${read_length},min_non_overlap:${read_length},discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
    -sr id:${sampleID},bam_file:${split_bam},back_distance:10,weight:1,min_mapping_threshold:20 \
    > lumpy.vcf
    svtyper --verbose -B $bam_file -T $genome_fasta -i lumpy.vcf | bgzip -c > ${sampleID}_lumpy.gt.vcf.gz
    """
}

process cnvnator {
    label 'lowcores'
    publishDir "$output_cnvnator", mode: 'copy', pattern: "*.root"

    input:
        tuple val(sampleID), file(bam_file), file(bai_file), file(split_bam), file(disc_bam)
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