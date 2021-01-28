nextflow.enable.dsl=2

/*
In case running as separate workflow you need an input file.
Input file is a tab-separated text file without header and columns
1. SAMPLEID
2. BAM FILE LOCATION
*/
params.input = 'input_file.txt'
params.mode = "WGS"
params.process_vcf = false
params.merge_gvcfs = false

params.out_folder = "VCF"
params.ref="/well/gel/HICF2/ref/genomes/GRCh38/bwakit/hs38DH.fa"

output_vcfs = file("${params.out_folder}")
output_htmls = file("${params.out_folder}/deepvariant_html")
output_htmls.mkdirs()

if (params.process_vcf) {
    include { process_deepvariant_vcf } from "./NF_module.process_deepvariant_VCF.nf" addParams(out_folder: "${params.out_folder}/processed_vcfs")
}
if (params.merge_gvcfs) {
    include { merge_gvcfs } from "./NF_module.merge_gvcfs.nf" addParams(out_folder: "${params.out_folder}/merged_vcf")
}

workflow {
    log.info """\
    Deepvariant - N F   P I P E L I N E    
    ===================================
    This pipeline works for hg38 only

    mode            (--mode)       : ${params.mode}
    input file      (--input)      : ${params.input}
    genome          (--ref)        : ${params.ref}
    outdir          (--out_folder) : ${params.out_folder}
    gvcf merge      (--merge_gvcfs): ${params.merge_gvcfs}
    vcf processing  (--process_vcf): ${params.process_vcf}
    ====================================
    """
    .stripIndent()
    
    designFile = file(params.input)
    project = Channel
                .from(designFile)
                .splitCsv(header: false, sep: '\t')
                .map{row -> return tuple(row[0], file(row[1]), file(row[1]+".bai"))}
    ref_genome = tuple(file("${params.ref}"), file("${params.ref}.fai"))

    main:
        deepvariant(ref_genome,project)

        if (params.merge_gvcfs) {
            merge_gvcfs(deepvariant.out.vcfs)
            vcfs_to_process = merge_gvcfs.out
        } else {
            vcfs_to_process = deepvariant.out.vcfs
        }
        
        if (params.process_vcf) {
            process_deepvariant_vcf(vcfs_to_process, ref_genome)
        }
}

process deepvariant {
    label 'highcores_AVX'
    publishDir "$output_vcfs", mode: 'copy', pattern: "*.vcf.gz*"
    publishDir "$output_htmls", mode: 'copy', pattern: "*.html"

    input:
        tuple file(genome), file(genome_fai)
        tuple val(sampleID), file(bam_file), file(bai_file)

    output:
        path "*.g.vcf.gz", emit: gvcfs
        tuple val("${sampleID}"), file("${sampleID}.vcf.gz"), file("${sampleID}.vcf.gz.tbi"), emit: vcfs
        path "*.html", emit: htmls

    script:
    STDCHRS = (1..22).collect { "chr$it"}.join(' ')
    STDCHRS = "$STDCHRS chrX chrY chrM"

    """
    mkdir temp_dir && \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${params.mode} --ref="$genome" \
    --reads="$bam_file" \
    --output_vcf="${sampleID}.vcf.gz" \
    --output_gvcf="${sampleID}.g.vcf.gz" \
    --intermediate_results_dir="temp_dir" \
    --num_shards=10 \
    --regions="$STDCHRS"
    """
}