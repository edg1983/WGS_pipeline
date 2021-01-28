nextflow.enable.dsl=2

/*
In case running as separate workflow you need an input file.
Input file is a text file with a list of g.vcf.gz to be merged, one per line
*/
params.input = 'input_file.txt'
params.cohort_id = 'family'
params.mode = "WGS"

params.out_folder = "VCF"

output_vcf = file("${params.out_folder}")
output_vcf.mkdirs()

workflow {
    log.info """\
    GLNexus gVCFs merge - N F   P I P E L I N E    
    ===================================
    This pipeline works for hg38 only

    mode        (--mode)       : ${params.mode}
    input file  (--input)      : ${params.input}
    cohort id   (--cohort_id)  : ${params.cohort_id}
    outdir      (--out_folder) : ${params.out_folder}
    ====================================
    """
    .stripIndent()
    
    designFile = file(params.input)
    gvcf_files = Channel
                .from(designFile)
                .eachLine() { line -> 
                    return tuple(val(file(line).getSimpleName()), file(line), file(line+".tbi"))
                    }
    merge_gvcfs(gvcf_files.collect())
}

process merge_gvcfs {
    label 'highcores'
    publishDir "$output_vcf", mode: 'copy', pattern: "${params.cohort_id}.vcf.gz*"

    input:
        file(gvcf_files)

    output:
        tuple val("${params.cohort_id}"), file("${params.cohort_id}.vcf.gz"), file("${params.cohort_id}.vcf.gz.csi")

    script:
    def input_vcfs = gvcf_files.collect{ "$it" }.join(" ")

    """
    glnexus_cli -m 64 -t 10 \
    -c DeepVariant${params.mode} \
    --dir GLnexusDB \
    $input_vcfs \
    | bcftools view - \
    | bgzip -@ 10 -c > ${params.cohort_id}.vcf.gz

    tabix -p vcf -m 12 --csi ${params.cohort_id}.vcf.gz
    """
}