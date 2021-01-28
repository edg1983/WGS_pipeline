nextflow.enable.dsl=2

/*
In case running as separate workflow you need an input file.
Input file is a text file with .vcf.gz files one per line
*/
params.input = 'input_file.txt'
params.out_folder = 'processed_vcfs'
params.ref="/well/gel/HICF2/ref/genomes/GRCh38/bwakit/hs38DH.fa"

output_vcffiles = file("${params.out_folder}")
output_reports = file("${params.out_folder}/vcf_reports")
output_reports.mkdirs()

workflow {
    log.info """\
    Deepvariant - N F   P I P E L I N E    
    ===================================
    This pipeline works for hg38 only

    input file  (--input)      : ${params.input}
    genome      (--ref)        : ${params.ref}
    outdir      (--out_folder) : ${params.vcf_folder}
    ====================================
    """
    .stripIndent()
    
    designFile = file(params.input)
    vcf_files = Channel
                .from(designFile)
                .eachLine() { line -> 
                    return tuple(val(file(line).getSimpleName()), file(line), file(line+".tbi"))
                    }
    ref_genome = tuple(file("${params.ref}"), file("${params.ref}.fai"))
    process_deepvariant_vcf(vcf_files, ref_genome)
}

workflow process_deepvariant_vcf {    
    take:
        vcf_files /* a tuple [val("prefix for output"), file(vcf_file), file(vcf_index)] */
        ref_genome  /* a tuple with [file(ref.fasta), file(ref.fasta.fai)] */ 

    main:
        clean_and_filter_vcf(vcf_files,ref_genome)
        multiqc_vcf(clean_and_filter_vcf.out.stats.collect())

    emit:
        pass_vcf = clean_and_filter_vcf.out.pass_vcf
        raw_vcf = clean_and_filter_vcf.out.processed_vcf
}

process clean_and_filter_vcf {
    label 'singlecore'
    
    publishDir "$output_reports", mode: 'copy', pattern: "*.stats"
    publishDir "$output_vcffiles", mode: 'copy', pattern: "${prefix}*.vcf.gz*"

    input:
        tuple val(prefix), file(vcf_file), file(vcf_index)
        tuple file(genome), file(genome_fai)

    output:
        tuple file("${prefix}.PASS.vcf.gz"), file("${prefix}.PASS.vcf.gz.csi"), emit: pass_vcf
        tuple file("${prefix}.processed.vcf.gz"), file("${prefix}.processed.vcf.gz.csi"), emit: processed_vcf
        path "${prefix}.PASS.vcf.stats", emit: stats

    script:
    """
    vt decompose -s $vcf_file \
    | vt normalize -r $genome - \
    | bgzip -c > normalized.vcf.gz

    tabix -p vcf -m 12 --csi normalized.vcf.gz
    
    VCF_fixes \
    -v normalized.vcf.gz -o raw.vcf.gz \
    --fixGQ --setMissing --fixhalfGT --computeAB
    tabix -p vcf -m 12 --csi raw.vcf.gz

    bcftools filter -m + -s LowQual \
    -e 'QUAL < 20 || N_PASS(GQ >= 20) == 0' \
    -Oz -o ${prefix}.processed.vcf.gz \
    raw.vcf.gz
    tabix -p vcf -m 12 --csi ${prefix}.processed.vcf.gz

    zgrep "#\\|PASS" ${prefix}.processed.vcf.gz \
    | bgzip -c > ${prefix}.PASS.vcf.gz
    tabix -p vcf -m 12 --csi ${prefix}.PASS.vcf.gz

    bcftools stats -s- ${prefix}.PASS.vcf.gz > ${prefix}.PASS.vcf.stats
    """
}

process multiqc_vcf {
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