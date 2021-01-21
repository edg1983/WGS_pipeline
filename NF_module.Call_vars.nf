nextflow.enable.dsl=2

/*
In case running as separate workflow you need an input file.
Input file is a tab-separated text file without header and columns
1. SAMPLEID
2. BAM FILE LOCATION
*/
params.input = 'input_file.txt'
params.famid = 'family'
params.mode = "WGS"

params.vcf_folder = "VCF"
params.ref="/well/gel/HICF2/ref/genomes/GRCh38/bwakit/hs38DH.fa"

output_singlevcf = file("${params.vcf_folder}/single_VCF")
output_cohortvcf = file("${params.vcf_folder}/multi_VCF")
output_vcfreport = file("${output_cohortvcf}/reports")

output_singlevcf.mkdirs()
output_vcfreport.mkdirs()
output_cohortvcf.mkdirs()

workflow {
    log.info """\
    Deepvariant - N F   P I P E L I N E    
    ===================================
    This pipeline works for hg38 only

    mode        (--mode)       : ${params.mode}
    input file  (--input)      : ${params.input}
    genome      (--ref)        : ${params.ref}
    outdir      (--vcf_folder) : ${params.vcf_folder}
    ====================================
    """
    .stripIndent()
    
    designFile = file(params.input)
    project = Channel
                .from(designFile)
                .splitCsv(header: false, sep: '\t')
                .map{row -> return tuple(row[0], row[1], file(row[1]+".bai"))}
    small_vars_call(project)
}

workflow small_vars_call {
    take:
        /* a tuple with [sampleID, file(bam_file), file(bai_file)] */ 
        project
        
    main:
        deepvariant(file("$params.ref"),file("${params.ref}.fai"),project)
        multisample_vcf(deepvariant.out.gvcfs.collect())
        normalize_vcf(multisample_vcf.out,file("$params.ref"),file("${params.ref}.fai"))
        clean_deepvar_vcf(normalize_vcf.out)
        filter_vcf(clean_deepvar_vcf.out)
        
        bcftools_stats(filter_vcf.out.pass_calls)
        multiqc_vcf(bcftools_stats.out.collect())

    emit:
        pass_vcf = filter_vcf.out.pass_calls
        raw_vcf = filter_vcf.out.raw_calls
}

process deepvariant {
    label 'highcores_AVX'
    container '/well/gel/HICF2/software/singularity/deepvariant-1.0.0.sif'
    publishDir "$output_singlevcf", mode: 'copy', pattern: "*.vcf.gz*"
    publishDir "$output_singlevcf", mode: 'copy', pattern: "*.html"

    input:
        file(genome)
        file(genome_fai)
        tuple val(sampleID), file(bam_file), file(bai_file)

    output:
        path "*.g.vcf.gz", emit: gvcfs
        path "*.vcf.gz*", emit: vcfs
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

process multisample_vcf {
    label 'highcores'
    container '/well/gel/HICF2/software/singularity/glnexus-v1.2.6.simg'

    input:
        file(gvcf_files)

    output:
        file 'multisample.vcf.gz'

    script:
    def input_vcfs = gvcf_files.collect{ "$it" }.join(" ")

    """
    /usr/local/bin/glnexus_cli -m 64 -t 10 \
    -c DeepVariant${params.mode} \
    --dir GLnexusDB \
    $input_vcfs \
    | bcftools view - \
    | bgzip -@ 10 -c > multisample.vcf.gz
    """
}

process normalize_vcf {
    label 'singlecore'
    
    input:
        file(vcf_file)
        file(genome)
        file(genome_fai)

    output:
        path 'normalized_vcf.vcf.gz'

    script:
    """
    /well/gel/HICF2/software/vt/vt decompose -s $vcf_file \
    | /well/gel/HICF2/software/vt/vt normalize -r $genome - \
    | /apps/well/htslib/1.8/bin/bgzip -c > normalized_vcf.vcf.gz
    """
}

process clean_deepvar_vcf {
    label 'singlecore'
    beforeScript 'export C_INCLUDE_PATH=$C_INCLUDE_PATH:/well/gel/HICF2/software/additional_libs/include/;export LIBRARY_PATH=$LIBRARY_PATH:/well/gel/HICF2/software/additional_libs/lib/;export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/well/gel/HICF2/software/additional_libs/lib/'

    input:
        file(normalized_vcf)

    output:
        tuple file("raw.vcf.gz"), file("raw.vcf.gz.csi")

    script:
    """
    /apps/well/htslib/1.8/bin/tabix -p vcf -m 12 --csi $normalized_vcf
    
    /well/gel/HICF2/software/BRC_tools/VCF_fixes/VCF_fixes \
    -v $normalized_vcf -o raw.vcf.gz \
    --fixGQ --setMissing --fixhalfGT --computeAB

    /apps/well/htslib/1.8/bin/tabix -p vcf -m 12 --csi raw.vcf.gz
    """
}

process filter_vcf {
    label 'singlecore'
    publishDir "$output_cohortvcf", mode: 'copy', pattern: '*.raw.vcf.gz*'
    publishDir "$output_cohortvcf", mode: 'copy', pattern: '*.PASS.vcf.gz*'

    input:
        tuple file(normalized_vcf), file(vcf_index)

    output:
        tuple file("${params.famid}.raw.vcf.gz"), file("${params.famid}.raw.vcf.gz.csi"), emit: raw_calls
        tuple file("${params.famid}.PASS.vcf.gz"), file("${params.famid}.PASS.vcf.gz.csi"), emit: pass_calls

    script:
    """
    /well/gel/HICF2/software/bcftools/bcftools-1.10.2/bcftools filter -m + -s LowQual \
    -e 'QUAL < 20 || N_PASS(GQ >= 20) == 0' \
    -Oz -o ${params.famid}.raw.vcf.gz \
    $normalized_vcf
    /apps/well/htslib/1.8/bin/tabix -p vcf -m 12 --csi ${params.famid}.raw.vcf.gz

    zgrep "#\\|PASS" ${params.famid}.raw.vcf.gz \
    | /apps/well/htslib/1.8/bin/bgzip -c > ${params.famid}.PASS.vcf.gz
    /apps/well/htslib/1.8/bin/tabix -p vcf -m 12 --csi ${params.famid}.PASS.vcf.gz
    """
}

process bcftools_stats {
    label 'singlecore'
    publishDir "$output_vcfreport", mode: 'copy'
    
    input:
        tuple file(vcf_file), file(csi_index)

    output:
        file "${vcf_file}.stats"

    script:
    """
    /well/gel/HICF2/software/bcftools/bcftools-1.10.2/bcftools stats -s- $vcf_file > ${vcf_file}.stats
    """
}

process multiqc_vcf {
    container '/well/gel/HICF2/software/singularity/QCbam-v1.0.sif'
    label 'singlecore'
    publishDir "$output_vcfreport", mode: 'copy'

    input:
    file('*')

    output:
    file 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}