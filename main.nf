nextflow.enable.dsl=2

// Define allowed options
def allowed_ops         = ['align', 'call_variants']
def allowed_modes       = ['WGS','WES']

// Default values for main output params are defined in config

// Print help message when --help is used
params.help = false
if (params.help) {
println "help"

println """\
    HICF2 full pipeline - PARAMETERS    
    ===========================================
    NB. This pipeline works for hg38 only

    --operation             :   Operation to run. Values: ${allowed_ops.join(', ')}
    --input input_file.txt  :   A tab-delimited text file describing input files
                                Exact format depend on the operation
    --ped ped_file.ped      :   A standard ped file describing all samples
    --ref genome.fa         :   Genome reference file
                                Note that bwa index files are supposed to be available
                                in the same folder with default extensions
    --mode WGS/WES          :   Set running mode for WES or WGS
    --outdir outdir         :   Output folder. A series of sub-folder will be created
    --cohort_id cohort      :   A cohort id that is used in merged files
    """
    .stripIndent()

exit 1
}

if (!params.cohort_id && params.operation == "call_variants") {
    exit 1, "Cohort ID not set, please use --cohort_id"
}

log.info """\
    =============================================
    WGS analysis pipeline - N F   P I P E L I N E    
    =============================================
    This pipeline works for hg38 only

    operation   : ${params.operation}
    mode        : ${params.mode}
    input file  : ${params.input}
    ped file    : ${params.ped}
    cohort id   : ${params.cohort_id}
    genome      : ${params.ref}
    outdir      : ${params.outdir}
    =============================================
    """
    .stripIndent()

// Check requested operation/mode/samplesize is allowed
if (!allowed_ops.contains(params.operation)) {
    exit 1, "Invalid operation: ${params.operation}. Valid options: ${allowed_ops.join(', ')}"
}

if (!allowed_modes.contains(params.mode)) {
    exit 1, "Invalid mode: ${params.mode}. Valid options: ${allowed_modes.join(', ')}"
}

// Check input file exists 
file(params.input, checkIfExists: true)

// Make main output directory
output_dir = file(params.outdir)
output_dir.mkdirs()

// Get number of samples in the input file and set params accordingly
projectFile = file(params.input)

samplesize = projectFile
    .splitCsv(header: false, sep: '\t')
    .collect{ it[0] }
    .unique()
    .size()

if (params.singlesamples) {
    params.samplesize = 1
    log.warn """\
    You requested a single sample analysis
    All samples will be analyzed separately
    =============================================
    """
    .stripIndent()
}

params.samplesize = params.singlesamples ? 1 : samplesize 

log.info """\
    Information from input file: ${params.input}
    $params.samplesize samples found in your input file
    =============================================
    """
    .stripIndent()

// If calling variants with mode WES give a warning about unreliable SV/ROH detection
if (params.operation == "call_variants" && params.mode == "WES") {
    log.info """\
    WARNING: WES mode selected
    ROH and SV detection are not accurate on WES data
    Consider these results with caution
    =============================================
    """
    .stripIndent()
}

// Load workflows as requested by operation
if (params.operation == "align") {
    include { ALIGN_DEDUP } from "./modules/align_dedup/align_dedup" addParams(outdir: "${params.bam_folder}")
    include { QC_BAM } from "./modules/qc_bam/qc_bam" addParams(outdir: "${params.qc_folder}")
}

if (params.operation == "call_variants") {
    include { SMALLVARS_DEEPVAR } from "./modules/deepvariant/smallvars_deepvariant" addParams(outdir: params.vcf_folder, mode: params.mode)
    include { QC_VCF } from "./modules/qc_vcf/qc_vcf" addParams(outdir: "${params.vcf_folder}/vcf_QC")

    include { SV_CALL_SVTOOLS } from "./modules/sv_svtools/sv_svtools" addParams(outdir: "${params.sv_folder}")
    include { SV_FILTER_SMALLSAMPLE; SV_FILTER_LARGESAMPLE } from "./modules/sv_svtools/sv_filter" addParams(outdir: "${params.sv_folder}")

    include { ROH_DETECTION } from "./modules/roh_detection/detect_roh" addParams(outdir: "${params.roh_folder}")
    include { EXPANSION_HUNTER_CALLS } from "./modules/expansion_hunter/expansion_hunter" addParams(outdir: "${params.exphunter_folder}")
}

workflow {
    // ================
    // ALIGN OPERATION 
    // ================
    if (params.operation == 'align') {
        // Check input path parameters to see if they exist
        checkPathParamList = [
            params.ref, params.ped,
            params.regions, params.somalier_data
        ]
        for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

        // Set input channels
        fastqfiles = Channel
                   .from(projectFile)
                   .splitCsv(header: false, sep: '\t')
                   .map{row -> return tuple(row[0], file(row[1]), file(row[2]))}

        genome_file = file(params.ref)
        genome_base = genome_file.getBaseName()
        genome_fasta = genome_file.getName()
        genome_folder = genome_file.getParent()
        genome_data = tuple("$genome_fasta", file("${genome_folder}/${genome_base}.*"))
        ref_genome = tuple(file("${params.ref}"), file("${params.ref}.fai"))
        
        //WORKFLOW
        ALIGN_DEDUP(fastqfiles, genome_data)
        ALIGN_DEDUP.out.bam_paths
            .collectFile(name: 'bam_files.txt', newLine: true, storeDir: "${params.outdir}")
        QC_BAM(ALIGN_DEDUP.out.main_bam, file(params.regions), file(params.somalier_data), file(params.ped), ref_genome)    
    }

    // ================
    // CALL VARIANTS OPERATION 
    // ================  
    if (params.operation == 'call_variants') {
        // Check input path parameters to see if they exist
        checkPathParamList = [
            params.input, params.ref, params.ped,
            params.sv_exclude_regions, params.sv_lowcomplexity, params.sv_mei_regions,
            params.sv_mei_regions, params.sv_training_bedpe,
            params.chrs_folder, params.somalier_data,
            params.roh_AFfile, params.roh_gmaps,
            params.exphunter_catalog
        ]
        for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    
        //Set input channels
        pedFile = file(params.ped)
        chrs_folder = file(params.chrs_folder)
        ref_genome = tuple(file("${params.ref}"), file("${params.ref}.fai"))
        regular_bam = Channel
                   .from(projectFile)
                   .splitCsv(header: false, sep: '\t')
                   .map{row -> return tuple(row[0], file(row[1]), file(row[1]+".bai"))}
        discsplit_bam = Channel
                   .from(projectFile)
                   .splitCsv(header: false, sep: '\t')
                   .map{row -> return tuple(row[0], file(row[1]), file(row[1]+".bai"), file(row[2]), file(row[3]))}
        AF_file = tuple(file("${params.roh_AFfile}"), file("${params.roh_AFfile}.csi"))
        genetic_maps = file("${params.roh_gmaps}/*.txt")
        variant_catalog = file(params.exphunter_catalog)
    
        exclude_regions = file("${params.sv_exclude_regions}")
        lowcomplexity_regions = tuple(file("${params.sv_lowcomplexity}"), file("${params.sv_lowcomplexity}.tbi"))
        mei_regions = tuple(file("${params.sv_mei_regions}"), file("${params.sv_mei_regions}.tbi"))
        training_bedpe = file("${params.sv_training_bedpe}")
        
        sex_definitions = Channel
                .from(pedFile)
                .splitCsv(header: false, sep: '\t')
                .map{row -> return row[1]+"\t"+row[4]}

        //WORKFLOW
        //Small variants
        SMALLVARS_DEEPVAR(ref_genome, regular_bam)
        QC_VCF(SMALLVARS_DEEPVAR.out.pass_vcf)

        //Structural variants
        SV_CALL_SVTOOLS(discsplit_bam, regular_bam, ref_genome, exclude_regions, chrs_folder)
        
        if (params.samplesize >= 30) {
            SV_FILTER_LARGESAMPLE(SV_CALL_SVTOOLS.out,sex_definitions,lowcomplexity_regions,mei_regions)
        } else {
            SV_FILTER_SMALLSAMPLE(SV_CALL_SVTOOLS.out,sex_definitions,lowcomplexity_regions,mei_regions,training_bedpe)
        }

        // ROH regions
        ROH_DETECTION(SMALLVARS_DEEPVAR.out.pass_vcf, AF_file, genetic_maps)

        // STR repeats
        exp_hunter_sex = Channel
                    .from(pedFile)
                    .splitCsv(header: false, sep: '\t')
                    .map{row -> return tuple(row[1], row[4].replace("1", "male").replace("2","female"))}
        exp_hunter_input = exp_hunter_sex.combine(regular_bam, by: 0)
        EXPANSION_HUNTER_CALLS(exp_hunter_input, ref_genome, variant_catalog)  
    } 
}

    //    .groupTuple(by: [0])
    //    .map { it ->  [ it[0], it[1].flatten() ] }
    
workflow.onComplete { 
	log.info ( workflow.success ? """\
    
    PIPELINE COMPLETED!    
    ============================================
    Aligned BAM files   --> ${params.bam_folder}
    BAM QC files        --> ${params.qc_folder}
    Small variants      --> ${params.vcf_folder}
    SV variants         --> ${params.sv_folder}
    ExpHunter results   --> ${params.exphunter_folder}
    ROH regions         --> ${params.roh_folder}
    """.stripIndent() : "Oops .. something went wrong" )
}