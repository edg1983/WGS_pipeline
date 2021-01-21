nextflow.enable.dsl=2

params.input = 'input_file.txt'
params.famid = 'family'
params.ped = 'input.ped'
params.output = 'output'
params.qc_folder = "${params.output}/QC"
params.bam_folder = "${params.output}/BAM"
params.vcf_folder = "${params.output}/VCF"
params.sv_folder = "${params.output}/SV"
params.roh_folder = "${params.output}/ROH"
params.exphunter_folder = "${params.output}/ExpHunter"
params.mode = "WGS"

params.help = false
if (params.help) {
println """\
    HICF2 full pipeline - PARAMETERS    
    ===========================================
    NB. This pipeline works for hg38 only

    --input input_file.txt  :   A tab-delimited text file describing samples
                                sampleID, fastq_R1.fastq.gz, fastq_R2.fastq.gz
    --ped ped_file.ped      :   A standard ped files describing all samples
    --ref genome.fa         :   Genome reference file
                                Note that bwa index files are supposed to be available
                                in the same folder with default extensions
    --mode WGS/WES          :   Set running mode for WES or WGS
    --output outdir         :   Output folder. A series of sub-folder will be created
    """
    .stripIndent()
exit 1
}

log.info """\
    HICF2 full pipeline - N F   P I P E L I N E    
    ===========================================
    This pipeline works for hg38 only

    mode        : ${params.mode}
    input file  : ${params.input}
    ped file    : ${params.ped}
    genome      : ${params.ref}
    chrs folder : ${params.chrs_folder}
    outdir      : ${params.output}
    ===========================================
    """
    .stripIndent()

include { align_dedup } from "./NF_module.Align_dedup.nf"
include { QC_bam } from "./NF_module.QC_bam.nf"
include { small_vars_call } from "./NF_module.Call_vars.nf"
include { SV_vars_call } from "./NF_module.SV_calls.nf"
include { ROH_detection } from "./NF_module.bcftools_ROH.nf"
include { ExpansionHunter } from "./NF_module.ExpansionHunter.nf"

output_dir = file(params.output)
output_dir.mkdirs()

workflow {
    pedfile = file(params.ped)
    
    designFile = file(params.input)
    fastqfiles = Channel
                   .from(designFile)
                   .splitCsv(header: false, sep: '\t')
                   .map{row -> return tuple(row[0], file(row[1]), file(row[2]))}
    fam_structure = Channel
                   .from(pedfile)
                   .splitCsv(header: false, sep: '\t')
                   .map{row -> return tuple(row[1], row[0])}
    sex_definitions = Channel
                   .from(pedfile)
                   .splitCsv(header: false, sep: '\t')
                   .map{row -> return row[1]+"\t"+row[4]}
    chrs_folder = file("${params.chrs_folder}")
    genome_file = file("${params.ref}")
    genome_base = genome_file.getBaseName()
    genome_fasta = genome_file.getName()
    genome_folder = genome_file.getParent()
    genome_data = tuple("$genome_fasta", file("${genome_folder}/${genome_base}.*"))
    ref_genome = tuple(file("${params.ref}"), file("${params.ref}.fai"))
    AF_file = tuple(file("${params.roh_AFfile}"), file("${params.roh_AFfile}.csi"))
    genetic_maps = file("${params.roh_gmaps}/*.txt")
    variant_catalog = file(params.exphunter_catalog)

    align_dedup(fastqfiles, genome_data)
    QC_bam(align_dedup.out.bamfiles, file(params.regions), file(params.somalier_data), ref_genome)
    
    small_vars_call(align_dedup.out.bamfiles)

    SV_vars_call(align_dedup.out.allbams, ref_genome, chrs_folder, sex_definitions)

    ROH_detection(small_vars_call.out.pass_vcf, AF_file, genetic_maps)

    exp_hunter_sex = Channel
                .from(pedfile)
                .splitCsv(header: false, sep: '\t')
                .map{row -> return tuple(row[1], row[4].replace("1", "male").replace("2","female"))}
    exp_hunter_input = exp_hunter_sex.combine(align_dedup.out.bamfiles, by: 0)
    ExpansionHunter(exp_hunter_input, ref_genome, variant_catalog)             
}

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