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
A vcf file named params.famid will be created containing SV across the family members
Based on BAM_FILES_PREFIX the pipeline search for .bam, .disc.bam, .split.bam 
*/
params.input = 'input_file.txt'

params.famid = "Family"
params.sv_folder = "SV"
params.ref="/well/gel/HICF2/ref/genomes/GRCh38/bwakit/hs38DH.fa"

output_lumpy = file("${params.sv_folder}/raw_data/lumpy_VCF")
output_cnvnator = file("${params.sv_folder}/raw_data/cnvnator_root")
output_famvcf = file("${params.sv_folder}/family_VCF")
output_singlesamples = file("${params.sv_folder}/samples_VCF")

output_lumpy.mkdirs()
output_cnvnator.mkdirs()
output_famvcf.mkdirs()
output_singlesamples.mkdirs()

workflow {
    log.info """\
    SV calls (lumpy + cnvnator) - N F   P I P E L I N E    
    ===================================================
    This pipeline works for hg38 only

    input file  (--input)      : ${params.input}
    genome      (--ref)        : ${params.ref}
    outdir      (--sv_folder)  : ${params.sv_folder}
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
    
    SV_vars_call(project, ref_genome, chrs_folder, sex_definitions)
}

workflow SV_vars_call {
    take:
        /* a tuple with [sampleID, file(bam_file), file(bai_file), file(split_bam), file(disc_bam)] */ 
        project
        ref_genome
        chrs_folder
        sex_definitions /* a list of values 'sampleID\tsex' where sex=1(M) / 2(F) */
        
    main:
        sex_file=sex_definitions.collectFile(name: 'sex.txt', newLine: true)
        
        estimate_readlength(project)
        estimate_libvalues(estimate_readlength.out) 
        lumpy(estimate_libvalues.out, ref_genome)
        cnvnator(project,chrs_folder)

        vcf_root_files = lumpy.out.bam_files.join(cnvnator.out)
        /* sampleID, bam, bai, root */

        concat_lumpy_calls(lumpy.out.vcfgt_files.collect())
        merge_lumpy_calls(concat_lumpy_calls.out)

        regenotype_family(vcf_root_files, merge_lumpy_calls.out)
        copynumber_family(regenotype_family.out, merge_lumpy_calls.out)

        make_family_vcf(copynumber_family.out.collect(), merge_lumpy_calls.out)
        prune_family_vcf(make_family_vcf.out)
        fix_family_vcf(prune_family_vcf.out)

        svtools_classify(fix_family_vcf.out, sex_file)
        quality_filter_sv(svtools_classify.out)
        convert_to_bed(quality_filter_sv.out)
        lowcomplexity_filter(quality_filter_sv.out,convert_to_bed.out)
        
    emit:
        lowcomplexity_filter.out
}

process estimate_readlength {
    label 'singlecore'
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'

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
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'

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
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'
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
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'
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

process concat_lumpy_calls {
    label 'singlecore'
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'

    input:
        file(lumpy_vcf)
    
    output:
        file('sorted.vcf')

    script:
    def lumpy_files = lumpy_vcf.collect{ "$it" }.join(" ")
    """
    svtools lsort $lumpy_files > sorted.vcf
    """
}

process merge_lumpy_calls {
    label 'singlecore'
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'

    input:
        file(sorted_vcf)
    
    output:
        tuple file('merged.vcf'), file('coordinates')

    script:
    """
    svtools lmerge -i $sorted_vcf -f 20 | grep -v "KN\\|JTFH\\|HLA" > merged.vcf
    create_coordinates -i merged.vcf -o coordinates
    """
}

process regenotype_family {
    label 'lowcores'
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'

    input:
        tuple val(sampleID), file(bam_file), file(bai_file), file(root_file)
        tuple file(merged_vcf), file(coordinates)

    output:
        tuple val(sampleID), file("${params.famid}_${sampleID}.gt.vcf"), file(root_file)
        
    script:
    """
    (grep "^#" $merged_vcf && grep -v "^#" $merged_vcf | awk '{OFS="\\t"};{ \$6="."; print }') \
    | svtools genotype -B $bam_file -l ${sampleID}.bam.json \
    | sed 's/PR...=[0-9\\.e,-]*\\(;\\)\\{0,1\\}\\(\\t\\)\\{0,1\\}/\\2/g' - \
    > ${params.famid}_${sampleID}.gt.vcf
    """
}

process copynumber_family {
    label 'lowcores'
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'
    publishDir "$output_singlesamples", mode: 'copy', pattern: '*.cn.vcf'

    input:
        tuple val(sampleID), file(lumpy_gt_file), file(root_file)
        tuple file(merged_vcf), file(coordinates)

    output:
        file("${params.famid}_${sampleID}.cn.vcf")

    script:
    """
    svtools copynumber --cnvnator cnvnator -w 100 \
    -s $sampleID \
    -r $root_file \
    -c $coordinates \
    -i $lumpy_gt_file \
    > ${params.famid}_${sampleID}.cn.vcf
    """
}

process make_family_vcf {
    label 'singlecore'
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'

    input: 
        file(cn_vcf)
        tuple file(merged_vcf), file(coordinates)
    
    output:
        file 'merged.sv.gt.cn.vcf'

    script:

    """
    ls *.cn.vcf > cnlist.txt
    svtools vcfpaste -m ${merged_vcf} -f cnlist.txt -q > merged.sv.gt.cn.vcf
    """
}

process prune_family_vcf {
    label 'singlecore'
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'

    input:
        file(family_vcf)
    
    output:
        file 'merged.sv.pruned.vcf.gz'

    script:

    """
    cat $family_vcf | svtools afreq | svtools vcftobedpe | svtools bedpesort \
    | svtools prune -s -d 100 -e \"AF\" | svtools bedpetovcf \
    | bgzip -c > temp.sv.pruned.vcf.gz
    (zgrep "#" temp.sv.pruned.vcf.gz && zgrep -v "#" temp.sv.pruned.vcf.gz \
    | awk -F"\t" '\$9 != "GT" && \$9 != "GT:CN" {print;}') \
    | bgzip -c > merged.sv.pruned.vcf.gz
    """  
}

process fix_family_vcf {
    label 'singlecore'
    container '/well/gel/HICF2/software/singularity/SV.annot_filter.sif'
    publishDir "$output_famvcf", mode: 'copy', pattern: '*.raw.vcf.gz'

    input:
        file(pruned_vcf)
    
    output:
        file "${params.famid}.raw.vcf.gz"

    script:
    """
    FixGQfield.py -v $pruned_vcf -o temp_GQfix.vcf
    sed 's/:\\./:0/g' temp_GQfix.vcf | bgzip -c > ${params.famid}.raw.vcf.gz
    """
}

process svtools_classify {
    label 'singlecore'
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'
    publishDir "$output_famvcf", mode: 'copy', pattern: '*.classify.vcf.gz'

    input:
        file(pruned_vcf)
        file(sex_file)
    
    output:
        file "${params.famid}.filters.classify.vcf.gz"

    script:

    """
    zcat $pruned_vcf \
    | svtools vcftobedpe  \
    | svtools varlookup -a stdin -b /resources/HICF2_training_vars.bedpe.gz -c HQ -d 50 \
    | svtools bedpetovcf \
    | svtools vcfsort \
    | bcftools view -i "INFO/HQ_AF > 0" \
    | bgzip -c > training.vars.vcf.gz

    zcat $pruned_vcf \
    | svtools classify \
    -g $sex_file \
    -a /resources/repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.GRCh38.sorted.bed.gz \
    -m large_sample \
    | bgzip -c > ${params.famid}.filters.classify.vcf.gz
    """     
}

process quality_filter_sv {
    label 'singlecore'
    container '/well/gel/HICF2/software/singularity/SV.annot_filter.sif'

    input:
        file(pruned_vcf)
    
    output:
        file 'temp_5.vcf'

    script:
    """
    zcat $pruned_vcf > temp_0.vcf
    
    bcftools filter -m + -s LowConfSmallDEL \
    -e 'INFO/SVTYPE == "DEL" && INFO/SVLEN >= -1000 && N_PASS(AS > 0) == 0' \
    temp_0.vcf > temp_1.vcf

    bcftools filter -m + -s INV_LowMSQ \
    -e 'INFO/SVTYPE == "INV" && INFO/MSQ < 150' \
    temp_1.vcf > temp_2.vcf

    bcftools filter -m + -s INV_LowLumpyEvidence \
    -e 'INFO/SVTYPE == "INV" && ((INFO/SR / INFO/SU) < 0.1 || (INFO/PE / INFO/SU) < 0.1)' \
    temp_2.vcf > temp_3.vcf

    bcftools filter -m + -s BND_LowMSQ \
    -e 'INFO/SVTYPE == "BND" && INFO/MSQ < 250' \
    temp_3.vcf > temp_4.vcf

    bcftools filter -m + -s VerySmallSV \
    -e '(INFO/SVTYPE == "DEL" && INFO/SVLEN > -50) || (INFO/SVTYPE == "DUP" && INFO/SVLEN < 50)' \
    temp_4.vcf > temp_5.vcf
    """    
}

process convert_to_bed {
    label 'singlecore'
    container '/well/gel/HICF2/software/singularity/SVpipeline.sif'

    input:
        file(temp_vcf)
    
    output:
        file 'temp_5.noBND.forOverlap.bed'

    script:
    """
    svtools vcftobedpe -i $temp_vcf -o temp_5.bedpe
    svtools bedpetobed12 -i temp_5.bedpe -o temp_5.bed
    grep -v "BND" temp_5.bed | tail -n+2 | cut -f1-4 | tr ";" "\t" \
    | cut -f1-5 | sed 's/ID=//g' > temp_5.noBND.forOverlap.bed
    """
}

process lowcomplexity_filter {
    label 'singlecore'
    container '/well/gel/HICF2/software/singularity/SV.annot_filter.sif'
    publishDir "$output_famvcf", mode: 'copy', pattern: '*.filters.vcf.gz'

    input:
        file(temp_vcf)
        file(bed_file)
    
    output:
        file "${params.famid}.filters.vcf.gz"

    script:
    """
    bedtools coverage -wa \
    -a $bed_file \
    -b /resources/hg38_CNV_ExcludeRegions+LowComplexity.bed \
    > overlapExcludeRegions.txt
    bgzip $temp_vcf
    
    FilterByValue.py \
    -v ${temp_vcf}.gz \
    -b overlapExcludeRegions.txt \
    -c 9 -i 5 -t LowComplexity -x 0.7 -l 1000 -o ${params.famid}.filters.raw.vcf
    sed 's/:\\./:0/g' ${params.famid}.filters.raw.vcf | bgzip -c > ${params.famid}.filters.vcf.gz
    """
}

