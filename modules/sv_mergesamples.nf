nextflow.enable.dsl=2

/*
Additional params from config
chrs_folder: a folder containing chrs.fa sequences (like chr1.fa ...)

In case running as separate workflow you need an input file describing
samples included in the family / cohort.
Input file is a tab-separated text file without header and columns
1. SAMPLEID
2. SEX (1=male, 2=female)
3. BAM FILE (.bam)
4. LUMPY GENOTYPED VCF (.vcf.gz)
5. CNVNATOR ROOT FILE (.root)
A vcf file named params.cohort_id will be created containing SV across the cohort members
*/
params.input = 'input_file.txt'

params.cohort_id = "Family"
params.out_folder = "SV"

output_famvcf = file("${params.out_folder}")
output_singlesamples = file("${params.out_folder}/single_gtVCF")

output_famvcf.mkdirs()
output_singlesamples.mkdirs()

workflow {
    log.info """\
    SV calls (lumpy + cnvnator) - N F   P I P E L I N E    
    ===================================================
    This pipeline works for hg38 only

    input file  (--input)      : ${params.input}
    cohort id   (--cohort_id)  : ${params.cohort_id}
    outdir      (--out_folder) : ${params.out_folder}
    ====================================================
    """
    .stripIndent()
    
    designFile = file(params.input)
    root_files = Channel
                .from(designFile)
                .splitCsv(header: false, sep: '\t')
                .map{row -> return tuple(
                    row[0], 
                    file(row[2]), 
                    file(row[2]+".bai"),
                    file(row[4]))}
    vcfgt_files = Channel
                .from(designFile)
                .splitCsv(header: false, sep: '\t')
                .map{row -> return file(row[3])}
    sex_definitions = Channel
                .from(designFile)
                .splitCsv(header: false, sep: '\t')
                .map{row -> return row[0]+"\t"+row[1]}
    
    lowcomplexity_regions = tuple(file("${params.sv_lowcomplexity}"), file("${params.sv_lowcomplexity}.tbi"))
    mei_regions = tuple(file("${params.sv_mei_regions}"), file("${params.sv_mei_regions}.tbi"))
    training_bedpe = file("${params.sv_training_bedpe}")
    
    SV_mergesamples(root_files, vcfgt_files, sex_definitions,lowcomplexity_regions,mei_regions,training_bedpe)
}

workflow SV_mergesamples {
    take:
        /* root_files: sampleID, bam, bai, root */
        /* vcfgt_files: genotyped lumpy files .vcf.gz */
        root_files
        vcfgt_files
        sex_definitions /* a list of values 'sampleID\tsex' where sex=1(M) / 2(F) */
        lowcomplexity_regions
        mei_regions
        training_bedpe

    main:
        sex_file=sex_definitions.collectFile(name: 'sex.txt', newLine: true)

        concat_lumpy_calls(vcfgt_files.collect())
        merge_lumpy_calls(concat_lumpy_calls.out)

        regenotype_family(root_files, merge_lumpy_calls.out)
        copynumber_family(regenotype_family.out, merge_lumpy_calls.out)

        make_family_vcf(copynumber_family.out.collect(), merge_lumpy_calls.out)
        prune_family_vcf(make_family_vcf.out)
        fix_family_vcf(prune_family_vcf.out)

        svtools_classify(fix_family_vcf.out, sex_file, training_bedpe, mei_regions)
        quality_filter_sv(svtools_classify.out)
        convert_to_bed(quality_filter_sv.out)
        lowcomplexity_filter(quality_filter_sv.out,convert_to_bed.out, lowcomplexity_regions)
        
    emit:
        lowcomplexity_filter.out
}

process concat_lumpy_calls {
    label 'singlecore'

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

    input:
        tuple val(sampleID), file(bam_file), file(bai_file), file(root_file)
        tuple file(merged_vcf), file(coordinates)

    output:
        tuple val(sampleID), file("${params.cohort_id}_${sampleID}.gt.vcf"), file(root_file)
        
    script:
    """
    (grep "^#" $merged_vcf && grep -v "^#" $merged_vcf | awk '{OFS="\\t"};{ \$6="."; print }') \
    | svtools genotype -B $bam_file -l ${sampleID}.bam.json \
    | sed 's/PR...=[0-9\\.e,-]*\\(;\\)\\{0,1\\}\\(\\t\\)\\{0,1\\}/\\2/g' - \
    > ${params.cohort_id}_${sampleID}.gt.vcf
    """
}

process copynumber_family {
    label 'lowcores'
    publishDir "$output_singlesamples", mode: 'copy', pattern: '*.cn.vcf'

    input:
        tuple val(sampleID), file(lumpy_gt_file), file(root_file)
        tuple file(merged_vcf), file(coordinates)

    output:
        file("${params.cohort_id}_${sampleID}.cn.vcf")

    script:
    """
    svtools copynumber --cnvnator cnvnator -w 100 \
    -s $sampleID \
    -r $root_file \
    -c $coordinates \
    -i $lumpy_gt_file \
    > ${params.cohort_id}_${sampleID}.cn.vcf
    """
}

process make_family_vcf {
    label 'singlecore'

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

    input:
        file(family_vcf)
    
    output:
        file 'merged.sv.pruned.vcf.gz'

    script:

    """
    cat $family_vcf | svtools afreq | svtools vcftobedpe | svtools bedpesort \
    | svtools prune -s -d 100 -e \"AF\" | svtools bedpetovcf \
    | bgzip -c > merged.sv.pruned.vcf.gz
    """  
}

process fix_family_vcf {
    label 'singlecore'
    publishDir "$output_famvcf", mode: 'copy', pattern: '*.raw.vcf.gz'

    input:
        file(pruned_vcf)
    
    output:
        file "${params.cohort_id}.raw.vcf.gz"

    script:
    """
    Fix_svtools_VCF.py -v $pruned_vcf -o ${params.cohort_id}.raw.vcf
    bgzip ${params.cohort_id}.raw.vcf
    """
}

process svtools_classify {
    label 'singlecore'
    publishDir "$output_famvcf", mode: 'copy', pattern: '*.classify.vcf.gz'

    input:
        file(pruned_vcf)
        file(sex_file)
        file(training_bedpe)
        tuple file(mei_regions), file(mei_index)
    
    output:
        file "${params.cohort_id}.classify.vcf.gz"

    script:

    """
    zcat $pruned_vcf \
    | svtools vcftobedpe  \
    | svtools varlookup -a stdin -b $training_bedpe -c HQ -d 50 \
    | svtools bedpetovcf \
    | svtools vcfsort \
    | bcftools view -i "INFO/HQ_AF > 0" \
    | bgzip -c > training.vars.vcf.gz

    zcat $pruned_vcf \
    | svtools classify \
    -g $sex_file \
    -a $mei_regions \
    -m large_sample \
    | bgzip -c > ${params.cohort_id}.classify.vcf.gz
    """     
}

process quality_filter_sv {
    label 'singlecore'

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
    publishDir "$output_famvcf", mode: 'copy', pattern: '*.filters.vcf.gz'

    input:
        file(temp_vcf)
        file(bed_file)
        tuple file(lowcomplexity_regions), file(lowcomplexity_index)
    
    output:
        file "${params.cohort_id}.classify.filters.vcf.gz"

    script:
    """
    bedtools coverage -wa \
    -a $bed_file \
    -b $lowcomplexity_regions \
    > overlapExcludeRegions.txt
    bgzip $temp_vcf
    
    FilterSVByValue.py \
    -v ${temp_vcf}.gz \
    -b overlapExcludeRegions.txt \
    -c 9 -i 5 -t LowComplexity -x 0.7 -l 1000 -o ${params.cohort_id}.classify.filters.vcf
    sed 's/:\\./:0/g' ${params.cohort_id}.classify.filters.vcf | bgzip -c > ${params.cohort_id}.classify.filters.vcf.gz
    """
}

