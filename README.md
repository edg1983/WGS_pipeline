# WGS analysis pipeline
WGS analysis pipeline. Can handle both WGS and WES data.

The whole pipeline use singularity images and will pull images from singularity library when needed. Singularity recipes used are provided in `singularity` folder for reference.

## How to run
The pipeline can be run directly using Nextflow >= v20.10.

```
nextflow WGS_analysis.nf -profile cluster --operation align --input input_file.txt --mode WGS --ped ped_file.ped --ref genome.fa --cohort_id cohort_name --outdir results 
```

The pipeline automatically infer the number of samples in the cohort from your input file and adjust the filtering accordingly. When more than one sample is present, small variants and structural variants from all samples are merged in cohort wide VCF files. 

Eventually update `singularity_cachedir` variable in `nextflow.config` to point to a proper folder where singularity images are stored / will be downloaded

**NB.** Various supporting files are needed and expected in the `resources` folder. This path can be configured by changing the parameters in `config/resources_GRCh37/38.conf`. All files needed are provided in a Zenodo repository. Please refer to the README file in the resources folder. 

### Arguments

```
operation   :   align or call_variants
mode        :   WGS only supported at the moment
ref         :   fasta file for the genome. Note that .fai and bwa index are expected in the same location
input       :   tab-separated file describing input files. 
                The exact format depends on operation requested (see below)
ped         :   standard PED file containing all samples
cohort_id   :   a arbitrary name for the cohort files generated
outdir      :   output folder for results
```

Use `--operation align/call_variants --help` for more explanations. 

## Input files format
### PED file
A standard tab-separated PED file without header, describing all samples provided in the input file. All sample IDs must match between ped and input file. All samples must have sex defined.

```
family_ID   individual_ID   father_ID   mother_ID   sex(1=M,2=F)    status(1=unaff,2=aff,0=unknown)
```

### input file

Note that all files need to be specified using **absolute paths**

#### Operation: align
A 3 columns tab-separated file without header

```
sampleID1   s1_lane1_R1.fastq.gz    s1_lane1_R2.fastq.gz
sampleID1   s1_lane2_R1.fastq.gz    s1_lane2_R2.fastq.gz
sampleID2   s2_lane2_R1.fastq.gz    s2_lane2_R2.fastq.gz
```

Note that if a sample has been sequenced with multiple pairs of fastq files you need to add multiple lines for each pair of fastq files using the same sampleID. The pipeline will take care of the merge.

#### Operation: call_variants
A 5 columns tab-separated file without header.
This file is automatically generated in the output folder when using `--operation align` (bam_files.txt) 

```
sampleID1   main_bam.bam    disc.bam    split.bam
sampleID2   main_bam.bam    disc.bam    split.bam
sampleID3   main_bam.bam    disc.bam    split.bam
```

`disc` and `split` BAM files are files containing only discordant pair and split reads like the 
ones that can be obtained using Samblaster

## Output

The pipeline generates a reach set of outputs including
- aligned deduplicated BAM files
- disc/split BAM files
- Extensive QC of alignements, which includes mapping stats, coverage, relatedness, ancestry
- Multi sample and single sample VCFs of small variants and structural variants (variants are provided as raw calls and filtered calls)
- Variants QC report for small variants
- ROH regions
- Repeat expansions by Expansion Hunter

## Pipeline components
1. Alignement and duplicate marking
    - BWA-MEM + samblaster + samtools
2. QC and coverage from BAM files
    - fastqc: reads stats
    - mosdepth: coverage
    - samtools flagstat / mapstat: alignment stats
    - somalier: ancestry, relatedness, sex check reports 
    - multiqc: interactive report
3. small variants
    - deepvariant: single sample calls
    - glnexus: gvcf merge 
4. structural variants
    - lumpy: structural variants events
    - CNVnator: CNV estimation
    - svtools: combine, merge and classify
5. repeat expansion detection
    - expansion hunter
6. ROH regions
    - bcftools ROH

## Future developments

- [ ] Support for WES
- [ ] Update SV pipeline to Manta / dysgu 
- [ ] Add duphold for SV quality check
- [ ] Variant annotation
- [ ] Segregation analysis