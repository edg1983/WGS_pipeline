# WGS analysis pipeline
WGS analysis pipeline as used in the HICF2 project. Can handle both WGS and WES data.

The whole pipeline use singularity images and will pull images from singularity libraries when needed

## How to run
The pipeline can be run directly using Nextflow >= v20.10.

Eventually update `singularity_cachedir` variable in `nextflow.config` to point to a proper folder where singularity images are stored / will be downloaded

A bash script `Run_pipeline.sh` is provided for convenience to set environmental variables for singularity cache and temp directories. The `singularity_dir` variable in the script must match `singularity_cachedir` variable in `nextflow.config`.

### Using Nextflow directly

```
nextflow WGS_analysis.nf -profile cluster --input input_file.txt --mode WGS --ped ped_file.ped --cohort_id cohort_name --outdir results
```

### Using the supporting script

```
Run_pipeline.sh -profile cluster --input input_file.txt --mode WGS --ped ped_file.ped --cohort_id cohort_name --outdir results
```

### Arguments

```
mode        : WGS or WES
input file  : tab-separated file describing input FASTQs
ped file    : standard PED file containing all samples
cohort id   : a arbitrary name for the cohort files generated
outdir      : output folder for results
```

## Input files format
### input_file.txt
A 3 columns tab-separated file without header

```
sampleID1   s1_lane1_R1.fastq.gz    s1_lane1_R2.fastq.gz
sampleID1   s1_lane2_R1.fastq.gz    s1_lane2_R2.fastq.gz
sampleID2   s2_lane2_R1.fastq.gz    s2_lane2_R2.fastq.gz
```

Note that if a sample has been sequenced with multiple pairs of fastq files you need to add multiple lines for each pair of fastq files using the same sampleID. The pipeline will take care of the merge.

### PED file
A standard tab-separated PED file without header, describing all samples provided in the input file. All sample IDs must match between ped and input file. All samples must have sex defined.

```
family_ID   individual_ID   father_ID   mother_ID   sex(1=M,2=F)    status(1=unaff,2=aff,0=unknown)
```

## Pipeline components
1. Alignement and duplicate marking
    - BWA + samblastersamtools
2. QC and coverage from BAM files
    - fastqc + mosdepth + samtools flagstat / mapstat
    - somalier (ancestry, relatedness, sex check reports) 
    - multiqc report generated
3. small variants
    - deepvariant
    - glnexus (gvcf merge)
4. structural variants
    - lumpy + CNVnator
    - svtools (merge and classify)
5. repeat expansion detection
    - expansion hunter
6. ROH regions
    - bcftools ROH