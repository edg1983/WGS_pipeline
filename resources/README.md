# WGS resource files
These are supporting files needed during various steps of the pipeline
Expected folders are:

- SV_pipeline_files: files used with Lumpy and CNVnator like exlude regions, repeat elements, etc
- chrs: fasta (.fa) files for each chromosome of the genome (used by CNVnator)
- exp_hunter: ExpansionHunter variant catalogs
- geneanno: BED files (.bed.gz) for exons regions (these are used for additional coverage computation)
- ROH_resources: AF files and genetic maps used by bcftools roh detection  
- somalier_data: supporting data for somalier analysis (1000G samples profile, variant catalogs)

The expected locations can be changed in Nextflow configuration (config/resources_GRCh37/38.conf).

Supporting files for GRCh37 and GRCh38 can be downloaded from [Zenodo repository](https://zenodo.org/record/5878875).
Just grab the various archives and decompress them in this folder, this will generate the above folders with necessary data.
