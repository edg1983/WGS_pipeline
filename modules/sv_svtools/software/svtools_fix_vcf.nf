params.outdir = 'SV'
output_dir = file(params.outdir)
output_dir.mkdirs()

process FIX_SVTOOLS_VCF {
    label 'singlecore'
    publishDir "$params.outdir", mode: 'copy'

    input:
        file(input_vcf)
    
    output:
        file "*.raw.vcf.gz"

    script:
    def outprefix = input_vcf.getName().replaceFirst(/\.vcf\.gz/, "")
    """
    Fix_svtools_VCF.py -v $input_vcf -o ${outprefix}.raw.vcf
    bgzip ${outprefix}.raw.vcf
    """
}