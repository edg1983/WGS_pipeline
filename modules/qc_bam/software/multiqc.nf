params.outdir = 'reports'

process MULTIQC {
    label 'singlecore'
    publishDir "$params.outdir", mode: 'copy'

    input:
    file('*')

    output:
    file 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}