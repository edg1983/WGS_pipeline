input=$1 #usually merged.sv.pruned.vcf.gz, more than 30 samples
output=$2 #prefix for .bedpe.gz output
mean_cn_pl="/well/gel/HICF2/software/svtools/scripts/mean_cn.pl"
zjoin="/well/gel/HICF2/software/svtools/scripts/zjoin.py"

echo "Calculate means and percentile" 
zcat $input \
| vawk --header '{if((I$SVTYPE=="DEL" || I$SVTYPE=="DUP" || I$SVTYPE=="MEI") && I$AF>0.01 && $1!="X" && $1!="Y") print $0}' \
| perl $mean_cn_pl 1>per_site.means.txt 2>overall_percentiles.txt

echo "Extracting training vars"
cat per_site.means.txt \
| cut -f -8 \
| $zjoin -a stdin -b <(cat overall_percentiles.txt | cut -f -8 ) -1 2 -2 1 \
| awk '{if($5>$11 && $5<$12 && $8>$15 && $8<$16) print $0}' \
| cut -f 1 \
| $zjoin -a <(zgrep -v "#" $input) -b stdin -1 3 -2 1 \
| cat <(zgrep "#" $input) - \
| svtools vcftobedpe \
| bgzip -c > ${output}.bedpe.gz

#Then to use classify on a single/small sample you need to compute AFs (afreq)
#svtools afreq input_file.vcf.gz | svtools vcftobedpe | svtools varlookup -a stdin -b /well/gel/HICF2/software/NF_pipeline/HICF2/resources/SV/HICF2_300samples_GRCh38.trainingVars.bedpe.gz -c HQ -d 50 | svtools bedpetovcf | svtools vcfsort | bcftools view -i "INFO/HQ_AF > 0" | bgzip -c > training.vars.vcf.gz

#Or if VCF already contains AF INFO
#zcat input_file.vcf.gz | svtools vcftobedpe | svtools varlookup -a stdin -b /well/gel/HICF2/software/NF_pipeline/HICF2/resources/SV/HICF2_300samples_GRCh38.trainingVars.bedpe.gz -c HQ -d 50 | svtools bedpetovcf | svtools vcfsort | bcftools view -i "INFO/HQ_AF > 0" | bgzip -c > training.vars.vcf.gz
