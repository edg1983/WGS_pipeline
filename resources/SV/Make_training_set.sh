input=$1 #usually merged.sv.pruned.vcf.gz, more than 30 samples
output=$2 #prefix for .bedpe.gz output
mean_cn_pl="/well/gel/HICF2/software/svtools/scripts/mean_cn.pl"
zjoin="/well/gel/HICF2/software/svtools/scripts/zjoin.py"

echo "Calculate means and percentile"
if [ ! -f per_site.means.txt && ! -f overall_percentiles.txt ]
then 
zcat $input \
| vawk --header '{if((I$SVTYPE=="DEL" || I$SVTYPE=="DUP" || I$SVTYPE=="MEI") && I$AF>0.01 && $1!="X" && $1!="Y") print $0}' \
| perl $mean_cn_pl 1>per_site.means.txt 2>overall_percentiles.txt
else
	echo "per_site.means.txt and overall_percentiles.txt already exist"
fi

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
