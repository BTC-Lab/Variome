# ERGP_Variome
Tools, Scripts, and cmd-line used in the ERGP Variome project


## to account for datalake extraction issues ##
# fix compression and replace header with own, working version
ls *.vcf.gz | grep -v / | sed 's/\./ /' | awk '{print "(cat new_header.txt && zcat "$1".vcf.gz | tail -n +41) | bgzip -c > "$1"_comp.vcf.gz"}' > 00_run_recompression.sh

# sort the files
ls *comp.vcf.gz | grep -v / | sed 's/\./ /' | awk '{print "bcftools sort "$1".vcf.gz --output-type z -o "$1"_sorted.vcf.gz"}' > 00_run_sort.sh

# and fix the weird concatenation issue
ls *sorted.vcf.gz | grep -v / | sed 's/\./ /' | awk '{print "bcftools norm -m +any "$1".vcf.gz --output-type z -o "$1"_normed.vcf.gz"}' > 00_run_norm.sh
