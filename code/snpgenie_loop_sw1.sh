## run VCF-split_by_gene.sh beforehand
## be sure to change refseq in it...
## wd should be ~/GitHub/SNPGenie
## this script only analyzes HA
## this script only saves sliding window

## start time
SECONDS=0

## ECHO REFSEQ OF CHOICE
echo ""
echo ${1}
echo ""


## CALCULATE TOTAL N VCFS
cd /Users/hries/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared\ drives/TCF\ lab/Current\ Lab\ Members/Hunter_Ries/_bioinformatics/SNPGenie/segments/HA
declare -i x=0
for reads in *.vcf
do
x=$(( x + 1 ))
done
declare -i n=0


## loop
for reads in *_recall_fn_ann_HA.vcf
do
f=${reads%%_recall_fn_ann_HA.vcf}
echo ${f}
cd /Users/hries/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared\ drives/TCF\ lab/Current\ Lab\ Members/Hunter_Ries/_bioinformatics/SNPGenie


## genes
echo -ne "[$n/$x]: ${f} [--HA--] \r"
gene=HA
snpgenie.pl \
  --vcfformat=2 \
  --snpreport=./segments/${gene}/${f}_recall_fn_ann_${gene}.vcf \
  --fastafile=./seqs/${1}/${1}-${gene}.fasta \
  --gtffile=./seqs/${1}/${1}-${gene}.gtf \
  --slidingwindow=1
chmod 755 SNPGenie_Results/product_results.txt
mv SNPGenie_Results/sliding_window_length* ./${f}_${gene}_sliding_window.txt
rm SNPGenie_Results/product_results.txt
rm SNPGenie_Results/codon_results.txt
rm -r SNPGenie_Results

cd /Users/hries/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared\ drives/TCF\ lab/Current\ Lab\ Members/Hunter_Ries/_bioinformatics/SNPGenie
## echo
n=$(( n + 1 ))
done

## organize
mkdir output
mv *_sliding_window.txt output
cd output
mkdir sample_files

## compile sliding_window.txt
awk /file/ *_sliding_window.txt > header
head -n 1 header > cols
for reads in *_sliding_window.txt
do
f=${reads%%_sliding_window.txt}
# add each file's data to the header
awk /segments/ ${f}_sliding_window.txt >> cols  
mv ${f}_sliding_window.txt ./sample_files
done
mv cols snpgenie_compiled_sliding_window_1.tsv
rm header



## time end          
echo ""
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo ""
