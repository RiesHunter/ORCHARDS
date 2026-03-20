## run VCF-split_by_gene.sh beforehand
## be sure to change refseq in it...
## wd should be ~/GitHub/SNPGenie

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
  --slidingwindow=30
chmod 755 SNPGenie_Results/product_results.txt
mv SNPGenie_Results/product_results.txt ./${f}_${gene}_product_results.txt
mv SNPGenie_Results/codon_results.txt ./${f}_${gene}_codon_results.txt
mv SNPGenie_Results/sliding_window_length* ./${f}_${gene}_sliding_window.txt
rm -r SNPGenie_Results

echo -ne "[$n/$x]: ${f} [--HA-MP--] \r"
gene=MP
snpgenie.pl \
  --vcfformat=2 \
  --snpreport=./segments/${gene}/${f}_recall_fn_ann_${gene}.vcf \
  --fastafile=./seqs/${1}/${1}-${gene}.fasta \
  --gtffile=./seqs/${1}/${1}-${gene}.gtf
chmod 755 SNPGenie_Results/product_results.txt
mv SNPGenie_Results/product_results.txt ./${f}_${gene}_product_results.txt
rm -r SNPGenie_Results

echo -ne "[$n/$x]: ${f} [--HA-MP-NA--] \r"
gene=NA
snpgenie.pl \
  --vcfformat=2 \
  --snpreport=./segments/${gene}/${f}_recall_fn_ann_${gene}.vcf \
  --fastafile=./seqs/${1}/${1}-${gene}.fasta \
  --gtffile=./seqs/${1}/${1}-${gene}.gtf
chmod 755 SNPGenie_Results/product_results.txt
mv SNPGenie_Results/product_results.txt ./${f}_${gene}_product_results.txt
rm -r SNPGenie_Results

echo -ne "[$n/$x]: ${f} [--HA-MP-NA-NP--] \r"
gene=NP
snpgenie.pl \
  --vcfformat=2 \
  --snpreport=./segments/${gene}/${f}_recall_fn_ann_${gene}.vcf \
  --fastafile=./seqs/${1}/${1}-${gene}.fasta \
  --gtffile=./seqs/${1}/${1}-${gene}.gtf
chmod 755 SNPGenie_Results/product_results.txt
mv SNPGenie_Results/product_results.txt ./${f}_${gene}_product_results.txt
rm -r SNPGenie_Results

echo -ne "[$n/$x]: ${f} [--HA-MP-NA-NP-NS--] \r"
gene=NS
snpgenie.pl \
  --vcfformat=2 \
  --snpreport=./segments/${gene}/${f}_recall_fn_ann_${gene}.vcf \
  --fastafile=./seqs/${1}/${1}-${gene}.fasta \
  --gtffile=./seqs/${1}/${1}-${gene}.gtf
chmod 755 SNPGenie_Results/product_results.txt
mv SNPGenie_Results/product_results.txt ./${f}_${gene}_product_results.txt
rm -r SNPGenie_Results

echo -ne "[$n/$x]: ${f} [--HA-MP-NA-NP-NS-PA--] \r"
gene=PA
snpgenie.pl \
  --vcfformat=2 \
  --snpreport=./segments/${gene}/${f}_recall_fn_ann_${gene}.vcf \
  --fastafile=./seqs/${1}/${1}-${gene}.fasta \
  --gtffile=./seqs/${1}/${1}-${gene}.gtf
chmod 755 SNPGenie_Results/product_results.txt
mv SNPGenie_Results/product_results.txt ./${f}_${gene}_product_results.txt
rm -r SNPGenie_Results

echo -ne "[$n/$x]: ${f} [--HA-MP-NA-NP-NS-PA-PB1--] \r"
gene=PB1
snpgenie.pl \
  --vcfformat=2 \
  --snpreport=./segments/${gene}/${f}_recall_fn_ann_${gene}.vcf \
  --fastafile=./seqs/${1}/${1}-${gene}.fasta \
  --gtffile=./seqs/${1}/${1}-${gene}.gtf
chmod 755 SNPGenie_Results/product_results.txt
mv SNPGenie_Results/product_results.txt ./${f}_${gene}_product_results.txt
rm -r SNPGenie_Results

echo -ne "[$n/$x]: ${f} [--HA-MP-NA-NP-NS-PA-PB1-PB2] \r"
gene=PB2
snpgenie.pl \
  --vcfformat=2 \
  --snpreport=./segments/${gene}/${f}_recall_fn_ann_${gene}.vcf \
  --fastafile=./seqs/${1}/${1}-${gene}.fasta \
  --gtffile=./seqs/${1}/${1}-${gene}.gtf
chmod 755 SNPGenie_Results/product_results.txt
mv SNPGenie_Results/product_results.txt ./${f}_${gene}_product_results.txt
rm -r SNPGenie_Results

cd /Users/hries/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared\ drives/TCF\ lab/Current\ Lab\ Members/Hunter_Ries/_bioinformatics/SNPGenie
## echo
n=$(( n + 1 ))
done


## organize
mkdir output
mv *_product_results.txt output
mv *_codon_results.txt output
mv *_sliding_window.txt output
cd output
mkdir sample_files



## compile product_results.txt
awk /file/ ${f}*HA_product_results.txt > header
for reads in *_product_results.txt
do
f=${reads%%_product_results.txt}
# add each file's data to the header
awk /segments/ ${f}_product_results.txt >> header
mv ${f}_product_results.txt ./sample_files
done
sed 's/NA/Neuraminidase/g' header > snpgenie_compiled.tsv
rm header

## compile codon_results.txt
awk /file/ *HA_codon_results.txt > header
head -n 1 header > cols
for reads in *_codon_results.txt
do
f=${reads%%_codon_results.txt}
# add each file's data to the header
awk /segments/ ${f}_codon_results.txt >> cols
mv ${f}_codon_results.txt ./sample_files
done
mv cols snpgenie_compiled_sw.tsv
rm header

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
mv cols snpgenie_compiled_sliding_window.tsv
rm header



## time end          
echo ""
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo ""
