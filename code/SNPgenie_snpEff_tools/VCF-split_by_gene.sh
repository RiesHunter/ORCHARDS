## run before snpgenie
# change refseq at top and bottom


## CALCULATE TOTAL N VCFS
declare -i x=0
for reads in *.vcf
do
x=$(( x + 1 ))
done
declare -i n=0



for file in *.vcf
do
f=${file%%.vcf}


## clarify CHROM annotations so we don't get false positives
# 17-19 H1N1
#sed -i.bak 's|A/Michigan/45/2015|CHROM_REF_string|g' ${f}.vcf
#sed -i.bak 's/CHROM_REF_string|EPI_ISL_336680|/foobar/g' ${f}.vcf

# 17-18 H3N2
#sed -i.bak 's|A/Hong_Kong/4801/2014|CHROM_REF_string|g' ${f}.vcf
#sed -i.bak 's/CHROM_REF_string|EPI_ISL_233740|/foobar/g' ${f}.vcf

# 18-19 H3N2
sed -i.bak 's|A/Singapore/INFIMH-16-0019/2016|CHROM_REF_string|g' ${f}.vcf
sed -i.bak 's/CHROM_REF_string|EPI_ISL_296168|/foobar/g' ${f}.vcf


## genes
gene=HA
sed '/#CHROM/q' ${f}.vcf >  ${f}_${gene}.vcf
awk /foobar${gene}/ ${f}.vcf > ${f}_${gene}_iSNVs.vcf
sed -i.bak '/##contig=/d' ${f}_${gene}_iSNVs.vcf
cat ${f}_${gene}_iSNVs.vcf >> ${f}_${gene}.vcf
rm *.bak
rm *iSNVs.vcf

gene=MP
sed '/#CHROM/q' ${f}.vcf >  ${f}_${gene}.vcf   
awk /foobar${gene}/ ${f}.vcf > ${f}_${gene}_iSNVs.vcf 
sed -i.bak '/##contig=/d' ${f}_${gene}_iSNVs.vcf
cat ${f}_${gene}_iSNVs.vcf >> ${f}_${gene}.vcf
rm *.bak
rm *iSNVs.vcf

gene=NA
sed '/#CHROM/q' ${f}.vcf >  ${f}_${gene}.vcf   
awk /foobar${gene}/ ${f}.vcf > ${f}_${gene}_iSNVs.vcf 
sed -i.bak '/##contig=/d' ${f}_${gene}_iSNVs.vcf
cat ${f}_${gene}_iSNVs.vcf >> ${f}_${gene}.vcf
rm *.bak
rm *iSNVs.vcf

gene=NP
sed '/#CHROM/q' ${f}.vcf >  ${f}_${gene}.vcf   
awk /foobar${gene}/ ${f}.vcf > ${f}_${gene}_iSNVs.vcf 
sed -i.bak '/##contig=/d' ${f}_${gene}_iSNVs.vcf
cat ${f}_${gene}_iSNVs.vcf >> ${f}_${gene}.vcf
rm *.bak
rm *iSNVs.vcf

gene=NS
sed '/#CHROM/q' ${f}.vcf >  ${f}_${gene}.vcf   
awk /foobar${gene}/ ${f}.vcf > ${f}_${gene}_iSNVs.vcf 
sed -i.bak '/##contig=/d' ${f}_${gene}_iSNVs.vcf
cat ${f}_${gene}_iSNVs.vcf >> ${f}_${gene}.vcf
rm *.bak
rm *iSNVs.vcf

gene=PA
sed '/#CHROM/q' ${f}.vcf >  ${f}_${gene}.vcf   
awk /foobar${gene}/ ${f}.vcf > ${f}_${gene}_iSNVs.vcf 
sed -i.bak '/##contig=/d' ${f}_${gene}_iSNVs.vcf
cat ${f}_${gene}_iSNVs.vcf >> ${f}_${gene}.vcf
rm *.bak
rm *iSNVs.vcf

gene=PB1
sed '/#CHROM/q' ${f}.vcf >  ${f}_${gene}.vcf   
awk /foobar${gene}/ ${f}.vcf > ${f}_${gene}_iSNVs.vcf 
sed -i.bak '/##contig=/d' ${f}_${gene}_iSNVs.vcf
cat ${f}_${gene}_iSNVs.vcf >> ${f}_${gene}.vcf
rm *.bak
rm *iSNVs.vcf

gene=PB2
sed '/#CHROM/q' ${f}.vcf >  ${f}_${gene}.vcf
awk /foobar${gene}/ ${f}.vcf > ${f}_${gene}_iSNVs.vcf
sed -i.bak '/##contig=/d' ${f}_${gene}_iSNVs.vcf
cat ${f}_${gene}_iSNVs.vcf >> ${f}_${gene}.vcf
rm *.bak
rm *iSNVs.vcf


## change foobar${gene} back to full annotation
## 17-19 H1N1
#sed -i.bak 's|foobar|A/Michigan/45/2015foobar|g' ${f}*.vcf
#sed -i.bak 's/foobar/|EPI_ISL_336680|/g' ${f}*.vcf

## 17-18 H3N2
#sed -i.bak 's|foobar|A/Hong_Kong/4801/2014foobar|g' ${f}*.vcf
#sed -i.bak 's/foobar/|EPI_ISL_233740|/g' ${f}*.vcf

## 18-19 H3N2
sed -i.bak 's|foobar|A/Singapore/INFIMH-16-0019/2016foobar|g' ${f}*.vcf
sed -i.bak 's/foobar/|EPI_ISL_296168|/g' ${f}*.vcf


rm *.bak


## echo
n=$(( n + 1 ))
echo "[$n/$x]: $f"

done


## organize
mkdir segments
mkdir segments/HA; mv *HA.vcf segments/HA
mkdir segments/MP; mv *MP.vcf segments/MP
mkdir segments/NA; mv *NA.vcf segments/NA
mkdir segments/NP; mv *NP.vcf segments/NP
mkdir segments/NS; mv *NS.vcf segments/NS
mkdir segments/PA; mv *PA.vcf segments/PA
mkdir segments/PB1; mv *PB1.vcf segments/PB1
mkdir segments/PB2; mv *PB2.vcf segments/PB2

echo "Did you remember to change the refseq in this script?"
