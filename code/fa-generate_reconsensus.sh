### Re-consensus generation from VCF
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



cp ./18-19_H3N2_within-season_consensus.fasta ./segments/HA
cd ./segments/HA

##################################################
### Refseq
ref=*.fasta
echo "${ref}"
bwa index ${ref}
samtools faidx ${ref}

##################################################
### Remove SNVs below 50%
for reads in *_recall_fn_ann_HA.vcf
do
f=${reads%%_recall_fn_ann_HA.vcf}
f=${f##*/}
bcftools filter -o ./${f}_recall_reconsensus.vcf -e 'INFO/AF < 0.5' ./${f}_recall_fn_ann_HA.vcf
done
echo "---------- Consensus variants ----------"

##################################################
### CONSENSUS
## develop per-sample consensus from VCF
for reads in *_recall_reconsensus.vcf
do
f=${reads%%_recall_reconsensus.vcf}
f=${f##*/}
bgzip ./${f}_recall_reconsensus.vcf
bcftools index ./${f}_recall_reconsensus.vcf.gz
cat ${ref} | \
    bcftools consensus ./${f}_recall_reconsensus.vcf.gz \
    -o ./${f}_reconsensus.fa
done
echo ""; echo "---------- Consensus developed ----------"; echo ""

##################################################
### organize for export
rm *.fasta.*
rm *.vcf.gz
rm *.vcf.gz.csi

##################################################
## Concatenate all HA .fa files
for file in *.fa
do
f=${file%%.fa}
## HA
sed -n '/>*HA/p' ./${f}.fa > name
sed "s/HA/HA|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-HA
sed -n '/>*HA/,/>/p' ./${f}.fa > seq
sed '$ d' seq > seq.temp
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-HA
echo -ne "${f}\r"
echo ""
done

mv segmented_compiled-HA segmented_compiled-HA.fasta

rm seq
rm name
rm name.temp
rm seq.temp

### remove everything before HA
sed -i.bak 's|>A/Singapore/INFIMH-16-0019/2016|>|g' segmented_compiled-HA.fasta
sed -i.bak 's|EPI_ISL_296168||g' segmented_compiled-HA.fasta
sed -i.bak 's|HA||g' segmented_compiled-HA.fasta
sed -i.bak 's|_reconsensus||g' segmented_compiled-HA.fasta
sed -i.bak 's/|||//g' segmented_compiled-HA.fasta

### Align concatenated HA .fa file
#clustalw \
#  -ALIGN \
#  -INFILE=./segmented_compiled-HA.fasta \
#  -OUTFILE=./clustalw_HA.fasta \
#  -OUTPUT=FASTA
#
#### raxml-ng ML
## Create .phy
#raxml-ng --check --msa ./clustalw_*.fasta --model GTR+G
## Check for MSA
#raxml-ng --parse --msa ./clustalw_*.fasta --model GTR+G
## Infer the tree
#raxml-ng --msa ./clustalw_*.phy --model GTR+G --prefix clustalw --threads 2 --seed 920

## One-liner
#raxml-ng --check --msa ./clustalw_*.fasta --model GTR+G; raxml-ng --parse --msa ./clustalw_*.phy --model GTR+G; raxml-ng --msa ./clustalw_*.phy --model GTR+G --prefix clustalw --threads 2 --seed 920