#!/bin/bash

##################################################
### PARAMETERS
echo $state
ls -lh

## Refseq
#ref=*.fasta
echo "${ref}"
cat ${ref}
bwa index ${ref}
samtools faidx ${ref}

##################################################
### Remove SNVs below 50%
for reads in *_call.vcf
do
f=${reads%%_call.vcf}
f=${f##*/}
bcftools filter -o ./${f}_call_consensus.vcf -e 'INFO/AF < 0.5' ./${f}_call.vcf
done
cp ./${f}_call_consensus.vcf ./${f}_call_consensus-1.vcf
ls -lh
echo "---------- Consensus variants ----------"

##################################################
### CONSENSUS
## develop per-sample consensus from VCF
for reads in *_call_consensus.vcf
do
f=${reads%%_call_consensus.vcf}
f=${f##*/}
bgzip ./${f}_call_consensus.vcf
bcftools index ./${f}_call_consensus.vcf.gz
cat ${ref} | \
    bcftools consensus ./${f}_call_consensus.vcf.gz \
    -o ./${f}_consensus.fa
done
ls -lh
echo "---------- Consensus developed ----------"

##################################################
### organize for export
rm *.fasta.*
rm *.vcf.gz
rm *.vcf.gz.csi
rm *call.vcf
rm ./${f}_call_consensus-1.vcf

##################################################
echo "State:"
echo $state
echo ""
rm docker_stderror

echo "Process:"
echo $Process
echo ""

echo "f:";echo $f;echo ""

echo "Files moved:"
ls -lah *consensus.fa
