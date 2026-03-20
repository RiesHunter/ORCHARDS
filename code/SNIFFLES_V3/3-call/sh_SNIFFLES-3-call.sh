#!/bin/bash

##################################################
### PARAMETERS
echo $state
ls -lh

## Refseq
echo "${ref}"
bwa index ${ref}
samtools faidx ${ref}

## callvariants.sh
minPHRED_call=30
  # ignore variants with lower Phred score
mincov=100
  # ignore variants with lower coverage
minAF=0.0
  # ignore variants with lower AF

##################################################
### ALIGN
## align merged reads to refseq
for reads in *_normalized.fastq.gz
do
f=${reads%%_normalized.fastq.gz}
bwa mem -t 8 ${ref} ./${f}_normalized.fastq.gz > ./${f}_aligned.bam
done
echo "---------- Aligned ----------"

ls -lh

##################################################
### SORT
## sort BAM by coordinate
for reads in *_aligned.bam
do
f=${reads%%_aligned.bam}
samtools sort -o ./${f}_sorted.bam ./${f}_aligned.bam
# samtools index
done
echo "---------- Sorted ----------"

ls -lh

##################################################
### MERGE
samtools merge -o ./${state}_merged.bam ./${state}_sorted.bam ./${state}_R1_sorted.bam ./${state}_R2_sorted.bam
echo "---------- Merged ----------"
ls -lh

##################################################
### CALL——initial
for reads in ./${state}_merged.bam
do
f=${reads%%_merged.bam}
callvariants.sh \
    in=./${f}_merged.bam \
    ref=${ref} \
    minallelefraction=${minAF} \
    rarity=${minAF} \
    coverage=t \
    calldel=f \
    callins=f \
    callsub=t \
    mincov=${mincov} \
    minreads=2 \
    minvarcopies=2 \
    minscore=${minPHRED_call} \
    out=${f}_call.vcf

#lofreq call \
#    -f ${ref} \
#    -q ${minPHRED_call} \
#    -Q ${minPHRED_call} \
#    -C ${mincov} \
#    --force-overwrite \
#    --verbose \
#    -o ./${f}_call_lofreq.vcf \
#    ./${f}_sorted.bam
#bcftools filter -o ./${f}_call.vcf -e 'INFO/AF < 0.01' ./${f}_call.vcf
done
echo "---------- Called ----------"

ls -lh

##################################################
### organize for export
rm *.fasta.*
rm *sorted.bam
rm *aligned.bam
rm *normalized*
rm docker_stderror

##################################################
echo "State:"
echo $state
echo ""

echo "File name:"
echo $f
echo ""

echo "Reference:"
echo $ref
echo ""

echo "Files moved:"
ls *sorted*; ls *call*
