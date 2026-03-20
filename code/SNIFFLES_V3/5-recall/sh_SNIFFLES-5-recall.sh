#!/bin/bash

## This is the recall script, but it's the exact same as the call script--other than this!

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
    coverage=t \
    calldel=f \
    callins=f \
    callsub=t \
    mincov=${mincov} \
    minreads=2 \
    minvarcopies=2 \
    rarity=${minAF} \
    minallelefraction=${minAF} \
    minscore=${minPHRED_call} \
    out=${f}_recall.vcf

#lofreq call \
#    -f ${ref} \
#    -q ${minPHRED_call} \
#    -Q ${minPHRED_call} \
#    -C ${mincov} \
#    --force-overwrite \
#    --call-indels \
#    --verbose \
#    -o ./${f}_call.vcf \
#    ./${f}_sorted.bam
#bcftools filter -o ./${f}_call.vcf -e 'INFO/AF < 0.01' ./${f}_call.vcf
done
echo "---------- Called ----------"

ls -lh

##################################################
### organize for export
rm *.fasta.*
rm *aligned*
rm *_sorted.bam
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
