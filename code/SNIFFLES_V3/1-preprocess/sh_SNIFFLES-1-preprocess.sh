#!/bin/bash

##################################################
### Define parameters
## Trimmomatic
window=5
minPHRED=20
  # both are for sliding window q control

adapter=Nextera_XT_adapter.fa
echo ">Nextera_XT_adapter" > Nextera_XT_adapter.fa
echo "CTGTCTCTTATACACATCT" >> Nextera_XT_adapter.fa

## Reformat.sh
minlen=100
  # Reads shorter than this after trimming will be discarded
  # Pairs will be discarded only if both are shorter

##################################################
### TRIM
for reads in *_R1.fastq.gz
do
f=${reads%%_R1.fastq.gz}
java -jar /usr/bin/Trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE \
    ./${f}_R1.fastq.gz \
    ./${f}_R2.fastq.gz \
    ./${f}_trimmed_R1.fastq.gz \
    ./${f}_untrimmed_R1.fastq.gz \
    ./${f}_trimmed_R2.fastq.gz \
    ./${f}_untrimmed_R2.fastq.gz \
    ILLUMINACLIP:./${adapter}:2:30:10 \
    SLIDINGWINDOW:${window}:${minPHRED} \
    MINLEN:100
done
echo "---------- Trimmed ----------"

##################################################
### MERGE
for reads in *_trimmed_R1.fastq.gz
do
f=${reads%%_trimmed_*.fastq.gz}
f=${f##*/}
bbmerge.sh \
    in1=./${f}_trimmed_R1.fastq.gz \
    in2=./${f}_trimmed_R2.fastq.gz \
    out=./${f}_merged.fastq.gz \
    outu=./${f}_unmerged_R1.fastq.gz \
    outu2=./${f}_unmerged_R2.fastq.gz
done
echo "---------- Merged ----------"

##################################################
### FILTER
for reads in *_merged.fastq.gz
do
f=${reads%%_merged*.fastq.gz}
f=${f##*/}
reformat.sh \
    qtrim=t \
    minlength=${minlen} \
    in=./${f}_merged.fastq.gz \
    out=./${f}_filtered.fastq.gz
done
echo "---------- Filtered ----------"

##################################################
### FILTER (unmerged)
for reads in *_unmerged_R1.fastq.gz
do
f=${reads%%_unmerged_R1.fastq.gz}
f=${f##*/}
reformat.sh \
    qtrim=t \
    minlength=${minlen} \
    in=./${f}_unmerged_R1.fastq.gz \
    out=./${f}_R1_filtered.fastq.gz
done

for reads in *_unmerged_R2.fastq.gz
do
f=${reads%%_unmerged_R2.fastq.gz}
f=${f##*/}
reformat.sh \
    qtrim=t \
    minlength=${minlen} \
    in=./${f}_unmerged_R2.fastq.gz \
    out=./${f}_R2_filtered.fastq.gz
done
echo "---------- Filtered ----------"

##################################################
### Clean
## remove unneeded files
du -h
ls -lah
rm docker_stderror
rm *merged.fastq.gz
rm *trimmed*
rm *untrimmed*
rm *R1.fastq.gz
rm *R2.fastq.gz
rm *Nextera*
echo "---------- Unneeded files removed ----------"
du -h
ls -lah

##################################################
echo "f:";echo $f;echo ""

echo "Directory contents:"
ls -lh *
echo ""


