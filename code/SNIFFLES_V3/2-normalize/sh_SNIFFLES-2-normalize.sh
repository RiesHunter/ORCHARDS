#!/bin/bash

##################################################
### PARAMETERS
echo $state
kmersize=31
downsample=2000

##################################################
### NORMALIZE
for reads in *_filtered.fastq.gz
do
f=${reads%%_filtered.fastq.gz}
bbnorm.sh \
    in=./${f}_filtered.fastq.gz \
    out=./${f}_normalized.fastq.gz \
    target=${downsample} \
    k=${kmersize} \
    maxdepth=200 \
    threads=auto \
    -Xmx75g
done

echo "---------- Normalized ----------"

##################################################
rm docker_stderror
#rm *.fastq.gz

echo "State:"
echo $state
echo ""

echo "f:";echo $f;echo ""

echo "Directory contents:"
ls -lh *
echo ""
