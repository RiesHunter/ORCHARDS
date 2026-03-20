#!/bin/bash
# bam_stats.sh
# Description: This script processes .bam files to generate detailed statistics on the 
#   average coverage, coverage percentage, and number of covered sites for specific viral gene segments
# Dependencies: Requires samtools (for 'idxstats', 'flagstat', and 'coverage' commands).
# Usage: Run the script in a directory containing one or more BAM files. Ensure samtools is installed and in your PATH.
# Outputs: Creates three text files with statistical data for the BAM file:
#   1. bam_stats.tsv - A comprehensive summary file where each line corresponds to a segment (e.g., NP, MP, PA) from each BAM file
#   2. position_coverage directory - Contains individual .tsv files for each BAM file with detailed position-based coverage information
#   3. to_remove_missing_segment.txt - Files missing alignment on any segment
#   4. Coverage threshold files - Lists of BAM files that fall below certain average or percentage coverage thresholds across segments
# Notes:
# - The script will overwrite any existing output files with the same name in the current directory.
# - Customize output file names as necessary by modifying the 'outfile_prefix' variable.
# Author: Hunter Ries

## TIME START
SECONDS=0


## CALCULATE TOTAL N BAMS
declare -i x=0
for reads in *.bam
do
x=$(( x + 1 ))
done
declare -i n=0


## Sample stats for average coverage on segments
printf "FileGene\tAvgCov\tSitesg0\tSites\tPerCov\n" > header

for reads in *.bam
do
f=${reads%%.bam}
samtools sort ${f}.bam > ${f}_sort_temp.bam
samtools index ${f}_sort_temp.bam
samtools depth -a -d 0 ${f}_sort_temp.bam > temp
cat temp > ${f}_position_coverage.tsv

# NP
grep A.*NP ./temp > ./temp_gene
fgene=$(printf ${f}""_NP)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[#-------]\r"

# MP
grep A.*MP ./temp > ./temp_gene
fgene=$(printf ${f}""_MP)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[##------]\r"

# PA
grep A.*PA ./temp > ./temp_gene
fgene=$(printf ${f}""_PA)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[###-----]\r"

# HA
grep A.*HA ./temp > ./temp_gene
fgene=$(printf ${f}""_HA)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[####----]\r"

# NA
grep A.*NA ./temp > ./temp_gene
fgene=$(printf ${f}""_NA)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[#####---]\r"

# NS
grep A.*NS ./temp > ./temp_gene
fgene=$(printf ${f}""_NS)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[######--]\r"

# PB1
grep A.*PB1 ./temp > ./temp_gene
fgene=$(printf ${f}""_PB1)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[#######-]\r"

# PB2
grep A.*PB2 ./temp > ./temp_gene
fgene=$(printf ${f}""_PB2)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[########]\r"

rm ${f}_sort_temp.bam; rm *.bai; rm temp; rm temp_gene

n=$(( n + 1 ))
echo "[$n/$x]: $f"
done
mv header bam_stats.tsv


## REPORT TIME             
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "----- Done -----"
echo ""


# Alignments that didn't map at 1 or more segment
grep "\t0\t\t" bam_stats.tsv | cut -f1 -d "_" | sort -u > to_remove_missing_segment

# Alignments with less than 200x coverage at 1 or more segment
awk 'FNR==1 || ($2<300)' bam_stats.tsv | cut -f1 -d "_" | sort -u > to_remove_l300_coverage
awk 'FNR==1 || ($2<275)' bam_stats.tsv | cut -f1 -d "_" | sort -u > to_remove_l275_coverage
awk 'FNR==1 || ($2<250)' bam_stats.tsv | cut -f1 -d "_" | sort -u > to_remove_l250_coverage
awk 'FNR==1 || ($2<225)' bam_stats.tsv | cut -f1 -d "_" | sort -u > to_remove_l225_coverage
awk 'FNR==1 || ($2<200)' bam_stats.tsv | cut -f1 -d "_" | sort -u > to_remove_l200_coverage
awk 'FNR==1 || ($2<175)' bam_stats.tsv | cut -f1 -d "_" | sort -u > to_remove_l175_coverage
awk 'FNR==1 || ($2<150)' bam_stats.tsv | cut -f1 -d "_" | sort -u > to_remove_l150_coverage
awk 'FNR==1 || ($2<125)' bam_stats.tsv | cut -f1 -d "_" | sort -u > to_remove_l125_coverage
awk 'FNR==1 || ($2<100)' bam_stats.tsv | cut -f1 -d "_" | sort -u > to_remove_l100_coverage

# Alignments with coverage of less than 99% of full genome
cat bam_stats.tsv > bam_stats.temp
sed -i '' '/\t0\t\t/d' bam_stats.temp
awk '{sub(/_.*/,"",$1)} 1' bam_stats.temp | awk '{ a[$1]+=$5 }END{ for(i in a) print a[i],i }' | awk '{$1=$1/8} 1' | awk '{ print $2 " " $1 }' > percent_coverage.tsv
awk -F" " '$2<0.99' percent_coverage.tsv | sed 's/\|/ /' | awk '{print $1}' > to_remove_less_than_99percent_coverage
awk -F" " '$2<0.98' percent_coverage.tsv | sed 's/\|/ /' | awk '{print $1}' > to_remove_less_than_98percent_coverage
awk -F" " '$2<0.97' percent_coverage.tsv | sed 's/\|/ /' | awk '{print $1}' > to_remove_less_than_97percent_coverage
awk -F" " '$2<0.96' percent_coverage.tsv | sed 's/\|/ /' | awk '{print $1}' > to_remove_less_than_96percent_coverage
awk -F" " '$2<0.95' percent_coverage.tsv | sed 's/\|/ /' | awk '{print $1}' > to_remove_less_than_95percent_coverage
awk -F" " '$2<0.94' percent_coverage.tsv | sed 's/\|/ /' | awk '{print $1}' > to_remove_less_than_94percent_coverage
awk -F" " '$2<0.93' percent_coverage.tsv | sed 's/\|/ /' | awk '{print $1}' > to_remove_less_than_93percent_coverage
awk -F" " '$2<0.92' percent_coverage.tsv | sed 's/\|/ /' | awk '{print $1}' > to_remove_less_than_92percent_coverage
awk -F" " '$2<0.91' percent_coverage.tsv | sed 's/\|/ /' | awk '{print $1}' > to_remove_less_than_91percent_coverage
awk -F" " '$2<0.90' percent_coverage.tsv | sed 's/\|/ /' | awk '{print $1}' > to_remove_less_than_90percent_coverage

# HA coverage
cat bam_stats.tsv | grep HA > HA_coverage.tsv
awk 'FNR==1 || ($2>200)' HA_coverage.tsv | cut -f1 -d "_" | sort -u > HA_CoverageG200
awk 'FNR==1 || ($2<200)' HA_coverage.tsv | cut -f1 -d "_" | sort -u > HA_CoverageL200
awk 'FNR==1 || ($2>100)' HA_coverage.tsv | cut -f1 -d "_" | sort -u > HA_CoverageG100 
awk 'FNR==1 || ($2<100)' HA_coverage.tsv | cut -f1 -d "_" | sort -u > HA_CoverageL100

# Clean up
sort -u to_remove_missing_segment > to_remove_missing_segment.txt

sort -u to_remove_l300_coverage > to_remove_l300_coverage.tsv
sort -u to_remove_l275_coverage > to_remove_l275_coverage.tsv
sort -u to_remove_l250_coverage > to_remove_l250_coverage.tsv
sort -u to_remove_l225_coverage > to_remove_l225_coverage.tsv
sort -u to_remove_l200_coverage > to_remove_l200_coverage.tsv
sort -u to_remove_l175_coverage > to_remove_l175_coverage.tsv
sort -u to_remove_l150_coverage > to_remove_l150_coverage.tsv
sort -u to_remove_l125_coverage > to_remove_l125_coverage.tsv
sort -u to_remove_l100_coverage > to_remove_l100_coverage.tsv

sort -u to_remove_less_than_99percent_coverage > to_remove_less_than_99percent_coverage.tsv
sort -u to_remove_less_than_98percent_coverage > to_remove_less_than_98percent_coverage.tsv
sort -u to_remove_less_than_97percent_coverage > to_remove_less_than_97percent_coverage.tsv
sort -u to_remove_less_than_96percent_coverage > to_remove_less_than_96percent_coverage.tsv
sort -u to_remove_less_than_95percent_coverage > to_remove_less_than_95percent_coverage.tsv
sort -u to_remove_less_than_94percent_coverage > to_remove_less_than_94percent_coverage.tsv
sort -u to_remove_less_than_93percent_coverage > to_remove_less_than_93percent_coverage.tsv
sort -u to_remove_less_than_92percent_coverage > to_remove_less_than_92percent_coverage.tsv
sort -u to_remove_less_than_91percent_coverage > to_remove_less_than_91percent_coverage.tsv
sort -u to_remove_less_than_90percent_coverage > to_remove_less_than_90percent_coverage.tsv

rm bam_stats.temp to_remove_missing_segment to_remove_l300_coverage to_remove_l275_coverage to_remove_l250_coverage to_remove_l225_coverage to_remove_l200_coverage to_remove_l175_coverage to_remove_l150_coverage to_remove_l125_coverage to_remove_l100_coverage to_remove_less_than_99percent_coverage to_remove_less_than_98percent_coverage to_remove_less_than_97percent_coverage to_remove_less_than_96percent_coverage to_remove_less_than_95percent_coverage to_remove_less_than_94percent_coverage to_remove_less_than_93percent_coverage to_remove_less_than_92percent_coverage to_remove_less_than_91percent_coverage to_remove_less_than_90percent_coverage

mkdir position_coverage/; mv *_position_coverage.tsv position_coverage/
