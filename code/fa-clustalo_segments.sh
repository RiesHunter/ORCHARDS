mv segmented_compiled-Whole_Genome.fa compiled-Whole_Genome.fa

SECONDS=0

## NP
clustalo -i segmented_compiled-NP.fa -o clustalo-NP.fasta -v

## NS
clustalo -i segmented_compiled-NS.fa -o clustalo-NS.fasta -v

## MP
clustalo -i segmented_compiled-MP.fa -o clustalo-MP.fasta -v

## PA
clustalo -i segmented_compiled-PA.fa -o clustalo-PA.fasta -v

## PB1
clustalo -i segmented_compiled-PB1.fa -o clustalo-PB1.fasta -v

## PB2
clustalo -i segmented_compiled-PB2.fa -o clustalo-PB2.fasta -v

## NA
clustalo -i segmented_compiled-NA.fa -o clustalo-NA.fasta -v

## HA
clustalo -i segmented_compiled-HA.fa -o clustalo-HA.fasta -v
echo ""

touch consensus-HA.fa
touch consensus-MP.fa
touch consensus-NA.fa
touch consensus-NP.fa
touch consensus-NS.fa
touch consensus-PA.fa
touch consensus-PB1.fa
touch consensus-PB2.fa

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
