for file in *.fa
do
f=${file%%.fa}

## NP
sed -n '/>*NP/p' ./${f}.fa > name
sed -i '' -e "s/NP/NP|${f}/g" name
cat name >> segmented_compiled-NP
sed -n '/>*NP/,/>/p' ./${f}.fa > seq
sed -i '' -e '$ d' seq
sed -i '' -e '1d' seq
cat seq >> segmented_compiled-NP
echo -ne "[#--------] ${f}\r"

## NS
sed -n '/>*NS/p' ./${f}.fa > name
sed -i '' -e "s/NS/NS|${f}/g" name
cat name >> segmented_compiled-NS
sed -n '/>*NS/,/>/p' ./${f}.fa > seq
sed -i '' -e '$ d' seq
sed -i '' -e '1d' seq
cat seq >> segmented_compiled-NS
echo -ne "[##-------] ${f}\r"

## MP
sed -n '/>*MP/p' ./${f}.fa > name
sed -i '' -e "s/MP/MP|${f}/g" name
cat name >> segmented_compiled-MP
sed -n '/>*MP/,/>/p' ./${f}.fa > seq
sed -i '' -e '$ d' seq
sed -i '' -e '1d' seq
cat seq >> segmented_compiled-MP
echo -ne "[###------] ${f}\r"

## PA
sed -n '/>*PA/p' ./${f}.fa > name
sed -i '' -e "s/PA/PA|${f}/g" name
cat name >> segmented_compiled-PA
sed -n '/>*PA/,/>/p' ./${f}.fa > seq
sed -i '' -e '$ d' seq
sed -i '' -e '1d' seq
cat seq >> segmented_compiled-PA
echo -ne "[####-----] ${f}\r"

## PB1
sed -n '/>*PB1/p' ./${f}.fa > name
sed -i '' -e "s/PB1/PB1|${f}/g" name
cat name >> segmented_compiled-PB1
sed -n '/>*PB1/,/>/p' ./${f}.fa > seq
sed -i '' -e '$ d' seq
sed -i '' -e '1d' seq
cat seq >> segmented_compiled-PB1
echo -ne "[#####----] ${f}\r"

## PB2
sed -n '/>*PB2/p' ./${f}.fa > name
sed -i '' -e "s/PB2/PB2|${f}/g" name
cat name >> segmented_compiled-PB2
sed -n '/>*PB2/,/>/p' ./${f}.fa > seq
sed -i '' -e '$ d' seq
sed -i '' -e '1d' seq
cat seq >> segmented_compiled-PB2
echo -ne "[######---] ${f}\r"

## NA
sed -n '/>*NA/p' ./${f}.fa > name
sed -i '' -e "s/NA/NA|${f}/g" name
cat name >> segmented_compiled-NA
sed -n '/>*NA/,/>/p' ./${f}.fa > seq
sed -i '' -e '$ d' seq
sed -i '' -e '1d' seq
cat seq >> segmented_compiled-NA
echo -ne "[#######--] ${f}\r"

## HA
sed -n '/>*HA/p' ./${f}.fa > name
sed -i '' -e "s/HA/HA|${f}/g" name
cat name >> segmented_compiled-HA
sed -n '/>*HA/,/>/p' ./${f}.fa > seq
#sed -i '' -e '$ d' seq
sed -i '' -e '1d' seq
cat seq >> segmented_compiled-HA
echo -ne "[########-] ${f}\r"

## Whole Genome
echo ">Whole_Genome|${f}" > name
sed '/>/d' ./${f}.fa > seq
cat name >> segmented_compiled-Whole_Genome
cat seq >> segmented_compiled-Whole_Genome
echo -ne "[#########] ${f}\r"

echo ""
done


mv segmented_compiled-NP segmented_compiled-NP.fa
mv segmented_compiled-NS segmented_compiled-NS.fa
mv segmented_compiled-MP segmented_compiled-MP.fa
mv segmented_compiled-PA segmented_compiled-PA.fa
mv segmented_compiled-PB2 segmented_compiled-PB2.fa
mv segmented_compiled-PB1 segmented_compiled-PB1.fa
mv segmented_compiled-NA segmented_compiled-NA.fa
mv segmented_compiled-HA segmented_compiled-HA.fa
mv segmented_compiled-Whole_Genome segmented_compiled-Whole_Genome.fa

mkdir segmented_compiled_fastas
mv segmented_compiled-* segmented_compiled_fastas

rm name
rm seq
