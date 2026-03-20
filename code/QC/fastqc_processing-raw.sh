for reads in *_R1.fastq.gz
do
f=${reads%%_R1.fastq.gz}
echo -ne "${f}_R1 [-----] \r"
fastqc -q ${f}_R1.fastq.gz
echo -ne "${f}_R1 [#----] \r"
unzip -q ${f}*.zip
echo -ne "${f}_R1 [##---] \r"
cd ${f}_R1_fastqc
sed -n '/>>Per base N content/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_R1_perbaseN.txt
sed -n '/>>Per base sequence quality/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_R1_perbaseQ.txt
sed -n '/>>Per sequence quality scores/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_R1_perseqQ.txt
sed -n '/>>Sequence Length Distribution/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_R1_seqLengthDist.txt
echo -ne "${f}_R1 [###--] \r"
mv Images/per_base_n_content.png ../${f}_R1_perbaseN.png
mv Images/per_base_quality.png ../${f}_R1_perbaseQ.png
mv Images/per_sequence_quality.png ../${f}_R1_perseqQ.png
mv Images/sequence_length_distribution.png ../${f}_R1_seqLengthDist.png
echo -ne "${f}_R1 [####-] \r"
cd ../
rm -r ${f}_R1_fastqc
rm ${f}_R1_fastqc.zip
rm ${f}_R1_fastqc.html
echo -ne "${f}_R1 [#####] \r"
echo ""
done

for reads in *_R2.fastq.gz
do
f=${reads%%_R2.fastq.gz}
echo -ne "${f}_R2 [-----] \r"
fastqc -q ${f}_R2.fastq.gz
echo -ne "${f}_R2 [#----] \r"
unzip -q ${f}*.zip
echo -ne "${f}_R2 [##---] \r"
cd ${f}_R2_fastqc
sed -n '/>>Per base N content/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_R2_perbaseN.txt
sed -n '/>>Per base sequence quality/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_R2_perbaseQ.txt
sed -n '/>>Per sequence quality scores/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_R2_perseqQ.txt
sed -n '/>>Sequence Length Distribution/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_R2_seqLengthDist.txt
echo -ne "${f}_R2 [###--] \r"
mv Images/per_base_n_content.png ../${f}_R2_perbaseN.png
mv Images/per_base_quality.png ../${f}_R2_perbaseQ.png
mv Images/per_sequence_quality.png ../${f}_R2_perseqQ.png
mv Images/sequence_length_distribution.png ../${f}_R2_seqLengthDist.png
echo -ne "${f}_R2 [####-] \r"
cd ../
rm -r ${f}_R2_fastqc
rm ${f}_R2_fastqc.zip
rm ${f}_R2_fastqc.html
echo -ne "${f}_R2 [#####] \r"
echo ""
done



mkdir report-raw
cd report-raw

mkdir seqLengthDist-png
mkdir seqLengthDist-txt
mv ../*seqLengthDist.png ./seqLengthDist-png
mv ../*seqLengthDist.txt ./seqLengthDist-txt

mkdir perbaseN-png
mkdir perbaseN-txt
mv ../*perbaseN.png ./perbaseN-png
mv ../*perbaseN.txt ./perbaseN-txt

mkdir perseqQ-png
mkdir perseqQ-txt
mv ../*perseqQ.png ./perseqQ-png
mv ../*perseqQ.txt ./perseqQ-txt

mkdir perbaseQ-png
mkdir perbaseQ-txt
mv ../*perbaseQ.png ./perbaseQ-png
mv ../*perbaseQ.txt ./perbaseQ-txt


cd ../

echo "----- fastqc organization -----"
echo ""
