for reads in *_filtered.fastq.gz
do
f=${reads%%_filtered.fastq.gz}
echo -ne "${f} [-----] \r"
fastqc -q ${f}_filtered.fastq.gz
echo -ne "${f} [#----] \r"
unzip -q ${f}_filtered_fastqc.zip
echo -ne "${f} [##---] \r"
cd ${f}_filtered_fastqc
sed -n '/>>Per base N content/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_perbaseN.txt
sed -n '/>>Per base sequence quality/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_perbaseQ.txt
sed -n '/>>Per sequence quality scores/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_perseqQ.txt
sed -n '/>>Sequence Length Distribution/,/>>END_MODULE/p' ./fastqc_data.txt > ../${f}_seqLengthDist.txt
echo -ne "${f} [###--] \r"
mv Images/per_base_n_content.png ../${f}_perbaseN.png
mv Images/per_base_quality.png ../${f}_perbaseQ.png
mv Images/per_sequence_quality.png ../${f}_perseqQ.png
mv Images/sequence_length_distribution.png ../${f}_seqLengthDist.png
echo -ne "${f} [####-] \r"
cd ../
rm -r ${f}_filtered_fastqc
rm ${f}_filtered_fastqc.zip
rm ${f}_filtered_fastqc.html
echo -ne "${f} [#####] \r"
echo ""
done

mkdir report-filtered
cd report-filtered

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

mkdir seqLengthDist-png
mkdir seqLengthDist-txt
mv ../*seqLengthDist.png ./seqLengthDist-png
mv ../*seqLengthDist.txt ./seqLengthDist-txt

cd ../

echo "----- read statistics completed -----"
echo ""
