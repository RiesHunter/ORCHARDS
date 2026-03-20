mkdir singlets doublets

# for every R1, move the matching R2
for file in *_R1*gz
do
file="${file%%_R1*.fastq.gz}"
file="$file"_R2
mv $file*gz doublets
done

cd doublets

# for every R2, move the matching R1 in higher dir
for file in *_R2*gz
do
file="${file%%_R2*.fastq.gz}"
file="$file"_R1
mv ../$file*gz .
done

