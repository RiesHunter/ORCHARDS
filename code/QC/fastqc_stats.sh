## TIME START
SECONDS=0



## Sample average for perbaseN.txt
cd ./perbaseN-txt/
for reads in *_perbaseN.txt
do
f=${reads%%.txt}
printf "file,N_average\n" > perbaseN_compiled.csv
done
for reads in *_perbaseN.txt
do
f=${reads%%.txt}
sed -e '1,2d' < ${f}.txt | awk '{ total += $2; count++ } END { print total/count }' > temp.csv
a=$(cat temp.csv)
line=$(echo $f","$a)
echo $line >> perbaseN_compiled.csv
done






## REPORT TIME
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
