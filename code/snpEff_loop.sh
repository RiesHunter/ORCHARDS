## TIME START
SECONDS=0

## CALCULATE TOTAL N VCFS
declare -i x=0
for reads in *_recall.vcf
do
x=$(( x + 1 ))
done
declare -i n=0

## RENAME VCFS
for reads in *_recall.vcf
do
f=${reads%%_recall.vcf}
sed  "s/PASS/${f}/" ${f}_recall.vcf > ${f}_temp.vcf
sed '/\t\t/d' ${f}_temp.vcf > ${f}_recall_fn.vcf
rm ${f}_temp.vcf
done
echo ""
echo "---------- Renaming complete ----------"

## REPORT TIME
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

## TIME START
SECONDS=0

## ECHO DATABASE OF CHOICE
echo ""
echo $1
echo ""

## SNPEFF WITH DATABASE OF CHOICE
for reads in *_fn.vcf
do
f=${reads%%.vcf}
java -jar snpEff.jar ann $1 ${f}.vcf > ${f}_ann.vcf
rm snpEff_summary.html
mv snpEff_genes.txt ${f}_snpEff_genes.tsv
sed -i '' '/# The following table is formatted as tab separated values./d' ./${f}_snpEff_genes.tsv
sed -i '' 's/#//g' ${f}_snpEff_genes.tsv
n=$(( n + 1 ))
echo "[$n/$x]: $f"
done

## ORGANIZE FILES
mkdir output
cd output
mkdir recall
mkdir recall_fn
mkdir recall_fn_ann
mkdir snpEff_genes
mv ../*_recall.vcf recall
mv ../*_recall_fn.vcf recall_fn
mv ../*_recall_fn_ann.vcf recall_fn_ann
mv ../*_recall_fn_snpEff_genes.tsv snpEff_genes
echo ""
echo "---------- Completed -----------"

## REPORT TIME
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
