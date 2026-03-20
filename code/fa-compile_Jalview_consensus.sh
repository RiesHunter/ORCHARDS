touch compiled-All.fa
echo ""

for file in consensus-*.fa
do
f=${file%%.fa}
sed "s/>Consensus.*/>${f}/g" ${f}.fa > ${f}.temp



## 17-19 H1N1
#sed "s|consensus-|A/Michigan/45/2015|g" ${f}.temp > ${f}2.temp
#sed "s/2015/2015|EPI_ISL_336680|/g" ${f}2.temp > ${f}3.temp

## 17-18 H3N2
#sed "s|consensus-|A/Hong_Kong/4801/2014|g" ${f}.temp > ${f}2.temp
#sed "s/2014/2014|EPI_ISL_233740|/g" ${f}2.temp > ${f}3.temp

## 18-19 H3N2
sed "s|consensus-|A/Singapore/INFIMH-16-0019/2016|g" ${f}.temp > ${f}2.temp
sed "s/2016/2016|EPI_ISL_296168|/g" ${f}2.temp > ${f}3.temp



cat ${f}3.temp >> compiled-All.fa
echo $f
done

rm *.temp

#mkdir consensus_segments
#mv consensus-*.fa consensus_segments
