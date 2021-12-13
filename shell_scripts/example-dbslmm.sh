#DIR=~/tmp/DBSLMM
DIR=~/research/DBSLMM
dbslmm=${DIR}/scr/dbslmm
#dbslmm=${DIR}/scr/dbslmm
### Parameters for DBSLMM
let chr=1
DBSLMM=${DIR}/software/DBSLMM.R
summf=${DIR}/test_dat/summary_gemma_chr
outPath=${DIR}/test_dat/out/
plink=/usr/cluster/bin/plink-1.9
ref=${DIR}/test_dat/ref_chr

blockf=${DIR}/block_data/EUR/chr

m=`cat ${summf}${chr}.assoc.txt | wc -l` 
h2=0.5
nobs=`sed -n "2p" ${summf}${chr}.assoc.txt | awk '{print $5}'`
nmis=`sed -n "2p" ${summf}${chr}.assoc.txt | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)
## execute Rscript
Rscript ${DBSLMM} --summary ${summf}${chr}.assoc.txt --outPath ${outPath} \
  --plink ${plink} --dbslmm ${dbslmm} --ref ${ref}${chr} --n ${n} \
  --nsnp ${m} --type auto --model DBSLMM --block ${blockf}${chr}.bed \
  --h2 ${h2} \
  --training true 

