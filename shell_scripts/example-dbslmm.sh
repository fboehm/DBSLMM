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
  --h2 ${h2}  
### Predict
bfilete=${DIR}/test_dat/test_chr
est=${DIR}/test_dat/out/summary_gemma_chr
InterPred=${DIR}/test_dat/out/internal_pred_chr
## plink 1.9
plink=/usr/cluster/bin/plink-1.9
${plink} --bfile ${bfilete}${chr} --score ${est}${chr}.dbslmm.txt 1 2 4 sum --out ${InterPred}${chr}
## plink 2
plink=plink2
${plink} --bfile ${bfilete}${chr} --score ${est}${chr}.dbslmm.txt 1 2 4 cols=+scoresums --out ${InterPred}${chr}

