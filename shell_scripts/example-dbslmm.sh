### change file permission for dbslmm
dbslmm=../scr/dbslmm
### Parameters for DBSLMM
let chr=1
DBSLMM=../software/DBSLMM.R
summf=../test_dat/summary_gemma_chr
outPath=../test_dat/out/
plink=/usr/cluster/bin/plink-1.9
ref=../test_dat/ref_chr

blockf=../block_data/EUR/chr
m=`cat ${summf}${chr}.assoc.txt | wc -l` 
h2=0.5
nobs=`sed -n "2p" ${summf}${chr}.assoc.txt | awk '{print $5}'`
nmis=`sed -n "2p" ${summf}${chr}.assoc.txt | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)
## execute Rscript
Rscript ${DBSLMM} --summary ${summf}${chr}.assoc.txt --outPath ${outPath} --plink ${plink} --dbslmm ${dbslmm} --ref ${ref}${chr} --n ${n} --nsnp ${m} --type d --model LMM --block ${blockf}${chr}.bed --h2 ${h2}
### Predict
bfilete=../test_dat/test_chr
est=../test_dat/out/summary_gemma_chr
InterPred=../test_dat/out/internal_pred_chr
## plink 1.9
plink=/usr/cluster/bin/plink-1.9
${plink} --bfile ${bfilete}${chr} --score ${est}${chr}.dbslmm.txt 1 2 4 sum --out ${InterPred}${chr}
## plink 2
plink=plink2
${plink} --bfile ${bfilete}${chr} --score ${est}${chr}.dbslmm.txt 1 2 4 cols=+scoresums --out ${InterPred}${chr}

