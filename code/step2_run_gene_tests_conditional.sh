#!/usr/bin/bash
# 03-July-2023 Diogo Ribeiro
# Run REGENIE v3.2 step2, gene-based

bgenPrefix=$1 # e.g. all_chr.JAVIERRE_rare.bgen
predFile=$2 # e.g. ukb_step1_pred.list
annoFile=$3 # e.g. c22.JAVIERRE.CADD15
trait=$4 # e.g. p30120_i0
maskFile=$5 #e.g. javi.mask
condListFile=$6 #e.g. all_kanai2021_conditional.sites.txt
condFile=$7 #e.g. all_chr.conditional.pip0.8 
outFile=$8 #e.g. regenie_javierre
outFolder=$9 # e.g. :/regenie/step2/JAVIERRE
# Note: hardcoded path of phenotype and covariates

inputText="-iin=regenie/200k/regenie_phenotype_21July2023.txt -iin=regenie/200k/regenie_covariate_21July2023.txt -iin=$bgenPrefix.pgen -iin=$bgenPrefix.psam -iin=$bgenPrefix.pvar -iin=$predFile -iin=$annoFile.annotations.gz -iin=$annoFile.setlist.gz -iin=$maskFile -iin=$condListFile -iin=$condFile.pgen -iin=$condFile.psam -iin=$condFile.pvar"

cmdText="regenie \
--step 2   \
--pgen $(basename $bgenPrefix) \
--phenoFile regenie_phenotype_21July2023.txt \
--phenoColList $trait \
--covarFile regenie_covariate_21July2023.txt \
--covarColList p31,p21022,p22009_a1,p22009_a2,p22009_a3,p22009_a4,p22009_a5,p22009_a6,p22009_a7,p22009_a8,p22009_a9,p22009_a10 \
--catCovarList p31 \
--qt  \
--pred $(basename $predFile) \
--bsize 400 \
--minMAC 1 \
--apply-rint \
--out $outFile \
--vc-tests skato,acato-full \
--anno-file $(basename $annoFile).annotations.gz \
--set-list $(basename $annoFile).setlist.gz \
--mask-def $(basename $maskFile) \
--condition-list $(basename $condListFile) \
--condition-file pgen,$(basename $condFile) \
"

echo $cmdText

dx run swiss-army-knife $inputText -icmd="$cmdText" -iimage_file=Docker/regenie.v3.2.8.tar.gz -imount_inputs=false --tag "regenie.step2" --tag $outFile --destination $outFolder --instance-type mem2_ssd1_v2_x8 --brief -y --priority low

