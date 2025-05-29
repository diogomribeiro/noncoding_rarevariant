inputText="
-iin=/Output/SNP_array/200K/ukb_200K_sampleprune_autosome_ldprune.bed
-iin=/Output/SNP_array/200K/ukb_200K_sampleprune_autosome_ldprune.bim
-iin=/Output/SNP_array/200K/ukb_200K_sampleprune_autosome_ldprune.fam
-iin=regenie/200k/regenie_phenotype_21July2023.txt
-iin=regenie/200k/regenie_covariate_21July2023.txt
-iin=/Output/SNP_array/200K/200k_qc_pass.id
-iin=/Output/SNP_array/200K/200k_qc_pass.snplist
"

cmdText="
regenie \
  --step 1 \
  --bed ukb_200K_sampleprune_autosome_ldprune \
  --extract 200k_qc_pass.snplist \
  --keep 200k_qc_pass.id \
  --phenoFile regenie_phenotype_21July2023.txt \
  --phenoColList p30160_i0,p30220_i0,p30150_i0,p30210_i0,p30030_i0,p30020_i0,p30300_i0,p30290_i0,p30280_i0,p30120_i0,p30180_i0,p30050_i0,p30060_i0,p30040_i0,p30100_i0,p30260_i0,p30270_i0,p30130_i0,p30190_i0,p30140_i0,p30200_i0,p30170_i0,p30230_i0,p30080_i0,p30090_i0,p30110_i0,p30010_i0,p30070_i0,p30250_i0,p30240_i0,p30000_i0,p30620_i0,p30600_i0,p30610_i0,p30630_i0,p30640_i0,p30650_i0,p30710_i0,p30680_i0,p30690_i0,p30700_i0,p30720_i0,p30660_i0,p30730_i0,p30740_i0,p30750_i0,p30760_i0,p30770_i0,p30780_i0,p30790_i0,p30800_i0,p30810_i0,p30820_i0,p30830_i0,p30850_i0,p30840_i0,p30860_i0,p30870_i0,p30880_i0,p30670_i0,p30890_i0,p50_i0,p20015_i0,p21001_i0,p23099_i0,p23105_i0,p48_i0,p49_i0,p21002_i0,p20022_i0,p22191_i0,p4080_i0_a0,p4079_i0_a0,p102_i0_a0 \
  --covarFile regenie_covariate_21July2023.txt \
  --covarColList p31,p21022,p22009_a1,p22009_a2,p22009_a3,p22009_a4,p22009_a5,p22009_a6,p22009_a7,p22009_a8,p22009_a9,p22009_a10 \
  --catCovarList p31 \
  --qt \
  --apply-rint \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix regenie_tmp_preds \
  --out ukb_200k_step1
"

dx run swiss-army-knife $inputText -icmd="$cmdText" -iimage_file=:Docker/regenie.v3.2.8.tar.gz -imount_inputs=true --tag "regenie.step1" --destination regenie/200k/step1/ --instance-type mem2_ssd1_v2_x8 --priority low --brief -y

