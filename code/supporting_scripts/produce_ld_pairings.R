
library(data.table)

topVars = fread("dnanexus/NC.single.cases.FDR0.05.tsv")
topVars$TRAIT = data.table(unlist(lapply(topVars$TRAIT, function(x) unlist(strsplit(x,"p"))[2])))$V1

cojoData = fread("dnanexus/all_trait.conditional_hg38.cojo", header = F)
cojoData$TRAIT = data.table(unlist(lapply(cojoData$V5, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
cojoData$CHROM = data.table(unlist(lapply(cojoData$V4, function(x) unlist(strsplit(x,"[_]"))[1])))$V1

# Note: WGS data (top single variants) is in hg38, but SNP array data (conditional) are in hg19.
# Here we use a map between hg19 and hg38 for calculating distance between variants (1Mb), but still print the "mismatched" coordinates in the output

pairing = data.table()
for (trait in unique(topVars$TRAIT)){
  print(trait)
  tops = topVars[TRAIT == trait]
  cojo = cojoData[TRAIT == trait]
  for (var1 in tops$VAR){
    for (j in seq(nrow(cojo))){
      var1Chrom = strsplit(var1,"[:]")[[1]][1]
      var1Pos = as.numeric(strsplit(var1,"[:]")[[1]][2])
      cojoVar = cojo[j,]
      var2 = cojoVar$V4
      var2Chrom = cojoVar$V1
      var2Pos = cojoVar$V2 # hg38 coord to match var1
      diff = abs(var1Pos - var2Pos)

      if (var1Chrom == var2Chrom & diff < 100000){
        pairing = rbind(pairing, data.table(topVar = var1, cojoVar = var2))
      }
    }
  }
}

pairing

write.table(unique(pairing), "dnanexus/ld_pairings_NC_fdr0.05.tsv", quote=F , sep= "\t", row.names=F, col.names=F)

