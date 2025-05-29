
library(data.table)
data = fread("~/git/noncoding_rvat/revision/maf0.1/NC.CADD15.annotations.gz",header = F)
pvar = fread("~/git/noncoding_rvat/revision/maf0.1/ukb200k_all_chr.pvar")
samples = 166740*2

mergedData = merge(data, pvar, by.x = "V1", by.y = "ID")

#test = mergedData[`#CHROM` == 22]

mergedData$AC = unlist(lapply(mergedData$INFO, function(x) strsplit(x,"=")[[1]][2]))

mergedData$MAF = as.numeric(mergedData$AC) / samples * 100

mergedDataFilt = mergedData[MAF < 0.1]

freq = data.table(table(mergedDataFilt$V2))

summary(freq$N)

# MAF 0.1% 554 median

# freq2 = data.table(table(mergedData$V2))
# summary(freq2$N)
# # MAF 1% 563 median 