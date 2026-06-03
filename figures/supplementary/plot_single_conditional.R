# Script to classify gene-trait associations as single or group
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png/single_conditional_results.png",1800,900)

# NC
dataCond = fread("~/git/noncoding_rvat/revision/LD/before_conditional/single_conditional/nc_single_cond_results.regenie")
dataCheck = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/NC_SKATO_process_FDR0.05.tsv")

dataCond$TRAIT = data.table(unlist(lapply(dataCond$V1, function(x) unlist(strsplit(x,"[_]"))[3])))$V1
dataCond = dataCond[,.(TRAIT,V3,V12)]
dataCond$GENE = data.table(unlist(lapply(dataCond$V3, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
colnames(dataCond) = c("TRAIT","V3","LOG10P","GENE")
dataCheck = dataCheck[,.(TRAIT,GENE,LOG10P)]

mergedData = merge(dataCheck, dataCond, by = c("GENE","TRAIT"))

g1 = ggplot(data = mergedData, aes(x = LOG10P.x, y = LOG10P.y)) +
  geom_hline(yintercept = -log10(6.50E-05), color = "red", linetype = "dashed", size = 2) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) + 
  ylab("Single conditional (-log10)") +
  xlab("No conditional (-log10)") +
  labs(tag = "a") + 
  theme_minimal() +
  theme(text = element_text(size=36), legend.position = c(0.3,0.9),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)

## CDS
dataCond = fread("~/git/noncoding_rvat/revision/LD/before_conditional/single_conditional/cds_single_cond_results.regenie")
dataCheck = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/CDS_SKATO_process_FDR0.05.tsv")

dataCond$TRAIT = data.table(unlist(lapply(dataCond$V1, function(x) unlist(strsplit(x,"[_]"))[3])))$V1
dataCond = dataCond[,.(TRAIT,V3,V12)]
dataCond$GENE = data.table(unlist(lapply(dataCond$V3, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
colnames(dataCond) = c("TRAIT","V3","LOG10P","GENE")
dataCheck = dataCheck[,.(TRAIT,GENE,LOG10P)]

mergedData = merge(dataCheck, dataCond, by = c("GENE","TRAIT"))

g2 = ggplot(data = mergedData, aes(x = LOG10P.x, y = LOG10P.y)) +
  geom_hline(yintercept = -log10(4.60E-05), color = "red", linetype = "dashed", size = 2) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) + 
  ylab("Single conditional (-log10)") +
  xlab("No conditional (-log10)") +
  labs(tag = "b") + 
  theme_minimal() +
  theme(text = element_text(size=36), legend.position = c(0.3,0.9),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)


g1 + g2

dev.off()