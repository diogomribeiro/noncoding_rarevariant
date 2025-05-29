
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png/ld_conditional_results_v2.png",1500,2250)

nocond = fread("~/git/noncoding_rvat/revision/LD/before_conditional/NC.FDR0.05.single.txt")

cond = fread("~/git/noncoding_rvat/revision/LD/before_conditional/single_conditional/LD/NC.FDR5.LD.cond.regenie", sep = " ")
cond$VAR = lapply(cond$`EXTRA	VAR`, function(x) strsplit(x,"[.]")[[1]][3])
cond$VAR1 = lapply(cond$VAR, function(x) strsplit(x,"[_]")[[1]][1])
cond$VAR2 = lapply(cond$VAR, function(x) strsplit(x,"[_]")[[1]][2])
cond$VAR3 = lapply(cond$VAR, function(x) strsplit(x,"[_]")[[1]][3])
cond$VAR4 = lapply(cond$VAR, function(x) strsplit(x,"[_]")[[1]][4])
cond$VARIANT = paste(cond$VAR1,cond$VAR2,cond$VAR3,cond$VAR4,sep=":")
cond$TRAIT = unlist(lapply(cond$VAR, function(x) strsplit(x,"[_]")[[1]][5]))
cond = cond[ID == VARIANT]
cond = cond[,.(ID,A1FREQ,LOG10P,TRAIT)]

fullCond = fread("~/git/noncoding_rvat/revision/LD/before_conditional/single_conditional/LD/NC.single.cond.wanted.regenie")
fullCond$TRAIT = unlist(lapply(fullCond$V5, function(x) unlist(strsplit(x,"[_]"))[2]))
fullCond = fullCond[,.(V1,V3,TRAIT)]
colnames(fullCond) = c("VAR","FULL_LOG10P","TRAIT")

mergedData = merge(nocond, cond, by.x = c("VAR","TRAIT"), by.y = c("ID","TRAIT"))
mergedData = merge(mergedData, fullCond, by = c("VAR","TRAIT"))
mergedData = unique(mergedData[,.(VAR,TRAIT,LOG10P.x,LOG10P.y,FULL_LOG10P)])
length(unique(mergedData$VAR))
length(unique(mergedData$TRAIT))

## LD conditional
# When doing conditional, most p-values get weaker, as seen previously
100*nrow(mergedData[LOG10P.x > LOG10P.y])/nrow(mergedData)
# Perc of FDR 5% associations that remain after conditioning for variants
minLOG10 = -log10(6.50E-05) #min(mergedData$LOG10P.x)
100 - 100 * nrow(mergedData[LOG10P.x > minLOG10 & LOG10P.y > minLOG10]) / nrow(mergedData[LOG10P.x > minLOG10])
# Perc of 1e-9 associations that remain after conditioning for variants in LD
100 - 100 * nrow(mergedData[LOG10P.x > 9 & LOG10P.y > 9]) / nrow(mergedData[LOG10P.x > 9])

## FULL conditional
# When doing conditional, most p-values get weaker, as seen previously
100*nrow(mergedData[LOG10P.x > FULL_LOG10P])/nrow(mergedData)
# Perc of FDR 5% associations that remain after conditioning for variants
100 - 100 * nrow(mergedData[LOG10P.x > minLOG10 & FULL_LOG10P > minLOG10]) / nrow(mergedData[LOG10P.x > minLOG10])
# Perc of 1e-9 associations that remain after conditioning for variants
100 - 100 * nrow(mergedData[LOG10P.x > 9 & FULL_LOG10P > 9]) / nrow(mergedData[LOG10P.x > 9])

# No cond vs LD cond
test = cor.test(mergedData$LOG10P.x,mergedData$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
g1 = ggplot(data = mergedData, aes(LOG10P.x, LOG10P.y)) +
  geom_hline(yintercept = -log10(6.50E-05), color = "red", linetype = "dashed", size = 2) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 8, fontface = "bold"  ) +
  ylab("LD conditional (-log10)") +
  xlab("No conditional (-log10)") +
  labs(tag = "a") +
  ggtitle("LD conditional | P-values") +
  theme_minimal() +
  theme(text = element_text(size=36), legend.position = c(0.86,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)


# NO cond vs Full cond
test = cor.test(mergedData$LOG10P.x,mergedData$FULL_LOG10P, method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
g2 = ggplot(data = mergedData, aes(LOG10P.x, FULL_LOG10P)) +
  geom_hline(yintercept = -log10(6.50E-05), color = "red", linetype = "dashed", size = 2) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 8, fontface = "bold"  ) +
  ylab("Full conditional (-log10)") +
  xlab("No conditional (-log10)") +
  labs(tag = "b") +
  ggtitle("Full conditional | P-values") +
  theme_minimal() +
  theme(text = element_text(size=36), legend.position = c(0.86,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)


nocond$AF = as.numeric(formatC(nocond$AF, format = "e", digits = 0))
mergedData = merge(mergedData,unique(nocond[,.(VAR,AF)]), by = "VAR")
remain = mergedData[LOG10P.x > minLOG10 & LOG10P.y > minLOG10]
nolonger = mergedData[LOG10P.x > minLOG10 & LOG10P.y < minLOG10]
summary(remain$AF)
summary(nolonger$AF)
wilcox.test(remain$AF, nolonger$AF)
remain$conditional = "Significant"
nolonger$conditional = "Not significant"
afData = rbind(remain,nolonger)
g3 = ggplot(data = afData, aes(x = AF, fill = conditional)) +
  geom_density(color = "black", alpha = 0.8, adjust = 1.5) +
  xlab("Allele frequency") +
  scale_fill_brewer(palette = "Set2") + 
  ggtitle("LD conditional | Allele freq.") +
  labs(tag = "c") +
  theme_minimal() +
  theme(text = element_text(size=36), legend.position = c(0.4,0.9), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)

remain = mergedData[LOG10P.x > minLOG10 & FULL_LOG10P > minLOG10]
nolonger = mergedData[LOG10P.x > minLOG10 & FULL_LOG10P < minLOG10]
summary(remain$AF)
summary(nolonger$AF)
wilcox.test(remain$AF, nolonger$AF)
remain$conditional = "Significant"
nolonger$conditional = "Not significant"

afData = rbind(remain,nolonger)
g4 = ggplot(data = afData, aes(x = AF, fill = conditional)) +
  geom_density(color = "black", alpha = 0.8, adjust = 1.5) +
  xlab("Allele frequency") +
  scale_fill_brewer(palette = "Set2") + 
  ggtitle("Full conditional | Allele freq.") +
  labs(tag = "d") +
  theme_minimal() +
  theme(text = element_text(size=36), legend.position = c(0.5,0.8), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)


#########################
# NEW PANELS
#########################

minLOG10 = -log10(6.50E-05) #min(mergedData$LOG10P.x)

nocond = fread("~/git/noncoding_rvat/revision/LD/before_conditional/NC.FDR0.05.single.txt")

cond = fread("~/git/noncoding_rvat/revision/LD/before_conditional/single_conditional/LD/NC.FDR5.LD.cond.regenie", sep = " ")
cond$VAR = lapply(cond$`EXTRA	VAR`, function(x) strsplit(x,"[.]")[[1]][3])
cond$VAR1 = lapply(cond$VAR, function(x) strsplit(x,"[_]")[[1]][1])
cond$VAR2 = lapply(cond$VAR, function(x) strsplit(x,"[_]")[[1]][2])
cond$VAR3 = lapply(cond$VAR, function(x) strsplit(x,"[_]")[[1]][3])
cond$VAR4 = lapply(cond$VAR, function(x) strsplit(x,"[_]")[[1]][4])
cond$VARIANT = paste(cond$VAR1,cond$VAR2,cond$VAR3,cond$VAR4,sep=":")
cond$TRAIT = unlist(lapply(cond$VAR, function(x) strsplit(x,"[_]")[[1]][5]))
cond = cond[ID == VARIANT]
cond = cond[,.(ID,A1FREQ,LOG10P,TRAIT)]

fullCond = fread("~/git/noncoding_rvat/revision/LD/before_conditional/single_conditional/LD/NC.single.cond.wanted.regenie")
fullCond$TRAIT = unlist(lapply(fullCond$V5, function(x) unlist(strsplit(x,"[_]"))[2]))
fullCond = fullCond[,.(V1,V3,TRAIT)]
colnames(fullCond) = c("VAR","FULL_LOG10P","TRAIT")

mergedData = merge(nocond, cond, by.x = c("VAR","TRAIT"), by.y = c("ID","TRAIT"))
mergedData = merge(mergedData, fullCond, by = c("VAR","TRAIT"))
mergedData = unique(mergedData[,.(VAR,TRAIT,LOG10P.x,LOG10P.y,FULL_LOG10P)])
length(unique(mergedData$VAR))
length(unique(mergedData$TRAIT))

nocond$AF = as.numeric(formatC(nocond$AF, format = "e", digits = 0))
mergedData = merge(mergedData,unique(nocond[,.(VAR,AF)]), by = "VAR")

### Remain / no longer based on LD variant conditioning
# 127 variant-trait associations remaining significant at FDR < 5% and the 79 associations losing significance with LD variant conditioning
remain = mergedData[LOG10P.x > minLOG10 & LOG10P.y > minLOG10]
nolonger = mergedData[LOG10P.x > minLOG10 & LOG10P.y < minLOG10]

### Remain / no longer based on Full conditional
# 63 variant-trait associations remaining significance at FDR < 5% and the 143 associations losing significance when conditioning for all GWAS variants
remain2 = mergedData[LOG10P.x > minLOG10 & FULL_LOG10P > minLOG10]
nolonger2 = mergedData[LOG10P.x > minLOG10 & FULL_LOG10P < minLOG10]

remain$conditional = "Significant"
nolonger$conditional = "Not significant"
remain2$conditional = "Significant"
nolonger2$conditional = "Not significant"

##############################
# read GWAS variants per trait. Perform analysis per topVar and trait (only using LD for GWAS variants used for trait)
##############################

ldData = fread("~/git/noncoding_rvat/revision/LD/before_conditional/single_conditional/LD/supplementary_table_ld.tsv")
# This LD data includes LD for all combinations between top Vars and conditional vars (nearby), regardless of trait

## Read Conditional variants per trait
conditionalData = fread("~/git/noncoding_rvat/revision/LD/conditional_GWAS_variants/all_trait.conditional.cojo", header = F)
conditionalData[, V2 := gsub("(\\d+)\\.conditional\\.cojo", "p\\1", V2)]

ldData = merge(ldData,conditionalData, by.x = "condVar", by.y = "V1", all.x = T)
colnames(ldData) = c("condVar","topVar","R2","D_prime","???","trait")

remainLD = merge(ldData, remain, by.x = c("topVar","trait"), by.y = c("VAR","TRAIT"), all.y = T)
nolongerLD = merge(ldData, nolonger, by.x = c("topVar","trait"), by.y = c("VAR","TRAIT"), all.y = T)
remainLD2 = merge(ldData, remain2, by.x = c("topVar","trait"), by.y = c("VAR","TRAIT"), all.y = T)
nolongerLD2 = merge(ldData, nolonger2, by.x = c("topVar","trait"), by.y = c("VAR","TRAIT"), all.y = T)

### Calculating mean per variant-trait
result1 <- remainLD[, .(max_R2 = max(R2)), by = .(topVar,trait)]
result2 <- nolongerLD[, .(max_R2 = max(R2)), by = .(topVar,trait)]
summary(result1$max_R2)
summary(result2$max_R2)
w = wilcox.test(result1$max_R2,result2$max_R2)
w$p.value
ks.test(result1$max_R2,result2$max_R2)
result1$conditional = "Significant"
result2$conditional = "Not significant"
result = rbind(result1,result2)

g5 = ggplot(data = result, aes(x = conditional, y = max_R2, fill = conditional)) +
  geom_boxplot() +
  # geom_density(color = "black", alpha = 0.8, adjust = 2) +
  ylab("LD") +
  xlab("Conditional") +
  scale_fill_brewer(palette = "Set2") +
  ggtitle("LD conditional | LD") +
  labs(tag = "e") +
  theme_minimal() +
  theme(text = element_text(size=36), legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)

### Calculating mean per variant-trait
result1 <- remainLD2[, .(max_R2 = max(R2)), by = .(topVar,trait)]
result2 <- nolongerLD2[, .(max_R2 = max(R2)), by = .(topVar,trait)]
summary(result1$max_R2)
summary(result2$max_R2)
w = wilcox.test(result1$max_R2,result2$max_R2)
w$p.value
ks.test(result1$max_R2,result2$max_R2)
result1$conditional = "Significant"
result2$conditional = "Not significant"
result2 = rbind(result1,result2)

g6 = ggplot(data = result2, aes(x = conditional, y = max_R2, fill = conditional)) +
  geom_boxplot() +
  # geom_density(color = "black", alpha = 0.8, adjust = 2) +
  ylab("LD") +
  xlab("Conditional") +
  scale_fill_brewer(palette = "Set2") +
  ggtitle("Full conditional | LD") +
  labs(tag = "f") +
  theme_minimal() +
  theme(text = element_text(size=36), legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)


(g1 | g2) / (g3 | g4) / ( g5 | g6 )
dev.off()
