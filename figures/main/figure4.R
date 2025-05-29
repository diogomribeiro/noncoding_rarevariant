
library(data.table)
library(ggplot2)
library(patchwork)

# png("~/git/noncoding_rvat/figures/final/figure4.png",2000,2000)
pdf("~/git/noncoding_rvat/figures/final/figure4.pdf",28,28)

##############################
# Panel A
##############################

wantedTraits = fread("~/git/noncoding_rvat/supporting_data/42_traits.txt", header = F)

normalData = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/NC_SKATO_process_all.tsv.gz")
condData = fread("~/git/noncoding_rvat/revision/cond_rerun/NC_vcmax_SKATO_conditional_all.tsv.gz")
mergedData = merge(normalData[,.(CHROM,GENE,LOG10P,PVAL,FDR,TRAIT)],condData[,.(CHROM,GENE,LOG10P,PVAL,FDR,TRAIT)], by = c("CHROM","GENE","TRAIT"))
mergedData = mergedData[TRAIT %in% wantedTraits$V1]

mergedData = mergedData[FDR.x < 0.05]
# Gene-trait associations
# Before conditional
nrow(mergedData[FDR.x < 0.05])
nrow(mergedData[PVAL.x < 1e-9])
# After conditional
nrow(mergedData[FDR.y < 0.05])
nrow(mergedData[PVAL.y < 1e-9])

# % remaining
remainFDR5 = nrow(mergedData[FDR.x < 0.05][FDR.y < 0.05])
remainP1e9 = nrow(mergedData[PVAL.x < 1e-9][PVAL.y < 1e-9])
round(remainFDR5 * 100 / nrow(mergedData[FDR.x < 0.05]),2)
round(remainP1e9) * 100 / nrow(mergedData[PVAL.x < 1e-9])

# Plot
test = cor.test(mergedData$LOG10P.x,mergedData$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
g1 = ggplot(data = mergedData, aes(LOG10P.x, LOG10P.y)) + 
  geom_point() +
  geom_point(x = mergedData[GENE == "ENSG00000101439"][TRAIT == "p30720"]$LOG10P.x, y = mergedData[GENE == "ENSG00000101439"][TRAIT == "p30720"]$LOG10P.y, size=14, shape=1, color="#006d2c") + 
  geom_text(x = mergedData[GENE == "ENSG00000101439"][TRAIT == "p30720"]$LOG10P.x + 12, y = mergedData[GENE == "ENSG00000101439"][TRAIT == "p30720"]$LOG10P.y - 3, size=11, label = "CST3 / Cystatin C", color="#006d2c") + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  # annotate("text", x = Inf, y = Inf, label = text, hjust = 1.2, vjust = 1.5, size = 8, fontface = "bold"  ) +
  ylab("Conditional p-value (-log10)") +
  xlab("No conditional p-value (-log10)") +
  labs(tag = "a") + 
  theme_minimal() + 
  theme(text = element_text(size=36), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

# Source data export
colnames(mergedData) = c("CHROM","GENE","TRAIT","LOG10P_nocond","PVAL_nocond","FDR_nocond","LOG10P_cond","PVAL_cond","FDR_cond")
d1 = mergedData

##############################
# Panel B
##############################

normalData = fread("/home/dribeiro/git/noncoding_rvat/revision/no_cond_42_traits/CDS_SKATO_process_all.tsv.gz")
condData = fread("~/git/noncoding_rvat/revision/cond_rerun/CDS_vcmax_SKATO_conditional_all.tsv.gz")
mergedData = merge(normalData[,.(CHROM,GENE,LOG10P,PVAL,FDR,TRAIT)],condData[,.(CHROM,GENE,LOG10P,PVAL,FDR,TRAIT)], by = c("CHROM","GENE","TRAIT"))
mergedData = mergedData[TRAIT %in% wantedTraits$V1]

mergedData = mergedData[FDR.x < 0.05]
# Gene-trait associations
# Before conditional
nrow(mergedData[FDR.x < 0.05])
nrow(mergedData[PVAL.x < 1e-9])
# After conditional
nrow(mergedData[FDR.y < 0.05])
nrow(mergedData[PVAL.y < 1e-9])

# % remaining
remainFDR5 = nrow(mergedData[FDR.x < 0.05][FDR.y < 0.05])
remainP1e9 = nrow(mergedData[PVAL.x < 1e-9][PVAL.y < 1e-9])
round(remainFDR5 * 100 / nrow(mergedData[FDR.x < 0.05]),2)
round(remainP1e9) * 100 / nrow(mergedData[PVAL.x < 1e-9])

# Plot
test = cor.test(mergedData$LOG10P.x,mergedData$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
g2 = ggplot(data = mergedData, aes(LOG10P.x, LOG10P.y)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  # annotate("text", x = Inf, y = Inf, label = text, hjust = 1.2, vjust = 1.5, size = 8, fontface = "bold"  ) +
  ylab("Conditional p-value (-log10)") +
  xlab("No conditional p-value (-log10)") +
  labs(tag = "b") + 
  theme_minimal() + 
  theme(text = element_text(size=36), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

# Source data export
colnames(mergedData) = c("CHROM","GENE","TRAIT","LOG10P_nocond","PVAL_nocond","FDR_nocond","LOG10P_cond","PVAL_cond","FDR_cond")
d2 = mergedData

##############################
# Panel C, D
##############################

## Conditional
NCfdr5 = fread("~/git/noncoding_rvat/revision/LD/before_conditional/single_conditional/NC.cond.reclassification_FDR0.05.txt")
CDSfdr5 = fread("~/git/noncoding_rvat/revision/LD/before_conditional/single_conditional/CDS.cond.reclassification_FDR0.05.txt")

filterNCP1e9 = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/NC_SKATO_process_PVAL1e9.tsv", header = F)
filterNCP1e9$tag = paste(filterNCP1e9$V1, filterNCP1e9$V2)
NCp1e9 = NCfdr5[tag %in% filterNCP1e9$tag]

filterCDSP1e9 = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/CDS_SKATO_process_PVAL1e9.tsv", header = F)
filterCDSP1e9$tag = paste(filterCDSP1e9$V1, filterCDSP1e9$V2)
CDSp1e9 = CDSfdr5[tag %in% filterCDSP1e9$tag]


fdr5Data = data.table(values = c(nrow(NCfdr5[CLASS == "Single"]),nrow(NCfdr5[CLASS == "Group"]),nrow(CDSfdr5[CLASS == "Single"]),nrow(CDSfdr5[CLASS == "Group"])), annotation = c("NC","NC","CDS","CDS"), association = c("Single","Group","Single","Group"))
p1e9Data = data.table(values = c(nrow(NCp1e9[CLASS == "Single"]),nrow(NCp1e9[CLASS == "Group"]),nrow(CDSp1e9[CLASS == "Single"]),nrow(CDSp1e9[CLASS == "Group"])), annotation = c("NC","NC","CDS","CDS"), association = c("Single","Group","Single","Group"))

fdr5NCSum = sum(fdr5Data[1:2,]$values)
fdr5CDSSum = sum(fdr5Data[3:4,]$values)
p1e9NCSum = sum(p1e9Data[1:2,]$values)
p1e9CDSSum = sum(p1e9Data[3:4,]$values)

g3 = ggplot(data = fdr5Data, aes(x = annotation, y = values, fill = association)) +
  geom_bar(stat = "identity", color = "black", size = 1.5) +
  annotate(geom = "text", x = "CDS", y = fdr5Data[annotation == "CDS"][association == "Single"]$values, label = paste0(fdr5Data[annotation == "CDS"][association == "Single"]$values," (", round(100 * fdr5Data[annotation == "CDS"][association == "Single"]$values/fdr5CDSSum),"%)"), vjust = 2, size = 10) +
  annotate(geom = "text", x = "CDS", y = fdr5Data[annotation == "CDS"][association == "Single"]$values + 450, label = paste0(fdr5Data[annotation == "CDS"][association == "Group"]$values," (", round(100 * fdr5Data[annotation == "CDS"][association == "Group"]$values/fdr5CDSSum),"%)"), vjust = 2, size = 10) +
  annotate(geom = "text", x = "NC", y = fdr5Data[annotation == "NC"][association == "Single"]$values, label = paste0(fdr5Data[annotation == "NC"][association == "Single"]$values," (", round(100 * fdr5Data[annotation == "NC"][association == "Single"]$values/fdr5NCSum),"%)"), vjust = 2, size = 10) +
  annotate(geom = "text", x = "NC", y = fdr5Data[annotation == "NC"][association == "Single"]$values + 240, label = paste0(fdr5Data[annotation == "NC"][association == "Group"]$values," (", round(100 * fdr5Data[annotation == "NC"][association == "Group"]$values/fdr5NCSum),"%)"), vjust = 2, size = 10) +
  labs(tag = "c") + 
  scale_fill_manual(values = c("#74c476","#fee08b")) +
  ylab("gene-trait associations") +
  xlab("Annotation") +
  theme_minimal() +
  theme(text = element_text(size=36), legend.position = c(0.3,0.9),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)

g4 = ggplot(data = p1e9Data, aes(x = annotation, y = values, fill = association)) + 
  geom_bar(stat = "identity", color = "black", size = 1.5) +
  annotate(geom = "text", x = "CDS", y = p1e9Data[annotation == "CDS"][association == "Single"]$values + 5, label = paste0(p1e9Data[annotation == "CDS"][association == "Single"]$values," (", round(100 * p1e9Data[annotation == "CDS"][association == "Single"]$values/p1e9CDSSum),"%)"), vjust = 3, size = 10) +
  annotate(geom = "text", x = "CDS", y = p1e9Data[annotation == "CDS"][association == "Single"]$values + 155, label = paste0(p1e9Data[annotation == "CDS"][association == "Group"]$values," (", round(100 * p1e9Data[annotation == "CDS"][association == "Group"]$values/p1e9CDSSum),"%)"), vjust = 2, size = 10) +
  annotate(geom = "text", x = "NC", y = p1e9Data[annotation == "NC"][association == "Single"]$values + 5, label = paste0(p1e9Data[annotation == "NC"][association == "Single"]$values," (", round(100 * p1e9Data[annotation == "NC"][association == "Single"]$values/p1e9NCSum),"%)"), vjust = 3, size = 10) +
  annotate(geom = "text", x = "NC", y = p1e9Data[annotation == "NC"][association == "Single"]$values + 80, label = paste0(p1e9Data[annotation == "NC"][association == "Group"]$values," (", round(100 * p1e9Data[annotation == "NC"][association == "Group"]$values/p1e9NCSum),"%)"), vjust = 2, size = 10) +
  labs(tag = "d") + 
  scale_fill_manual(values = c("#74c476","#fee08b")) +
  ylab("gene-trait associations") +
  xlab("Annotation") +
  theme_minimal() + 
  theme(text = element_text(size=36), legend.position = c(0.3,0.9),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  




(g1 + g2) / (g3 + g4)

dev.off()


## Source data
library(openxlsx)
write.xlsx(list(fig4a = d1, fig4b = d2, fig4c = fdr5Data, fig4d = p1e9Data), "~/git/noncoding_rvat/figures/final/source_data/Ribeiro_SourceData_Fig4.xlsx")
