
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png/conditional_results.png",2000,1000)

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

###

normalData = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/CDS_SKATO_process_all.tsv.gz")
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

g1 + g2

dev.off()

