
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png/maf0.1_conditional.png",900,900)

wantedTraits = fread("~/git/noncoding_rvat/supporting_data/42_traits.txt", header = F)

normalData = fread("~/git/noncoding_rvat/revision/maf0.1/NC.AF0.001_SKATO_process_all.tsv.gz")
condData = fread("~/git/noncoding_rvat/revision/maf0.1/NC.conditional_AF0.001_SKATO_process_all.tsv.gz")

mergedData = merge(normalData[,.(CHROM,GENE,LOG10P,PVAL,FDR,TRAIT)],condData[,.(CHROM,GENE,LOG10P,PVAL,FDR,TRAIT)], by = c("CHROM","GENE","TRAIT"))
mergedData = mergedData[TRAIT %in% wantedTraits$V1]


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
test = cor.test(mergedData[FDR.x < 0.05]$LOG10P.x,mergedData[FDR.x < 0.05]$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
ggplot(data = mergedData[FDR.x < 0.05], aes(LOG10P.x, LOG10P.y)) + 
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 2, vjust = 1.5, size = 8, fontface = "bold"  ) +
  ylab("Conditional (-log10)") +
  xlab("No conditional (-log10)") +
  theme_minimal() + 
  theme(text = element_text(size=36), legend.position = c(0.36,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  


dev.off()
