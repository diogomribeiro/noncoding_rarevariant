
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png/intersect_union.png",1600,1600)

intersectDataBefore = fread("~/git/noncoding_rvat/revision/intersection/RESULTS/INTERSECT_SKATO_process_all.tsv")
atleast2DataBefore = fread("~/git/noncoding_rvat/revision/intersection/RESULTS/ATLEAST2_SKATO_process_all.tsv")

intersectDataAfter = fread("~/git/noncoding_rvat/revision/intersection/RESULTS/INTERSECT_SKATO_conditional_all.tsv")
atleast2DataAfter = fread("~/git/noncoding_rvat/revision/intersection/RESULTS/ATLEAST2_SKATO_conditional_all.tsv")

checkDataBefore = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/NC_SKATO_process_all.tsv.gz")
checkDataAfter = fread("~/git/noncoding_rvat/revision/cond_rerun/NC_vcmax_SKATO_conditional_all.tsv.gz")

intersectDataBefore$tag = paste(intersectDataBefore$GENE,intersectDataBefore$TRAIT,sep="|")
atleast2DataBefore$tag = paste(atleast2DataBefore$GENE,atleast2DataBefore$TRAIT,sep="|")
intersectDataAfter$tag = paste(intersectDataAfter$GENE,intersectDataAfter$TRAIT,sep="|")
atleast2DataAfter$tag = paste(atleast2DataAfter$GENE,atleast2DataAfter$TRAIT,sep="|")
checkDataBefore$tag = paste(checkDataBefore$GENE,checkDataBefore$TRAIT,sep="|")
checkDataAfter$tag = paste(checkDataAfter$GENE,checkDataAfter$TRAIT,sep="|")

checkDataBefore[FDR < 0.05]
intersectDataBefore[FDR < 0.05]
atleast2DataBefore[FDR < 0.05]
checkDataBefore[LOG10P > 9]
intersectDataBefore[LOG10P > 9]
atleast2DataBefore[LOG10P > 9]

checkDataAfter$keepFDR = 0
checkDataAfter[FDR < 0.05][tag %in% checkDataBefore[FDR < 0.05]$tag]$keepFDR = 1

intersectDataAfter$keepFDR = 0 
intersectDataAfter[FDR < 0.05][tag %in% intersectDataBefore[FDR < 0.05]$tag]$keepFDR = 1

atleast2DataAfter$keepFDR = 0
atleast2DataAfter[FDR < 0.05][tag %in% atleast2DataBefore[FDR < 0.05]$tag]$keepFDR = 1

checkDataAfter$keepP1e9 = 0
checkDataAfter[LOG10P > 9][tag %in% checkDataBefore[LOG10P > 9]$tag]$keepP1e9 = 1

intersectDataAfter$keepP1e9 = 0
intersectDataAfter[LOG10P > 9][tag %in% intersectDataBefore[LOG10P > 9]$tag]$keepP1e9 = 1

atleast2DataAfter$keepP1e9 = 0
atleast2DataAfter[LOG10P > 9][tag %in% atleast2DataBefore[LOG10P > 9]$tag]$keepP1e9 = 1

mergeData1 = merge(checkDataBefore,intersectDataBefore, by = c("GENE","TRAIT"))
mergeData2 = merge(checkDataBefore,atleast2DataBefore, by = c("GENE","TRAIT"))

logMax1 = max(mergeData1$LOG10P.x, mergeData1$LOG10P.y)
logMax2 = max(mergeData2$LOG10P.x, mergeData2$LOG10P.y)


test = cor.test(mergeData2$LOG10P.x,mergeData2$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,3), "P-value",format.pval(test$p.value,2))
g1 = ggplot(data = mergeData2, aes(LOG10P.x, LOG10P.y)) + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 1) +
  geom_hline(yintercept = -log10(1e-9), color = "#de2d26", linetype = "dashed", size = 1) +
  geom_vline(xintercept = -log10(1e-9), color = "#de2d26", linetype = "dashed", size = 1) +
  geom_point(size = 3) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 9, fontface = "bold"  ) +
  ggtitle("Partial intersect | Before conditional") + 
  ylab("Partial intersect (-log10)") +
  xlab("Union (-log10)") +
  xlim(0,logMax2) +
  ylim(0,logMax2) +
  labs(tag = "a") +
  theme_minimal() + 
  theme(text = element_text(size=28), legend.position = c(0.86,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  

test = cor.test(mergeData1$LOG10P.x,mergeData1$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,3), "P-value",format.pval(test$p.value,2))
g2 = ggplot(data = mergeData1, aes(LOG10P.x, LOG10P.y)) + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 1) +
  geom_hline(yintercept = -log10(1e-9), color = "#de2d26", linetype = "dashed", size = 1) +
  geom_vline(xintercept = -log10(1e-9), color = "#de2d26", linetype = "dashed", size = 1) +
  geom_point(size = 3) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 9, fontface = "bold"  ) +
  ylab("Full intersect (-log10)") +
  xlab("Union (-log10)") +
  ggtitle("Full intersect | Before conditional") + 
  xlim(0,logMax1) +
  ylim(0,logMax1) +
  labs(tag = "b") +
  theme_minimal() + 
  theme(text = element_text(size=28), legend.position = c(0.86,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  

mergeData1 = merge(checkDataAfter,intersectDataAfter, by = c("GENE","TRAIT"))
mergeData2 = merge(checkDataAfter,atleast2DataAfter, by = c("GENE","TRAIT"))

logMax1 = max(mergeData1$LOG10P.x, mergeData1$LOG10P.y)
logMax2 = max(mergeData2$LOG10P.x, mergeData2$LOG10P.y)

test = cor.test(mergeData2$LOG10P.x,mergeData2$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,3), "P-value",format.pval(test$p.value,2))
g3 = ggplot(data = mergeData2, aes(LOG10P.x, LOG10P.y)) + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 1) +
  geom_hline(yintercept = -log10(1e-9), color = "#de2d26", linetype = "dashed", size = 1) +
  geom_vline(xintercept = -log10(1e-9), color = "#de2d26", linetype = "dashed", size = 1) +
  geom_point(size = 3) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 9, fontface = "bold"  ) +
  ggtitle("Partial intersect | After conditional") + 
  ylab("Partial intersect (-log10)") +
  xlab("Union (-log10)") +
  xlim(0,logMax2) +
  ylim(0,logMax2) +
  labs(tag = "c") +
  theme_minimal() + 
  theme(text = element_text(size=28), legend.position = c(0.86,0.85), plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  

test = cor.test(mergeData1$LOG10P.x,mergeData1$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,3), "P-value",format.pval(test$p.value,2))
g4 = ggplot(data = mergeData1, aes(LOG10P.x, LOG10P.y)) + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 1) +
  geom_hline(yintercept = -log10(1e-9), color = "#de2d26", linetype = "dashed", size = 1) +
  geom_vline(xintercept = -log10(1e-9), color = "#de2d26", linetype = "dashed", size = 1) +
  geom_point(size = 3) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 9, fontface = "bold"  ) +
  ylab("Full intersect (-log10)") +
  xlab("Union (-log10)") +
  ggtitle("Full intersect | After conditional") + 
  xlim(0,logMax1) +
  ylim(0,logMax1) +
  labs(tag = "d") +
  theme_minimal() + 
  theme(text = element_text(size=28), legend.position = c(0.86,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  

g1 + g2 + g3 + g4

dev.off()


# # ADD TRAIT NAME AND GENE NAME # THIS CAN APPLY TRAIT FILTER
# traitMapping = fread("~/git/noncoding_rvat/supporting_data/42_traits_mapping.txt", header = F)
# traitMapping$TRAIT = data.table(unlist(lapply(traitMapping$V1, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
# geneMapping = fread("~/git/noncoding_rvat/supporting_data/gencode_id_mapping.tsv.gz")
# intersectData = merge(intersectData, traitMapping[,.(TRAIT,V2)], by = "TRAIT")
# intersectData = merge(intersectData, geneMapping[,.(V4,V6)], by.x = "GENE", by.y = "V4")
