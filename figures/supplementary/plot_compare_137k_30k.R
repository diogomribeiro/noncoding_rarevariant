
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png//replication_pval.png",1600,900)

cdsData137k = fread("zcat ~/burdenNC/data/results/137k_30k/CDS_137k_SKATO_process_all.tsv.gz")
cdsData30k = fread("zcat  ~/burdenNC/data/results/137k_30k/CDS_30k_SKATO_process_all.tsv.gz")
ncData137k = fread("zcat ~/burdenNC/data/results/137k_30k/NC_137k_SKATO_process_all.tsv.gz")
ncData30k = fread("zcat  ~/burdenNC/data/results/137k_30k/NC_30k_SKATO_process_all.tsv.gz")

wantedTraits = fread("~/git/noncoding_rvat/supporting_data/42_traits.txt", header = F)

cdsData137k = cdsData137k[TRAIT %in% wantedTraits$V1]
cdsData30k = cdsData30k[TRAIT %in% wantedTraits$V1]
ncData137k = ncData137k[TRAIT %in% wantedTraits$V1]
ncData30k = ncData30k[TRAIT %in% wantedTraits$V1]

###########
# NC
###########
mergedNC = merge(ncData30k,ncData137k,by=c("GENE","TRAIT"), all.x = T, all.y = T)
paste(nrow(unique(mergedNC[is.na(PVAL.x)])), "NC gene-traits not evaluated in 30k")
paste(nrow(unique(mergedNC[is.na(PVAL.y)])), "NC gene-traits not evaluated in 137k")
paste(nrow(unique(mergedNC[!is.na(PVAL.x)][!is.na(PVAL.y)])),"NC gene-traits compared")

mergedNC = mergedNC[!is.na(PVAL.x)][!is.na(PVAL.y)]
test = cor.test(-log10(mergedNC[-log10(PVAL.y) > 3]$PVAL.x),-log10(mergedNC[-log10(PVAL.y) > 3]$PVAL.y), method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
g1 = ggplot(data = mergedNC[-log10(PVAL.y) > 3], aes(-log10(PVAL.x), -log10(PVAL.y))) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 8, fontface = "bold"  ) +
  ylab("137k samples (-log10)") +
  xlab("30k samples (-log10)") +
  labs(tag = "a", title = "NC") + 
  theme_minimal() + 
  theme(text = element_text(size=36), legend.position = c(0.86,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

###########
## CDS
###########
mergedCDS = merge(cdsData30k,cdsData137k,by=c("GENE","TRAIT"), all.x = T, all.y = T)
paste(nrow(unique(mergedCDS[is.na(PVAL.x)])), "CDS gene-traits not evaluated in 30k")
paste(nrow(unique(mergedCDS[is.na(PVAL.y)])), "CDS gene-traits not evaluated in 137k")
paste(nrow(unique(mergedCDS[!is.na(PVAL.x)][!is.na(PVAL.y)])),"NC gene-traits compared")

mergedCDS = mergedCDS[!is.na(PVAL.x)][!is.na(PVAL.y)]
test = cor.test(-log10(mergedCDS[-log10(PVAL.y) > 3]$PVAL.x),-log10(mergedCDS[-log10(PVAL.y) > 3]$PVAL.y), method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
g2 = ggplot(data = mergedCDS[-log10(PVAL.y) > 3], aes(-log10(PVAL.x), -log10(PVAL.y))) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 8, fontface = "bold"  ) +
  ylab("137k samples (-log10)") +
  xlab("30k samples (-log10)") +
  labs(tag = "b", title = "CDS") + 
  theme_minimal() + 
  theme(text = element_text(size=36), legend.position = c(0.86,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  


g1 + g2

dev.off()
