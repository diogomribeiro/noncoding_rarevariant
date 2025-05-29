
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png/group_single_before_conditional.png",2000,1000)

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

fdr5Data
p1e9Data

g1 = ggplot(data = fdr5Data, aes(x = annotation, y = values, fill = association)) +
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

g2 = ggplot(data = p1e9Data, aes(x = annotation, y = values, fill = association)) + 
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

g1 + g2

dev.off()