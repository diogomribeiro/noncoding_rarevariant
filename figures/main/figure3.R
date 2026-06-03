
library(data.table)
library(ggplot2)
library(patchwork)

# png("~/git/noncoding_rvat/figures/final/figure3.png",1800,1800)
pdf("~/git/noncoding_rvat/figures/final/figure3.pdf",25,25)

################################
# Panel A
################################

scale <- function(x){(x-min(x))/(max(x)-min(x))}

wantedGene = "ENSG00000101439"
wantedChr = "chr20"

cdsAllData = fread("~/git/noncoding_rvat/revision/example/CDS_wanted_genes.bed.gz", header = F)
geneData = fread("~/git/noncoding_rvat/supporting_data/autosomal_coding_genes.bed")
abcData = fread("~/git/noncoding_rvat/revision/example/CST3/ABC_CST3", header = F)
crdData = fread("~/git/noncoding_rvat/revision/example/CST3/CRD_CST3", header = F)
hicData = fread("~/git/noncoding_rvat/revision/example/CST3/HIC_CST3", header = F)
ncData = fread("~/git/noncoding_rvat/revision/example/CST3/NC_CST3", header = F)
gwasData = fread("~/git/noncoding_rvat/revision/example/CST3/panUKB.p30720.CST3.hg38.bed", header = F)
singleBeforeData = fread("~/git/noncoding_rvat/revision/example/CST3/NC.single.before_cond.chr20.p30720.regenie", header = F)
singleAfterData = fread("~/git/noncoding_rvat/revision/example/CST3/NC.single.after_cond.chr20.p30720.regenie", header = F)
consData = fread("~/git/noncoding_rvat/revision/example/CST3/hg38.phastCons30way.CST3.bed.gz")

singleBeforeData = singleBeforeData[V2 %in% ncData$V2]
singleAfterData = singleAfterData[V2 %in% ncData$V2]

tssData = fread("~/git/noncoding_rvat/supporting_data/gencode_coordinates.tsv")
tssData = tssData[V5 == wantedGene]
tssData$tss = tssData$V2
tssData[V4 == "-"]$tss = tssData[V4 == "-"]$V3

# min = min(tssData$tss) - 1000000
# max = max(tssData$tss) + 1000000
min = min(tssData$tss) - 250000
max = max(tssData$tss) + 250000
# min = min(tssData$tss) - 25000
# max = max(tssData$tss) + 25000

cdsData = unique(cdsAllData[V4 == wantedGene])
cdsOtherData = unique(cdsAllData[V1 == unique(cdsData$V1)][V2 > min & V3 < max])
wantedGeneData = geneData[V4 == wantedGene]
otherGeneData = geneData[V1 == wantedChr & V2 > min & V3 < max]

gwasData = gwasData[!is.na(V4)]
gwasData = rbind(gwasData, data.table("test",0,0,-log10(1e-9)))
gwasData$scaled = scale(gwasData$V4)
gwasDataThreshold = gwasData[V1 == "test"]$scaled
gwasData = gwasData[V1 != "test"]

singleBeforeData = rbind(singleBeforeData, data.table("test",0,"test","test",0,-log10(1e-9),"test"))
singleBeforeData$scaled = scale(singleBeforeData$V6)
singleBeforeDataThreshold = singleBeforeData[V1 == "test"]$scaled
singleBeforeData = singleBeforeData[V1 != "test"]

singleAfterData = rbind(singleAfterData, data.table("test",0,"test","test",-log10(1e-9)))
singleAfterData$scaled = scale(singleAfterData$V5)
singleAfterDataThreshold = singleAfterData[V1 == "test"]$scaled
singleAfterData = singleAfterData[V1 != "test"]

consData = consData[V2 > min & V3 < max]
consData$scaled = scale(consData$V4)

#####################
numTracks = 6
midFactor = 0.33

g1 = ggplot() +
  geom_rect(data = cdsOtherData, aes(xmin = V2, xmax = V3), ymin = numTracks + midFactor, ymax = numTracks + midFactor * 2, fill = "#969696", color = "#969696") +
  geom_rect(data = cdsData, aes(xmin = V2, xmax = V3), ymin = numTracks + midFactor, ymax = numTracks + midFactor * 2, fill = "black", color = "black") +
  geom_rect(data = otherGeneData, aes(xmin = V2, xmax = V3), ymin = numTracks + midFactor + 0.16, ymax = numTracks + midFactor + 0.17, fill = "#969696", color = "#969696") +
  geom_rect(data = wantedGeneData, aes(xmin = V2, xmax = V3), ymin = numTracks + midFactor + 0.16, ymax = numTracks + midFactor + 0.17, fill = "black", color = "black") +
  geom_text(data = otherGeneData[V6 != "CST3"], aes(x = V2, label = V6), y = numTracks + midFactor*3, size = 7, color = "#969696") + 
  annotate(geom = "text", label = "CST3", x = tssData$tss, y = numTracks + midFactor*3, size = 7) + 
  geom_point(data = tssData, aes(x = tss), y = numTracks + 0.5, fill = "black", color = "black", shape = 9) +
  geom_rect(data = abcData, aes(xmin = V2, xmax = V3), ymin = numTracks-1, ymax = numTracks-1 + midFactor*2, fill = "#006d2c") +
  geom_rect(data = crdData, aes(xmin = V2, xmax = V3), ymin = numTracks-1 + midFactor, ymax = numTracks-1 + midFactor*2, fill = "#756bb1") +
  geom_rect(data = hicData, aes(xmin = V2, xmax = V3), ymin = numTracks-1 + midFactor*2, ymax = numTracks-1 + midFactor*3, fill = "#a50f15") +
  geom_point(data = gwasData, aes(x = V2, y = scaled + numTracks - 3 + midFactor*2), fill = "black", color = "black", shape = 1, size = 1) +
  geom_hline(yintercept = gwasDataThreshold + numTracks - 3 + midFactor*2 + 0.1, color = "black", linetype = "dashed") +
  geom_point(data = singleBeforeData, aes(x = V2, y = scaled + numTracks - 4 + midFactor*2 - 0.2), fill = "#feb24c", color = "#feb24c", shape = 1) +
  geom_hline(yintercept = singleBeforeDataThreshold + numTracks - 4 + midFactor*2 - 0.2, color = "#feb24c", linetype = "dashed") +
  geom_point(data = singleAfterData, aes(x = V2, y = scaled + numTracks - 5 + midFactor - 0.2), fill = "#41b6c4", color = "#41b6c4", shape = 1) +
  geom_hline(yintercept = singleAfterDataThreshold + numTracks - 5 + midFactor - 0.2, color = "#41b6c4", linetype = "dashed") +
  geom_point(data = consData, aes(x = V2, y = V4 + numTracks - 6.3), color = "#41ab5d", size = 1, shape = 16, alpha = 0.5) +
  annotate(geom = "text", label = "Genes", x = min-11000, y = numTracks + 0.7, hjust = 1, size = 7) + 
  annotate(geom = "text", label = "HIC", x = min-11000, y = numTracks-1 + midFactor*3 - 0.2 , hjust = 1, color = "#a50f15", size = 7) + 
  annotate(geom = "text", label = "CRD", x = min-11000, y = numTracks-1 + midFactor*2 - 0.2, hjust = 1, color = "#756bb1", size = 7) + 
  annotate(geom = "text", label = "ABC", x = min-11000, y = numTracks-1 + midFactor - 0.2, hjust = 1, color = "#006d2c", size = 7) + 
  annotate(geom = "text", label = "GWAS", x = min-11000, y = numTracks-2 + 0.2, hjust = 1, color = "black", size = 7) + 
  annotate(geom = "text", label = "Rare no\n cond.", x = min-11000, y = numTracks-3, hjust = 0.9, color = "black", size = 7) + 
  annotate(geom = "text", label = "Rare \ncond.", x = min-11000, y = numTracks-4-0.3, hjust = 1, color = "black", size = 7) + 
  annotate(geom = "text", label = "Conserv.", x = min-11000, y = numTracks-5-0.8, hjust = 0.85, color = "black", size = 7) + 
  geom_hline(yintercept = numTracks + 0.1, color = "grey") + 
  geom_hline(yintercept = numTracks-1 - 0.1, color = "grey") + 
  geom_hline(yintercept = numTracks-2 - midFactor - 0.1, color = "grey") + 
  geom_hline(yintercept = numTracks-3 - midFactor - 0.3, color = "grey") + 
  geom_hline(yintercept = numTracks-4 - midFactor - 0.7, color = "grey") + 
  labs(tag = "a") +
  xlim(c(min-17000,max-17000)) +
  ylim(c(-0.3,numTracks + 1)) +
  xlab("Chr20 position") +
  ylab("Track") + 
  theme_linedraw() +
  theme(axis.text.y = element_blank(), text = element_text(size=36), axis.ticks.y = element_blank(),
        panel.grid.minor.y =  element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5)
  )


################################
# Panel B
################################

nAssoc = fread("~/git/noncoding_rvat/revision/known_replication/sample_size_gene_trait_assoc.txt")

data = fread("~/git/noncoding_rvat/revision/known_replication/full.summary")
data$pval_cutoff = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[2])))$V1
data$annotation = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
data[annotation == "NC_shuffle"]$annotation = "EXPEC."
data[annotation == "PROMOTER"]$annotation = "PROM."
data[annotation == "JAVIERRE"]$annotation = "HIC"
colnames(data) = c("replication","percentage","V3","pval_cutoff","annotation")

data = merge(data, nAssoc, by = "annotation")


# Direct replication
g2 = ggplot(data[replication == "direct_replication"][pval_cutoff == "FDR0"], aes(x=annotation, y=percentage, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", label = paste("N =",data$FDR5), size = 8, fontface= "bold", x = data$annotation, y = -3, vjust = 0.5) +
  scale_x_discrete(limits = c("ABC","CRD","HIC","NC","CDS","PROM.","CONTROL","RANDOM","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "b") +
  guides(fill = "none") +
  ylim(c(-3,105)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#02818a","#d9d9d9","black","#fed976"), limits = c("ABC","CRD","HIC","NC","CDS","PROM.","CONTROL","RANDOM","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=36), legend.position = c(0.40,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))

################################
# Panel C
################################

# Window replication
g3 = ggplot(data[replication == "window_replication"][pval_cutoff == "FDR0"], aes(x=annotation, y=percentage, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", label = paste("N =",data$FDR5), size = 8, fontface= "bold", x = data$annotation, y = -3, vjust = 0.5) +
  scale_x_discrete(limits = c("ABC","CRD","HIC","NC","CDS","PROM.","CONTROL","RANDOM","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "c") +
  guides(fill = "none") +
  ylim(c(-3,105)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#02818a","#d9d9d9","black","#fed976"), limits = c("ABC","CRD","HIC","NC","CDS","PROM.","CONTROL","RANDOM","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=36), legend.position = c(0.40,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))

g1 / (g2 + g3) + plot_layout(heights = c(1.5, 1))

dev.off()



################################
# Source data
################################
# Merge regulatory regions (ABC, CRD, HIC)
panelA_regions <- rbind(
  data.table(type = "ABC", abcData),
  data.table(type = "CRD", crdData),
  data.table(type = "HIC", hicData)
)

# Merge GWAS and rare variant points
panelA_points <- rbind(
  data.table(type = "GWAS", gwasData[, .(V2, scaled)]),
  data.table(type = "Rare_no_cond", singleBeforeData[, .(V2, scaled)]),
  data.table(type = "Rare_cond", singleAfterData[, .(V2, scaled)])
)

# Conservation data
panelA_conservation <- consData[, .(V2, V4, scaled)]
setnames(panelA_conservation, c("position", "score", "scaled_score"))

## Source data
library(openxlsx)
write.xlsx(list(fig3a_regions = panelA_regions, fig3a_assoc = panelA_points, fig3a_conservation = panelA_conservation, 
                fig3b = data[replication == "direct_replication"][pval_cutoff == "FDR0"], fig3c = data[replication == "window_replication"][pval_cutoff == "FDR0"]), 
           "~git/noncoding_rvat/figures/final/source_data/Ribeiro_SourceData_Fig3.xlsx")
