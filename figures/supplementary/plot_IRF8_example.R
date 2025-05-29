
library(data.table)
library(ggplot2)
scale <- function(x){(x-min(x))/(max(x)-min(x))}

png("~/git/noncoding_rvat/figures/revision/png/example_IRF8.png",1200,600)

wantedGene = "ENSG00000140968"
wantedChr = "chr16"

cdsAllData = fread("~/git/noncoding_rvat/revision/example/CDS_wanted_genes.bed.gz", header = F)
geneData = fread("~/git/noncoding_rvat/supporting_data/autosomal_coding_genes.bed")
abcData = fread("~/git/noncoding_rvat/revision/example/IRF8/ABC_IRF8", header = F)
crdData = fread("~/git/noncoding_rvat/revision/example/IRF8/CRD_IRF8", header = F)
hicData = fread("~/git/noncoding_rvat/revision/example/IRF8/HIC_IRF8", header = F)
ncData = fread("~/git/noncoding_rvat/revision/example/IRF8/NC_IRF8", header = F)
gwasData = fread("~/git/noncoding_rvat/revision/example/IRF8/panUKB.p30130.IRF8.hg38.tsv", header = F)
singleBeforeData = fread("~/git/noncoding_rvat/revision/example/IRF8/NC.single.before_cond.chr16.p30130.regenie", header = F)
singleAfterData = fread("~/git/noncoding_rvat/revision/example/IRF8/NC.single.after_cond.chr16.p30130.regenie", header = F)
consData = fread("~/git/noncoding_rvat/revision/example/IRF8/hg38.phastCons30way.IRF8.bed")

singleBeforeData = singleBeforeData[V2 %in% ncData$V2]
singleAfterData = singleAfterData[V2 %in% ncData$V2]

tssData = fread("~/git/noncoding_rvat/supporting_data/gencode_coordinates.tsv")
tssData = tssData[V5 == wantedGene]
tssData$tss = tssData$V2
tssData[V4 == "-"]$tss = tssData[V4 == "-"]$V3

# min = min(tssData$tss) - 1000000
# max = max(tssData$tss) + 1000000
min = min(tssData$tss) - 800000
max = max(tssData$tss) + 800000
# min = min(tssData$tss) - 10000
# max = max(tssData$tss) + 10000

cdsData = unique(cdsAllData[V4 == wantedGene])
cdsOtherData = unique(cdsAllData[V1 == unique(cdsData$V1)][V2 > min & V3 < max])
wantedGeneData = geneData[V4 == wantedGene]
otherGeneData = geneData[V1 == wantedChr & V2 > min & V3 < max]

# gwasData[V4 > 308]$V4 = 308
gwasData = gwasData[!is.na(V4)]
gwasData = rbind(gwasData, data.table("test",0,0,-log10(5e-8)))
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

summary(consData$V4)
consData[V4 > 0.9]

#####################
numTracks = 6
midFactor = 0.33

ggplot() +
  geom_rect(data = cdsOtherData, aes(xmin = V2, xmax = V3), ymin = numTracks + midFactor, ymax = numTracks + midFactor * 2, fill = "#d9d9d9", color = "#d9d9d9") +
  geom_rect(data = cdsData, aes(xmin = V2, xmax = V3), ymin = numTracks + midFactor, ymax = numTracks + midFactor * 2, fill = "black", color = "black") +
  geom_rect(data = otherGeneData, aes(xmin = V2, xmax = V3), ymin = numTracks + midFactor + 0.16, ymax = numTracks + midFactor + 0.17, fill = "#d9d9d9", color = "#d9d9d9") +
  geom_rect(data = wantedGeneData, aes(xmin = V2, xmax = V3), ymin = numTracks + midFactor + 0.16, ymax = numTracks + midFactor + 0.17, fill = "black", color = "black") +
  # geom_text(data = otherGeneData[V6 != "IRF8"], aes(x = V2, label = V6), y = numTracks + midFactor*3, size = 5, color = "#d9d9d9") + 
  annotate(geom = "text", label = "IRF8", x = tssData$tss, y = numTracks + midFactor*3, size = 5) + 
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
  geom_point(data = consData, aes(x = V2, y = V4 + numTracks - 6.3), color = "#238b45", size = 1, shape = 16) +
  annotate(geom = "text", label = "Genes", x = min-11000, y = numTracks + 0.7, hjust = 1, size = 5) + 
  annotate(geom = "text", label = "HIC", x = min-11000, y = numTracks-1 + midFactor*3 - 0.2 , hjust = 1, color = "#a50f15", size = 5) + 
  annotate(geom = "text", label = "CRD", x = min-11000, y = numTracks-1 + midFactor*2 - 0.2, hjust = 1, color = "#756bb1", size = 5) + 
  annotate(geom = "text", label = "ABC", x = min-11000, y = numTracks-1 + midFactor - 0.2, hjust = 1, color = "#006d2c", size = 5) + 
  annotate(geom = "text", label = "GWAS", x = min-11000, y = numTracks-2 + 0.2, hjust = 1, color = "black", size = 5) + 
  annotate(geom = "text", label = "Rare no\n cond.", x = min-11000, y = numTracks-3, hjust = 0.9, color = "black", size = 5) + 
  annotate(geom = "text", label = "Rare \ncond.", x = min-11000, y = numTracks-4-0.3, hjust = 1, color = "black", size = 5) + 
  annotate(geom = "text", label = "Conserv.", x = min-11000, y = numTracks-5-0.8, hjust = 0.85, color = "black", size = 5) + 
  geom_hline(yintercept = numTracks + 0.1, color = "grey") + 
  geom_hline(yintercept = numTracks-1 - 0.1, color = "grey") + 
  geom_hline(yintercept = numTracks-2 - midFactor - 0.1, color = "grey") + 
  geom_hline(yintercept = numTracks-3 - midFactor - 0.3, color = "grey") + 
  geom_hline(yintercept = numTracks-4 - midFactor - 0.7, color = "grey") + 
  xlim(c(min-17000,max-17000)) +
  ylim(c(-0.3,numTracks + 1)) +
  xlab("Chr16 position") +
  ylab("Track") + 
  theme_linedraw() +
  theme(axis.text.y = element_blank(), text = element_text(size=22), axis.ticks.y = element_blank(),
        panel.grid.minor.y =  element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5)
  )

dev.off()
