
library(data.table)
library(ggplot2)
library(patchwork)
options(scipen = 999)

png("~/git/noncoding_rvat/figures/revision/png/variant_gene_distance.png",2000,1300)


geneCoord = fread("~/git/noncoding_rvat/supporting_data/gencode_coordinates.tsv")
geneCoord$tss = geneCoord$V2
geneCoord[V4 == "-"]$tss = geneCoord[V4 == "-"]$V3
colnames(geneCoord) = c("chro","start","end","strand","gene","tss")

# bins <- c("0-1kb", "1-5kb", "5-20kb", "20-50kb", "50-100kb", "100-200kb", "200-500kb", "500-1000kb")
# mergedData$bin = cut(mergedData$dist, breaks = c(0, 1000, 5000, 20000, 50000, 100000, 200000, 500000, 1000000), labels = bins)
# mergedData[dist == 0]$bin = "0-1kb"

plotHist = function(annotData, tag){
  # annotData$chro = sapply(strsplit(as.character(annotData$V1), ":"), function(x) x[1])
  # annotData$var_pos = as.numeric(sapply(strsplit(as.character(annotData$V1), ":"), function(x) x[2]))
  colnames(annotData) = c("chro","var_pos","ref","alt","gene","annot")
  mergedData = merge(annotData, geneCoord, by = c("gene"))
  mergedData = mergedData[chro.x == chro.y]
  mergedData$dist = mergedData$var_pos - mergedData$tss
  median = median(abs(mergedData$dist))
  mean = mean(abs(mergedData$dist))
  text = paste("Abs. Median:",median,"\nAbs. Mean:",round(mean,1))
  
  plot =  ggplot(mergedData, aes(x = dist)) +
    geom_histogram( fill = "white", color = "black", bins = 30) +
    annotate(geom = "text", x = Inf, y = Inf, label = text, hjust = 1.2, vjust = 1.3, size = 8) +
    labs(x="Distance to TSS (kb)", y = "# variants", tag = tag) +
    scale_x_continuous(breaks = seq(-1000000, 1000000, by = 200000),
                       labels = paste0(seq(-1000, 1000, by = 200)), limits = c(-1000000,1000000)) +
    # xlim(c(-1000000,1000000)) +
    theme_linedraw() +
    theme(text = element_text(size=30), plot.title = element_text(hjust = 0.5),  axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=2), plot.tag.position = c(0.0, 1.03), plot.tag = element_text(size = rel(1.1), hjust = 0))
  
  return(plot)
  
}

annotData = fread("zcat ~/git/noncoding_rvat/revision/no_cond_42_traits/annotation/ABC.CADD15.annotations.gz", header = F)
g1 = plotHist(annotData, "a - ABC")

gc()

annotData = fread("zcat ~/git/noncoding_rvat/revision/no_cond_42_traits/annotation/CRD.CADD15.annotations.gz", header = F)
g2 = plotHist(annotData, "b - CRD")

gc()

annotData = fread("zcat ~/git/noncoding_rvat/revision/no_cond_42_traits/annotation/JAVIERRE.CADD15.annotations.gz", header = F)
g3 = plotHist(annotData, "c - HIC")

gc()

annotData = fread("zcat ~/git/noncoding_rvat/revision/no_cond_42_traits/annotation/NC.CADD15.annotations.gz", header = F)
g4 = plotHist(annotData, "d - NC")

gc()

annotData = fread("zcat ~/git/noncoding_rvat/revision/no_cond_42_traits/annotation/CDS.CADD15.annotations.gz", header = F)
g5 = plotHist(annotData, "e - CDS")

gc()

annotData = fread("zcat ~/git/noncoding_rvat/revision/no_cond_42_traits/annotation/PROMOTER.CADD15.annotations.gz", header = F)
g6 = plotHist(annotData, "f - PROMOTER")

gc()

g1 + g2 + g3 + g4 + g5 + g6

dev.off()
