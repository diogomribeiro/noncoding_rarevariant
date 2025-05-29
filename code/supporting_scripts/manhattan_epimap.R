
library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)

# Create Manhattan plot
manhattanPlot = function(mergedData, color1 = "green", color2 = "blue", line2y = 5, tag = "", cap = 50, range_width = 1000000){
  mergedData = mergedData[order(chro,tss)]
  mergedData[LOG10P > 50]$LOG10P = cap
  mergedData$index = seq(nrow(mergedData))
  midindex = round(mergedData[, .(mean_value = mean(index)), by = chro],0)
  
  # select genes to plot label
  top_hits = subset(mergedData, LOG10P > 9)
  top_hits = top_hits[top_hits$LOG10P == ave(top_hits$LOG10P, top_hits$GENENAME, FUN = max), ]
  top_hits[, group := round(tss / range_width)]
  top_hits = top_hits[order(group, -LOG10P)]
  top_hits = top_hits[, .SD[1], by = group]
  
  gg = ggplot(mergedData, aes(x = index, y = LOG10P, color = factor(ifelse(chro %% 2 == 0, "Even", "Odd")))) +
    geom_point(size = 1) +
    geom_text_repel(data = top_hits, aes(label = ifelse(LOG10P > 9, sprintf(GENENAME))), vjust = -0.5, size = 5, fontface = "bold", max.overlaps = 20) +
    geom_hline(yintercept = 9, color = "black", linetype = "dashed") + 
    geom_hline(yintercept = line2y, color = "#2c7fb8", linetype = "dashed") + # TODO: change depending on annotation
    scale_color_manual(values = c("Even" = color1, "Odd" = color2)) +
    scale_x_continuous(breaks = midindex$mean_value, labels = midindex$chro) + 
    ylim(c(1, cap)) +
    labs(x="Chromosome", y = "-log10(p-value)", tag = tag) +
    theme_linedraw() +
    theme(text = element_text(size=22), plot.title = element_text(hjust = 0.5), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", plot.margin = margin(1, 0, 0, 0, "cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=2), plot.tag.position = c(0.0, 1.05), plot.tag = element_text(size = rel(1.2), hjust = 0))
  
  return(gg)
}

png("~/git/noncoding_rvat/revision/epimap/manhattan_epimap.png",1800,350*6)

### GENE COORDINATES
geneCoordinates = fread("~/git/noncoding_rvat/supporting_data/gencode_coordinates.tsv")
geneCoordinates$tss = geneCoordinates$V2
geneCoordinates[V4 == "-"]$tss = geneCoordinates[V4 == "-"]$V3

################
# T cell
################
mergedData1 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/EPI_BLOOD_T_CELL_SKATO_process_all.tsv.gz")
mergedData2 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/conditional/EPI_BLOOD_T_CELL_SKATO_conditional_all.tsv")
mergedData1 = mergedData1[PVAL > 1e-9]
mergedData1$FDR = NULL
mergedData = rbind(mergedData1, mergedData2)
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
g1 = manhattanPlot(mergedData, "#fec44f", "#fff7bc", 0, "a - T cell") # blue

################
# B cell
################
mergedData1 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/EPI_HSC_B_CELL_SKATO_process_all.tsv.gz")
mergedData2 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/conditional/EPI_HSC_B_CELL_SKATO_conditional_all.tsv")
mergedData1 = mergedData1[PVAL > 1e-9]
mergedData1$FDR = NULL
mergedData = rbind(mergedData1, mergedData2)
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
g2 = manhattanPlot(mergedData, "#bae4b3", "#006d2c", 0, "b - B cell") # blue


################
# Liver
################
mergedData1 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/EPI_LIVER_SKATO_process_all.tsv.gz")
mergedData2 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/conditional/EPI_LIVER_SKATO_conditional_all.tsv")
mergedData1 = mergedData1[PVAL > 1e-9]
mergedData1$FDR = NULL
mergedData = rbind(mergedData1, mergedData2)
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
g3 = manhattanPlot(mergedData, "#bcbddc", "#756bb1", 0, "c - Liver") # blue


################
# Kidney
################
mergedData1 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/EPI_KIDNEY_SKATO_process_all.tsv.gz")
mergedData2 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/conditional/EPI_KIDNEY_SKATO_conditional_all.tsv")
mergedData1 = mergedData1[PVAL > 1e-9]
mergedData1$FDR = NULL
mergedData = rbind(mergedData1, mergedData2)
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
g4 = manhattanPlot(mergedData, "#fcae91", "#a50f15", 0, "d - Kidney") # blue


################
# Pancreas
################
mergedData1 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/EPI_PANCREAS_SKATO_process_all.tsv.gz")
mergedData2 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/conditional/EPI_PANCREAS_SKATO_conditional_all.tsv")
mergedData1 = mergedData1[PVAL > 1e-9]
mergedData1$FDR = NULL
mergedData = rbind(mergedData1, mergedData2)
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
g5 = manhattanPlot(mergedData, "#9ecae1", "#2171b5", 0, "e - Pancreas") # blue


################
# Endo
################
mergedData1 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/EPI_ENDO_SKATO_process_all.tsv.gz")
mergedData2 = fread("/home/dribeiro/git/noncoding_rvat/revision/epimap/conditional/EPI_ENDO_SKATO_conditional_all.tsv")
mergedData1 = mergedData1[PVAL > 1e-9]
mergedData1$FDR = NULL
mergedData = rbind(mergedData1, mergedData2)
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
g6 = manhattanPlot(mergedData, "#d9d9d9", "#636363", 0, "f - Endothelial") # blue


g1 / g2 / g3 / g4 / g5 / g6
dev.off()
