
library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)

# png("~/git/noncoding_rvat/figures/final/figure2.png",1800,2450)
pdf("~/git/noncoding_rvat/figures/final/figure2.pdf", width = 25, height = 34.03)

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


### GENE COORDINATES
geneCoordinates = fread("~/git/noncoding_rvat/supporting_data/gencode_coordinates.tsv")
geneCoordinates$tss = geneCoordinates$V2
geneCoordinates[V4 == "-"]$tss = geneCoordinates[V4 == "-"]$V3

sourceData = data.table()

################
# ABC
################
mergedData = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/ABC_SKATO_process_all.tsv.gz")
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
p1 = mergedData[PVAL > 0.1][sample(nrow(mergedData)*0.1),]
p2 = mergedData[PVAL < 0.1]
mergedData = rbind(p1,p2)
mergedData$annotation = "ABC"
sourceData = rbind(sourceData, mergedData)
g1 = manhattanPlot(mergedData, "#74c476", "#006d2c", -log10(6.93E-06), "a - ABC") # green

################
# CRD
################
mergedData = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/CRD_SKATO_process_all.tsv.gz")
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
p1 = mergedData[PVAL > 0.1][sample(nrow(mergedData)*0.05),]
p2 = mergedData[PVAL < 0.1]
mergedData = rbind(p1,p2)
mergedData$annotation = "CRD"
sourceData = rbind(sourceData, mergedData)
g2 = manhattanPlot(mergedData, "#9e9ac8", "#54278f", -log10(3.68E-05), "b - CRD") # purple

################
# HIC
################
mergedData = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/JAVIERRE_SKATO_process_all.tsv.gz")
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
p1 = mergedData[PVAL > 0.1][sample(nrow(mergedData)*0.05),]
p2 = mergedData[PVAL < 0.1]
mergedData = rbind(p1,p2)
mergedData$annotation = "HIC"
sourceData = rbind(sourceData, mergedData)
g3 = manhattanPlot(mergedData, "#fb6a4a", "#a50f15", -log10(5.49E-05), "c - HIC") # red

################
# NC
################
mergedData = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/NC_SKATO_process_all.tsv.gz")
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
p1 = mergedData[PVAL > 0.1][sample(nrow(mergedData)*0.05),]
p2 = mergedData[PVAL < 0.1]
mergedData = rbind(p1,p2)
mergedData$annotation = "NC"
sourceData = rbind(sourceData, mergedData)
g4 = manhattanPlot(mergedData, "#6baed6", "#08519c", -log10(6.50E-05), "d - NC") # blue

################
# NC - conditional
################
mergedData = fread("~/git/noncoding_rvat/revision/cond_rerun/NC_vcmax_SKATO_conditional_all.tsv.gz")
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
p1 = mergedData[PVAL > 0.1][sample(nrow(mergedData)*0.05),]
p2 = mergedData[PVAL < 0.1]
mergedData = rbind(p1,p2)
mergedData$annotation = "NC_COND"
sourceData = rbind(sourceData, mergedData)
g5 = manhattanPlot(mergedData, "#6baed6", "#08519c", -log10(4.56e-06), "e - NC (cond)") # blue

################
# CDS
################
mergedData = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/CDS_SKATO_process_all.tsv.gz")
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
p1 = mergedData[PVAL > 0.1][sample(nrow(mergedData)*0.05),]
p2 = mergedData[PVAL < 0.1]
mergedData = rbind(p1,p2)
mergedData$annotation = "CDS"
sourceData = rbind(sourceData, mergedData)
g6 = manhattanPlot(mergedData, "#969696", "#252525", -log10(4.60E-05), "f - CDS") # grey

################
# CDS - conditional
################
mergedData = fread("~/git/noncoding_rvat/revision/cond_rerun/CDS_vcmax_SKATO_conditional_all.tsv.gz")
mergedData = merge(mergedData, geneCoordinates, by.x = "GENE", by.y = "V5")
mergedData$chro = mergedData$CHROM
p1 = mergedData[PVAL > 0.1][sample(nrow(mergedData)*0.05),]
p2 = mergedData[PVAL < 0.1]
mergedData = rbind(p1,p2)
mergedData$annotation = "CDS_COND"
sourceData = rbind(sourceData, mergedData)
g7 = manhattanPlot(mergedData, "#969696", "#252525", -log10(1.91e-05), "g - CDS (cond)") # grey


# write.table(sourceData[,.(TRAIT,CHROM,GENPOS,LOG10P,GENENAME,annotation)],"/home/dribeiro/git/noncoding_rvat/figures/final/source_data/Ribeiro_SourceData_Fig2.tsv",quote=F,sep="\t",row.names=F)


g1 / g2 / g3 / g4 / g5 / g6 / g7
# dev.off()


## Source data
library(openxlsx)

wb <- createWorkbook()
panels <- c("ABC", "CRD", "HIC", "NC", "NC_COND", "CDS", "CDS_COND")
sheet_names <- paste0("fig2", letters[1:length(panels)])

# For each panel, add a worksheet and write the data
for (i in seq_along(panels)) {
  panel <- panels[i]
  sheet <- sheet_names[i]
  panel_data <- sourceData[annotation == panel, .(TRAIT, CHROM, GENPOS, LOG10P, GENENAME, annotation)]
  addWorksheet(wb, sheetName = sheet)
  writeData(wb, sheet = sheet, panel_data)
}

saveWorkbook(wb, file = "~/git/noncoding_rvat/figures/final/source_data/Ribeiro_SourceData_Fig2.xlsx", overwrite = TRUE)
