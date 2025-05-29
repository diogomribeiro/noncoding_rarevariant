
library(data.table)
# Counting number of associated genes independent of other genes in window 
# data = fread("/home/dribeiro/git/noncoding_rvat/revision/no_cond_42_traits/NC_SKATO_process_FDR0.05.tsv")
# data = fread("/home/dribeiro/git/noncoding_rvat/revision/no_cond_42_traits/CDS_SKATO_process_FDR0.05.tsv")
# data = fread("/home/dribeiro/git/noncoding_rvat/revision/no_cond_42_traits/NC_SKATO_process_PVAL1e9.tsv")
data = fread("/home/dribeiro/git/noncoding_rvat/revision/no_cond_42_traits/CDS_SKATO_process_PVAL1e9.tsv")

coords = fread("/home/dribeiro/git/noncoding_rvat/supporting_data/gencode_coordinates.tsv")
coords$tss = coords$V2
coords[V4 == "-"]$tss = coords[V4 == "-"]$V3

mergedData = merge(data, coords[,.(V5,tss)], by.x = "GENE", by.y = "V5")

length(unique(mergedData[order(mergedData$tss), ]$GENE))

window <- 500000
alreadySeen <- c()
indep <- 0
for (gene in unique(mergedData[order(mergedData$tss), ]$GENE)) {
  entry <- head(mergedData[GENE == gene, ], 1)
  
  # Check if the gene is not in alreadySeen
  if (!(entry$GENE %in% alreadySeen)) {
    indep <- indep + 1
    
    # Find other genes within the window
    otherGenes <- mergedData[CHROM == entry$CHROM & 
                               tss > entry$tss - window & 
                               tss < entry$tss + window, ]
    
    # Add unique genes to alreadySeen
    alreadySeen <- c(alreadySeen, unique(otherGenes$GENE))
  }
}
indep

indep/length(unique(mergedData$GENE))
