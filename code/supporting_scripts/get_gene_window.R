options(scipen = 100)
library(data.table)
data = fread("gencode_coordinates.tsv")
wantedGenes = fread("autosomal_coding_genes.bed")
data = data[V5 %in% wantedGenes$V4]
data$tss = data$V2
data[V4 == "-"]$tss = data[V4=="-"]$V3
data$start = data$tss - 1000000
data$end = data$tss + 1000000
data[start < 0]$start = 0
data[end < 0]$end = 0
write.table(data[,.(V1,start,end,V5,V4)],"coding_gene_1MB_window.bed",quote=F,sep="\t",row.names=F,col.names=F)
