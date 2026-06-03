
library(data.table)
library(UpSetR)
# Note: patchwork does not work with these plots

png("~/git/noncoding_rvat/figures/revision/png/epimap_upset_plot_before.png",1600,900)

## Before conditional
bloodTCell = fread("~/git/noncoding_rvat/revision/epimap/EPI_BLOOD_T_CELL_SKATO_process_PVAL1e9.tsv")
hscBCell = fread("~/git/noncoding_rvat/revision/epimap/EPI_HSC_B_CELL_SKATO_process_PVAL1e9.tsv")
endo = fread("~/git/noncoding_rvat/revision/epimap/EPI_ENDO_SKATO_process_PVAL1e9.tsv")
kidney = fread("~/git/noncoding_rvat/revision/epimap/EPI_KIDNEY_SKATO_process_PVAL1e9.tsv")
liver = fread("~/git/noncoding_rvat/revision/epimap/EPI_LIVER_SKATO_process_PVAL1e9.tsv")
pancreas = fread("~/git/noncoding_rvat/revision/epimap/EPI_PANCREAS_SKATO_process_PVAL1e9.tsv")
nc = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/NC_SKATO_process_PVAL1e9.tsv")
cds = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/CDS_SKATO_process_PVAL1e9.tsv")

bloodTCell$cell_type = "T cells"
hscBCell$cell_type = "B cells"
endo$cell_type = "Endothelial"
kidney$cell_type = "Kidney"
liver$cell_type = "Liver"
pancreas$cell_type = "Pancreas"
nc$cell_type = "NC"
cds$cell_type = "CDS"

bloodTCell$tag = paste(bloodTCell$GENE, bloodTCell$TRAIT)
hscBCell$tag = paste(hscBCell$GENE, hscBCell$TRAIT)
endo$tag = paste(endo$GENE, endo$TRAIT)
kidney$tag = paste(kidney$GENE, kidney$TRAIT)
liver$tag = paste(liver$GENE, liver$TRAIT)
pancreas$tag = paste(pancreas$GENE, pancreas$TRAIT)
nc$tag = paste(nc$GENE, nc$TRAIT)
cds$tag = paste(cds$GENE, cds$TRAIT)

dataset = rbind(bloodTCell,hscBCell,endo,kidney,liver,pancreas,nc,cds)
binary_matrix <- dcast(dataset, tag ~ cell_type, fun.aggregate = length)
upset_data <- as.data.frame(binary_matrix)
upset_data$tag=NULL

# Create the UpSet plot
upset(
  upset_data,
  sets = colnames(upset_data),
  text.scale = 4)

dev.off()


## After conditional
png("~/git/noncoding_rvat/figures/revision/png/epimap_upset_plot_after.png",1600,900)

bloodTCell = fread("~/git/noncoding_rvat/revision/epimap/conditional/EPI_BLOOD_T_CELL_SKATO_conditional_PVAL1e9.tsv")
hscBCell = fread("~/git/noncoding_rvat/revision/epimap/conditional/EPI_HSC_B_CELL_SKATO_conditional_PVAL1e9.tsv")
endo = fread("~/git/noncoding_rvat/revision/epimap/conditional/EPI_ENDO_SKATO_conditional_PVAL1e9.tsv")
kidney = fread("~/git/noncoding_rvat/revision/epimap/conditional/EPI_KIDNEY_SKATO_conditional_PVAL1e9.tsv")
liver = fread("~/git/noncoding_rvat/revision/epimap/conditional/EPI_LIVER_SKATO_conditional_PVAL1e9.tsv")
pancreas = fread("~/git/noncoding_rvat/revision/epimap/conditional/EPI_PANCREAS_SKATO_conditional_PVAL1e9.tsv")
nc = fread("~/git/noncoding_rvat/revision/cond_rerun/NC_PVAL1e9_gene_trait_remaining.txt", header = F)
cds = fread("~/git/noncoding_rvat/revision/cond_rerun/CDS_PVAL1e9_gene_trait_remaining.txt", header = T)

bloodTCell$cell_type = "T cells"
hscBCell$cell_type = "B cells"
endo$cell_type = "Endothelial"
kidney$cell_type = "Kidney"
liver$cell_type = "Liver"
pancreas$cell_type = "Pancreas"
nc$cell_type = "NC"
cds$cell_type = "CDS"

bloodTCell$tag = paste(bloodTCell$GENE, bloodTCell$TRAIT)
hscBCell$tag = paste(hscBCell$GENE, hscBCell$TRAIT)
endo$tag = paste(endo$GENE, endo$TRAIT)
kidney$tag = paste(kidney$GENE, kidney$TRAIT)
liver$tag = paste(liver$GENE, liver$TRAIT)
pancreas$tag = paste(pancreas$GENE, pancreas$TRAIT)
nc$tag = paste(nc$V1, nc$V2)
cds$tag = paste(cds$GENE, cds$TRAIT)

dataset = rbind(bloodTCell[,.(tag,cell_type)],hscBCell[,.(tag,cell_type)],endo[,.(tag,cell_type)],kidney[,.(tag,cell_type)],liver[,.(tag,cell_type)],pancreas[,.(tag,cell_type)],nc[,.(tag,cell_type)],cds[,.(tag,cell_type)])
binary_matrix <- dcast(dataset, tag ~ cell_type, fun.aggregate = length)
upset_data <- as.data.frame(binary_matrix)
upset_data$tag=NULL

# Create the UpSet plot
upset(
  upset_data,
  sets = colnames(upset_data),
  text.scale = 4)

dev.off()