
library(data.table)

args = commandArgs(trailingOnly=TRUE)

# Gene info
wantedGenes = fread("wanted_genes_coding_autosomal_noMHC.txt", header = F)
coordData = fread("gencode_coordinates.tsv", header = F)
coordData$tss = coordData$V2
coordData[V4 == "-"]$tss = coordData[V4 == "-"]$V3
coordData = coordData[V5 %in% wantedGenes$V1]
window = 500000

## Process replication data
# repData = fread("data/genebass_75traits_significant_lof_missense.tsv")
# repData = fread("data/gwas_catalog_all_gene_75trait.tsv.gz")
repData = fread("data/genebass_gwas_catalog_75traits.tsv.gz")
repData = repData[gene_id %in% wantedGenes$V1]
repData$trait_id = data.table(unlist(lapply(repData$trait_id, function(x) unlist(strsplit(x,"[_]"))[1])))$V1

## Process burden results
burdenData = fread(args[1], header = T)

###########
burdenUniq = unique(burdenData[,.(GENE,TRAIT)])
paste("Processing entries:",nrow(burdenUniq))
finalDT = data.table()
for (tra in unique(burdenUniq$TRAIT)){
  genes = burdenUniq[TRAIT == tra]
  repGenes = repData[trait_id == tra]
  overl = nrow(genes[GENE %in% repGenes$gene_id])
  
  # direct overlap
  tp = overl
  fp = nrow(genes) - tp
  fn = nrow(repGenes) - tp
  tn = nrow(wantedGenes) - (tp + fp + fn)
  
  m = matrix(c(tp,fp,fn,tn), nrow = 2)
  fish = fisher.test(m)
  
  # overlap with expanded genes
  expandedGenes = c(repGenes$gene_id)
  for (g in unique(repGenes$gene_id)){
    g1 = coordData[V5 == g]
    newGenes = coordData[V1 == g1$V1][V2 > g1$tss - window][V2 < g1$tss + window]
    expandedGenes = c(expandedGenes, newGenes$V5)
  }
  expandedGenes = unique(expandedGenes)
  
  overlExp = nrow(genes[GENE %in% expandedGenes])
  # tp = overlExp
  # fp = nrow(genes) - tp
  # fn = length(expandedGenes) - tp
  # tn = nrow(wantedGenes) - (tp + fp + fn)
  # 
  # m2 = matrix(c(tp,fp,fn,tn), nrow = 2)
  # fish2 = fisher.test(m2)
  
  finalDT = rbind(finalDT, data.table(trait = tra, n_burden_genes = nrow(genes), n_rep_genes = nrow(repGenes), 
                                      overlap = overl, perc = overl * 100 / nrow(genes), 
                                      fisher_or = fish$estimate, fisher_pval = fish$p.value,
                                      n_expanded = length(expandedGenes),overlap_expanded = overlExp, perc_expanded = overlExp * 100 / nrow(genes))
                                      #fisher_or_expanded = fish2$estimate, fisher_pval_expanded = fish2$p.value)
  )
}

# Overall overlap
sum(finalDT$overlap) * 100 / sum(finalDT$n_burden_genes)
sum(finalDT$overlap_expanded) * 100 / sum(finalDT$n_burden_genes)

r = data.table(V1 = c("direct_replication","window_replication"),V2 = c(sum(finalDT$overlap) * 100 / sum(finalDT$n_burden_genes),sum(finalDT$overlap_expanded) * 100 / sum(finalDT$n_burden_genes)))

write.table(finalDT[order(n_burden_genes)],args[2], quote=F,sep="\t",row.names=F)
write.table(r,paste0(args[2],".summary"), quote=F,sep="\t",row.names=F, col.names=F)
