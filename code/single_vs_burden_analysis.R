
library(data.table)

singleData = fread("zcat ~/burdenNC/data/results/single/NC.single.processed.regenie.gz")
# singleData = fread("zcat zcat ~/burdenNC/data/results/single//CDS.single.processed.regenie.gz")
singleData = singleData[AF < 0.01 | AF > 0.99] # Note: single variant analysis was run regardless of MAF, whereas burden only considered MAF < 1%

# burdenData = fread("~/burdenNC/data/results/200k/NC_SKATO_process_FDR0.05.tsv")
burdenData = fread("~/burdenNC/data/results/200k/NC_SKATO_process_PVAL1e9.tsv")
# burdenData = fread("~/burdenNC/data/results/200k/CDS_SKATO_process_FDR0.05.tsv")
# burdenData = fread("~/burdenNC/data/results/200k/CDS_SKATO_process_PVAL1e9.tsv")

mapping = fread("zcat ~/burdenNC/data/annotations/NC.CADD15.annotations.gz", header = F)
# mapping = fread("zcat ~/burdenNC/data/annotations/CDS.CADD15.annotations.gz", header = F)

burdenData$LOG10P = -log10(burdenData$PVAL)

singleCount = 0
burdenCount = 0
singleRes = data.table()
burdenRes = data.table()
for (trait in unique(burdenData$TRAIT)){
  print(trait)
  traitBurdenData = burdenData[TRAIT == trait]
  traitSingleData = singleData[TRAIT == trait]
  
  for (gene in unique(traitBurdenData$GENE)){
    burden = traitBurdenData[GENE == gene]
    if (nrow(burden) > 1){print("PROBLEM: multiple gene-trait entries", burden)}
    burdenPval = burden$LOG10P
    
    map = mapping[V2 == gene]
    single = traitSingleData[VAR %in% map$V1]
    if (nrow(single) == 0){print("PROBLEM: no variants for gene")}
    singleCases = single[LOG10P >= burdenPval]
    if (nrow(singleCases) > 0){
      singleCount = singleCount + 1
      singleCases$GENE = gene
      singleCases$TRAIT = trait
      singleCases$SKATO = burdenPval
      singleRes = rbind(singleRes, data.table(singleCases))
    } else{ 
      burdenCount = burdenCount + 1
      burdenCases = single
      burdenCases$GENE = gene
      burdenCases$TRAIT = trait
      burdenCases$SKATO = burdenPval
      burdenRes = rbind(burdenRes, data.table(burdenCases))
    }
    
  }
  
}

singleCount # number of gene-trait associations where single variant has stronger or equal association as SKAT-O test
burdenCount # number of gene-trait associations where single variant has weaker association as SKAT-O test


sing = unique(singleRes[,.(TRAIT, GENE)])
burd = unique(burdenRes[,.(GENE,TRAIT)])
sing$CLASS = "Single"
burd$CLASS = "Group"
mergedData = rbind(sing,burd)
mergedData$tag = paste0(mergedData$GENE,"|",mergedData$TRAIT)

# number of 'burden cases' where at least 1 single variant passes genome-wide significance
cutoff = -log10(5e-8)
singleGroup = unique(burdenRes[LOG10P > cutoff][,.(GENE,TRAIT)])
singleGroup$tag = paste0(singleGroup$GENE,"|",singleGroup$TRAIT)
mergedData[tag %in% singleGroup$tag]$CLASS = "Group-Single"
mergedData$tag = NULL


traitMapping = fread("~/git/burdenNC/src/regenie/analysis/wanted_traits.txt", header = F)
traitMapping$TRAIT = data.table(unlist(lapply(traitMapping$V1, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
geneMapping = fread("~/git/burdenNC/src/saige/annotations/gencode_id_mapping.tsv")

mergedData = merge(mergedData, traitMapping[,.(TRAIT,V2)], by = "TRAIT")
mergedData = merge(mergedData, geneMapping[,.(V4,V6)], by.x = "GENE", by.y = "V4")
colnames(mergedData)[4] = "TRAIT_NAME"
colnames(mergedData)[5] = "GENE_NAME"

# write.table(mergedData,"~/git/burdenNC/src/regenie/analysis/single/NC.classification.FDR0.05.tsv",quote=F,row.names=F,sep="\t")
# write.table(singleRes,"dnanexus/NC.single.cases.FDR0.05.tsv",quote=F,row.names=F,sep="\t")
# write.table(burdenRes[LOG10P > 3],"dnanexus/NC.burden.cases.FDR0.05.tsv",quote=F,row.names=F,sep="\t")
write.table(mergedData,"~/git/burdenNC/src/regenie/analysis/single/NC.classification.P1e9.tsv",quote=F,row.names=F,sep="\t")
# write.table(singleRes,"dnanexus/NC.single.cases.P1e9.tsv",quote=F,row.names=F,sep="\t")
# write.table(burdenRes[LOG10P > 3],"dnanexus/NC.burden.cases.P1e9.tsv",quote=F,row.names=F,sep="\t")

# write.table(mergedData,"~/git/burdenNC/src/regenie/analysis/single/CDS.classification.FDR0.05.tsv",quote=F,row.names=F,sep="\t")
# write.table(singleRes,"dnanexus/CDS.single.cases.FDR0.05.tsv",quote=F,row.names=F,sep="\t")
# write.table(burdenRes[LOG10P > 3],"dnanexus/CDS.burden.cases.FDR0.05.tsv",quote=F,row.names=F,sep="\t")
# write.table(mergedData,"/home/dribeiro/git/burdenNC/src/regenie/analysis/single/CDS.classification.P1e9.tsv",quote=F,row.names=F,sep="\t")
# write.table(singleRes,"dnanexus/CDS.single.cases.P1e9.tsv",quote=F,row.names=F,sep="\t")
# write.table(burdenRes[LOG10P > 3],"dnanexus/CDS.burden.cases.P1e9.tsv",quote=F,row.names=F,sep="\t")

# check AC of associated single variants (could be interesting to report)
summary(singleRes$AF)
summary(singleData$AF)
cor.test(singleRes$LOG10P, singleRes$SKATO)


# # # CODE USED TO PROCESS RAW SINGLE FILE
# singleData = fread("zcat dnanexus/single/CDS.single.wanted.regenie.gz")
# colnames(singleData) = c("VAR","AF","LOG10P","TAG")
# singleData$TRAIT = data.table(unlist(lapply(singleData$TAG, function(x) unlist(strsplit(x,"[.]"))[2])))$V1
# singleData$TRAIT = data.table(unlist(lapply(singleData$TRAIT, function(x) unlist(strsplit(x,"[_]"))[2])))$V1
# write.table(singleData[,.(VAR,AF,LOG10P,TRAIT)],"dnanexus/CDS.single.processed.regenie",col.names = T, row.names=F, sep = "\t", quote=F)

# Note: NC.single.wanted.regenie.gz is a file based on NC.single.regenie.txt.gz but filtered for having only the variants associated with significant genes in NC
# If not using NC, this input is not valid

