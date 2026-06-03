
library(data.table)
library(ggplot2)

png("~/git/noncoding_rvat/figures/revision/png/coding_conditional.png",800,800)

# ADD TRAIT NAME AND GENE NAME # THIS CAN APPLY TRAIT FILTER
traitMapping = fread("~/git/noncoding_rvat/supporting_data/42_traits_mapping.txt", header = F)
traitMapping$TRAIT = data.table(unlist(lapply(traitMapping$V1, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
geneMapping = fread("~/git/noncoding_rvat/supporting_data/gencode_id_mapping.tsv.gz")

cdsCondData = fread("~/git/noncoding_rvat/revision/coding_conditional/NC_CDScond.all_chr.regenie.gz")
cdsCondData$GENE = data.table(unlist(lapply(cdsCondData$ID, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
cdsCondData$TRAIT = data.table(unlist(lapply(cdsCondData$TRAIT, function(x) unlist(strsplit(x,"[.]"))[4])))$V1
cdsCondData$TRAIT = data.table(unlist(lapply(cdsCondData$TRAIT, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
cdsCondData = cdsCondData[TEST == "ADD-SKATO"][,.(CHROM,GENPOS,GENE,LOG10P,TRAIT)]

cdsCondData = merge(cdsCondData, traitMapping[,.(TRAIT,V2)], by = "TRAIT")
cdsCondData = merge(cdsCondData, geneMapping[,.(V4,V6)], by.x = "GENE", by.y = "V4")

checkData = fread("~/git/noncoding_rvat/revision/coding_conditional/test/check2/NC.all_chr.regenie.gz")
checkData$GENE = data.table(unlist(lapply(checkData$ID, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
checkData$TRAIT = data.table(unlist(lapply(checkData$TRAIT, function(x) unlist(strsplit(x,"[.]"))[4])))$V1
checkData$TRAIT = data.table(unlist(lapply(checkData$TRAIT, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
checkData = checkData[TEST == "ADD-SKATO"][,.(CHROM,GENPOS,GENE,LOG10P,TRAIT)]

checkData = merge(checkData, traitMapping[,.(TRAIT,V2)], by = "TRAIT")
checkData = merge(checkData, geneMapping[,.(V4,V6)], by.x = "GENE", by.y = "V4")

## How does the before and after CDS conditional compare

mergeData2 = merge(checkData,cdsCondData, by = c("GENE","TRAIT"))
plot(mergeData2$LOG10P.x, mergeData2$LOG10P.y, xlab = "GWAS conditional only", ylab = "GWAS + CDS conditional")
cor.test(mergeData2$LOG10P.x, mergeData2$LOG10P.y)

wantedAssoc = fread("~/git/noncoding_rvat/revision/cond_rerun/NC_PVAL1e9_gene_trait_remaining.txt", header = F)
mergedData2 = mergeData2[GENE %in% wantedAssoc$V1 & TRAIT %in% wantedAssoc$V2]

test = cor.test(mergeData2$LOG10P.x,mergeData2$LOG10P.y, method = "spearman")
text = paste("Spearman R:",round(test$estimate,3), "P-value",format.pval(test$p.value,2))
ggplot(data = mergeData2, aes(LOG10P.x, LOG10P.y)) + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  geom_point(size = 3) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 8, fontface = "bold"  ) +
  ylab("GWAS + CDS conditional (-log10)") +
  xlab("GWAS conditional (-log10)") +
  xlim(0,34) +
  ylim(0,34) +
  theme_minimal() + 
  theme(text = element_text(size=32), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

dev.off()