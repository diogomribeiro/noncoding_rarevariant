
library(data.table)

dataset = "NC" # ABC CRD JAVIERRE CDS CONTROL

regenieData = fread(paste0("zcat ~/burdenNC/data/results/200k/",dataset,".all_chr.regenie.gz"), sep = " ")
regenieData$GENE = data.table(unlist(lapply(regenieData$ID, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
regenieData$TRAIT = data.table(unlist(lapply(regenieData$TRAIT, function(x) unlist(strsplit(x,"[.]"))[2])))$V1
regenieData$TRAIT = data.table(unlist(lapply(regenieData$TRAIT, function(x) unlist(strsplit(x,"[_]"))[2])))$V1
regenieDataFilt = regenieData[TEST == "ADD-SKATO"][,.(CHROM,GENPOS,GENE,LOG10P,TRAIT)]

# ADD TRAIT NAME AND GENE NAME # THIS CAN APPLY TRAIT FILTER
# traitMapping = fread("~/git/burdenNC/src/saige/saige_step1/phenotype_list_mapping2", header = F)
traitMapping = fread("~/burdenNC/data/supporting_data/wanted_traits.txt", header = F)
traitMapping$TRAIT = data.table(unlist(lapply(traitMapping$V1, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
geneMapping = fread("~/burdenNC/data/supporting_data/gencode_id_mapping.tsv")

regenieDataFilt = merge(regenieDataFilt, traitMapping[,.(TRAIT,V2)], by = "TRAIT")
regenieDataFilt = merge(regenieDataFilt, geneMapping[,.(V4,V6)], by.x = "GENE", by.y = "V4")

colnames(regenieDataFilt)[6] = "TRAITNAME"
colnames(regenieDataFilt)[7] = "GENENAME"

# FDR
regenieDataFilt$PVAL = 10^(-regenieDataFilt$LOG10P)

## For creating shuffled dataset
# scramble_set <- function(input_dt, colname) { set(input_dt, j = colname, value = sample(input_dt[[colname]])) }
# regenieDataFilt = scramble_set(regenieDataFilt, "PVAL")

regenieDataFilt$FDR = p.adjust(regenieDataFilt$PVAL, method = "BH")

# # Write all associations
write.table(regenieDataFilt,paste0("~/burdenNC/data/results/200k/",dataset,"_SKATO_process_all.tsv"), quote = F, sep = "\t", row.names = F)
# Write significant FDR
write.table(regenieDataFilt[FDR < 0.05],paste0("~/burdenNC/data/results/200k/",dataset,"_SKATO_process_FDR0.05.tsv"), quote = F, sep = "\t", row.names = F)
# Write significant 1e-9
write.table(regenieDataFilt[PVAL < 1e-9],paste0("~/burdenNC/data/results/200k/",dataset,"_SKATO_process_PVAL1e9.tsv"), quote = F, sep = "\t", row.names = F)

