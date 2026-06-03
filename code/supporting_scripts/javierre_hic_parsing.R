
library(data.table)
library(tidyr)
library(dplyr)

inFile = "javierre2016/PCHiC_peak_matrix_cutoff5.tsv"
data = fread( inFile, stringsAsFactors = FALSE, header = T, sep="\t")

unpackData = data %>% 
  mutate(baitName = strsplit(as.character(baitName), ";")) %>% 
  unnest(baitName)

unpackData = data.table(unpackData)
# write.table(unique(unpackData$baitName),"gene_list_split_unpack.txt",quote=F,col.names = F,row.names = F)

nrow(data)
nrow(unpackData)

geneMapping = fread( "gprofiler_gene_id_convert_simple.txt", stringsAsFactors = FALSE, header = F, sep="\t")
geneMapping = geneMapping[V2 != "None"][V1 != "Y_RNA"]

mergedData = merge(unpackData,geneMapping, by.x = "baitName", by.y = "V1")

# Only retain interactions in same chromosome
mergedData = mergedData[baitChr == oeChr]

# Only considering contacts within 1MB of gene
mergedData = mergedData[dist > -1000000][dist < 1000000]

# Only considering autosomes
mergedData = mergedData[oeChr != "X"][oeChr != "Y"]

write.table(mergedData[,.(oeChr,oeStart,oeEnd,V2)][order(oeChr,oeStart)],"javierre2016_chicago_above5_all_celltypes_no_interchrom_autosome_1mb_limit.tsv",sep="\t",quote=F,col.names=F,row.names=F)

