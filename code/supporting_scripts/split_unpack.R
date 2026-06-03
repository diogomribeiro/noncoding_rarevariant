library(data.table)
library(tidyr)
library(dplyr)
args=commandArgs(trailingOnly=TRUE)
inFile = args[1]
data = fread( inFile, stringsAsFactors = FALSE, header = F, sep="\t")
newData = data %>% mutate(V4 = strsplit(as.character(V4), ",")) %>% unnest(V4)
write.table(data.table(newData),args[2],col.names=F,quote=F,sep="\t",row.names=F)

