
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png/variant_per_annot.png",1200,800)

data1 = fread("zcat ~/burdenNC/data/annotations/ABC.CADD15.annotations.gz", header = F)
data2 = fread("zcat ~/burdenNC/data/annotations/CRD.CADD15.annotations.gz", header = F)
data3 = fread("zcat ~/burdenNC/data/annotations/JAVIERRE.CADD15.annotations.gz", header = F)
data4 = fread("zcat ~/burdenNC/data/annotations/NC.CADD15.annotations.gz", header = F)
data5 = fread("zcat ~/burdenNC/data/annotations/CDS.CADD15.annotations.gz", header = F)
data6 = fread("zcat ~/burdenNC/data/annotations/PROMOTER.CADD15.annotations.gz", header = F)

dataLab = "ABC" 
freq = data.table(table(data1$V2))
nGenes = length(unique(data1$V2))
nVars = length(unique(data1$V1))
title = paste(dataLab) #, "\n# Genes:", nGenes, "# Vars:", round(nVars/1000000,2),"M")
text = paste("Mean:", round(mean(freq$N),2), "\nMedian:", median(freq$N))
g1 = ggplot(data = freq[N < 2000], aes(N)) + 
  geom_histogram(color = "black", fill = "grey", bins = 100) +
  ggtitle(title) + 
  ylab("# Genes") +
  xlab("# Variants") +
  labs(tag = "a") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.1, vjust = 1.5, size = 6, fontface = "bold"  ) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 24), text = element_text(size=22), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )

dataLab = "CRD"
freq = data.table(table(data2$V2))
nGenes = length(unique(data2$V2))
nVars = length(unique(data2$V1))
title = paste(dataLab) #, "\n# Genes:", nGenes, "# Vars:", round(nVars/1000000,2),"M")
text = paste("Mean:", round(mean(freq$N),2), "\nMedian:", median(freq$N))
g2 = ggplot(data = freq[N < 2000], aes(N)) + 
  geom_histogram(color = "black", fill = "grey", bins = 100) +
  ggtitle(title) + 
  ylab("# Genes") +
  xlab("# Variants") +
  labs(tag = "b") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.1, vjust = 1.5, size = 6, fontface = "bold"  ) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 24), text = element_text(size=22), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )



dataLab = "HIC"
freq = data.table(table(data3$V2))
nGenes = length(unique(data3$V2))
nVars = length(unique(data3$V1))
title = paste(dataLab) #, "\n# Genes:", nGenes, "# Vars:", round(nVars/1000000,2),"M")
text = paste("Mean:", round(mean(freq$N),2), "\nMedian:", median(freq$N))
g3 = ggplot(data = freq[N < 2000], aes(N)) + 
  geom_histogram(color = "black", fill = "grey", bins = 100) +
  ggtitle(title) + 
  ylab("# Genes") +
  xlab("# Variants") +
  labs(tag = "c") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.1, vjust = 1.5, size = 6, fontface = "bold"  ) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 24), text = element_text(size=22), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )


dataLab = "NC"
freq = data.table(table(data4$V2))
nGenes = length(unique(data4$V2))
nVars = length(unique(data4$V1))
title = paste(dataLab) #, "\n# Genes:", nGenes, "# Vars:", round(nVars/1000000,2),"M")
text = paste("Mean:", round(mean(freq$N),2), "\nMedian:", median(freq$N))
g4 = ggplot(data = freq[N < 2000], aes(N)) + 
  geom_histogram(color = "black", fill = "grey", bins = 100) +
  ggtitle(title) + 
  ylab("# Genes") +
  xlab("# Variants") +
  labs(tag = "d") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.1, vjust = 1.5, size = 6, fontface = "bold"  ) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 24), text = element_text(size=22), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )



dataLab = "CDS"
freq = data.table(table(data5$V2))
nGenes = length(unique(data5$V2))
nVars = length(unique(data5$V1))
title = paste(dataLab) #, "\n# Genes:", nGenes, "# Vars:", round(nVars/1000000,2),"M")
text = paste("Mean:", round(mean(freq$N),2), "\nMedian:", median(freq$N))
g5 = ggplot(data = freq[N < 2000], aes(N)) + 
  geom_histogram(color = "black", fill = "grey", bins = 100) +
  ggtitle(title) + 
  ylab("# Genes") +
  xlab("# Variants") +
  labs(tag = "e") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.1, vjust = 1.5, size = 6, fontface = "bold"  ) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 24), text = element_text(size=22), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )


dataLab = "PROMOTER"
freq = data.table(table(data6$V2))
nGenes = length(unique(data6$V2))
nVars = length(unique(data6$V1))
title = paste(dataLab) #, "\n# Genes:", nGenes, "# Vars:", round(nVars/1000000,2),"M")
text = paste("Mean:", round(mean(freq$N),2), "\nMedian:", median(freq$N))
g6 = ggplot(data = freq[N < 2000], aes(N)) + 
  geom_histogram(color = "black", fill = "grey", bins = 100) +
  ggtitle(title) + 
  ylab("# Genes") +
  xlab("# Variants") +
  labs(tag = "e") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.1, vjust = 1.5, size = 6, fontface = "bold"  ) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 24), text = element_text(size=22), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )


g1 + g2 + g3 + g4 + g5 + g6

dev.off()
