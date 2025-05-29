
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png/known_replication_fdr0.05.png",2000,800)

nAssoc = fread("~/git/noncoding_rvat/revision/known_replication/sample_size_gene_trait_assoc.txt")

data = fread("~/git/noncoding_rvat/revision/known_replication/full.summary")
data$pval_cutoff = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[2])))$V1
data$annotation = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
data[annotation == "NC_shuffle"]$annotation = "EXPEC."
data[annotation == "PROMOTER"]$annotation = "PROM."
data[annotation == "JAVIERRE"]$annotation = "HIC"
colnames(data) = c("replication","percentage","V3","pval_cutoff","annotation")

data = merge(data, nAssoc, by = "annotation")


# Direct replication
g1 = ggplot(data[replication == "direct_replication"][pval_cutoff == "FDR0"], aes(x=annotation, y=percentage, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", label = paste("N =",data$FDR5), size = 8, fontface= "bold", x = data$annotation, y = -3, vjust = 0.5) +
  scale_x_discrete(limits = c("ABC","CRD","HIC","NC","CDS","PROM.","CONTROL","RANDOM","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "b") +
  guides(fill = "none") +
  ylim(c(-3,105)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#02818a","#d9d9d9","black","#fed976"), limits = c("ABC","CRD","HIC","NC","CDS","PROM.","CONTROL","RANDOM","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=36), legend.position = c(0.40,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))

# Window replication
g2 = ggplot(data[replication == "window_replication"][pval_cutoff == "FDR0"], aes(x=annotation, y=percentage, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", label = paste("N =",data$FDR5), size = 8, fontface= "bold", x = data$annotation, y = -3, vjust = 0.5) +
  scale_x_discrete(limits = c("ABC","CRD","HIC","NC","CDS","PROM.","CONTROL","RANDOM","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "c") +
  guides(fill = "none") +
  ylim(c(-3,105)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#02818a","#d9d9d9","black","#fed976"), limits = c("ABC","CRD","HIC","NC","CDS","PROM.","CONTROL","RANDOM","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=36), legend.position = c(0.40,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))

g1 + g2

dev.off()
