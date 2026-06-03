
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png/epimap_known_replication.png",2000,1400)


nAssoc = fread("~/git/noncoding_rvat/revision/epimap/known_replication/number_assoc.txt")

data = fread("~/git/noncoding_rvat/revision/epimap/known_replication/full.summary")
data$pval_cutoff = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[2])))$V1
data$annotation = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
data[annotation == "EPI_BLOOD_T_CELL"]$annotation = "T cell"
data[annotation == "EPI_ENDO"]$annotation = "Endothelial"
data[annotation == "EPI_HSC_B_CELL"]$annotation = "B cell"
data[annotation == "EPI_KIDNEY"]$annotation = "Kidney"
data[annotation == "EPI_LIVER"]$annotation = "Liver"
data[annotation == "EPI_PANCREAS"]$annotation = "Pancreas"
colnames(data) = c("replication","percentage","V3","pval_cutoff","annotation")

data = data[pval_cutoff == "P1e9"]
data = merge(data, nAssoc, by = "annotation")

# Direct replication
g1 =  ggplot(data[replication == "direct_replication"], aes(x=annotation, y=percentage, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", label = paste("N =",data$no_cond), size = 8, fontface= "bold", x = data$annotation, y = -3, vjust = 0.5) +
  scale_x_discrete(limits = c("T cell","B cell", "Kidney","Liver","Pancreas","Endothelial")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "a", title = "Before conditional | Direct") +
  guides(fill = "none") +
  ylim(c(-3,105)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#d9d9d9"), limits = c("T cell","B cell", "Kidney","Liver","Pancreas","Endothelial")) +
  theme_linedraw() +
  theme(text = element_text(size=36), legend.position = c(0.40,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))

# Window replication
g2 =  ggplot(data[replication == "window_replication"], aes(x=annotation, y=percentage, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", label = paste("N =",data$no_cond), size = 8, fontface= "bold", x = data$annotation, y = -3, vjust = 0.5) +
  scale_x_discrete(limits = c("T cell","B cell", "Kidney","Liver","Pancreas","Endothelial")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "b", title = "Before conditional | Window") +
  guides(fill = "none") +
  ylim(c(-3,105)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#d9d9d9"), limits = c("T cell","B cell", "Kidney","Liver","Pancreas","Endothelial")) +
  theme_linedraw() +
  theme(text = element_text(size=36), legend.position = c(0.40,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))



data = fread("~/git/noncoding_rvat/revision/epimap/known_replication/full.conditional.summary")
data$pval_cutoff = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[2])))$V1
data$annotation = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
data[annotation == "EPI_BLOOD_T_CELL"]$annotation = "T cell"
data[annotation == "EPI_ENDO"]$annotation = "Endothelial"
data[annotation == "EPI_HSC_B_CELL"]$annotation = "B cell"
data[annotation == "EPI_KIDNEY"]$annotation = "Kidney"
data[annotation == "EPI_LIVER"]$annotation = "Liver"
data[annotation == "EPI_PANCREAS"]$annotation = "Pancreas"
colnames(data) = c("replication","percentage","V3","pval_cutoff","annotation")

data = data[pval_cutoff == "P1e9"]
data = merge(data, nAssoc, by = "annotation")

data = rbind(data, data.table(annotation = "Endothelial",replication = "direct_replication",percentage = 0, V3 = "bla",pval_cutoff = "P1e9", no_cond = 9, cond = 0))
data = rbind(data, data.table(annotation = "Endothelial",replication = "window_replication",percentage = 0, V3 = "bla",pval_cutoff = "P1e9", no_cond = 9, cond = 0))


# Direct replication
g3 =  ggplot(data[replication == "direct_replication"], aes(x=annotation, y=percentage, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", label = paste("N =",data$cond), size = 8, fontface= "bold", x = data$annotation, y = -3, vjust = 0.5) +
  scale_x_discrete(limits = c("T cell","B cell", "Kidney","Liver","Pancreas","Endothelial")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "c", title = "After conditional | Direct") +
  guides(fill = "none") +
  ylim(c(-3,105)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#d9d9d9"), limits = c("T cell","B cell", "Kidney","Liver","Pancreas","Endothelial")) +
  theme_linedraw() +
  theme(text = element_text(size=36), legend.position = c(0.40,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))

# Window replication
g4 =  ggplot(data[replication == "window_replication"], aes(x=annotation, y=percentage, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  annotate(geom = "text", label = paste("N =",data$cond), size = 8, fontface= "bold", x = data$annotation, y = -3, vjust = 0.5) +
  scale_x_discrete(limits = c("T cell","B cell", "Kidney","Liver","Pancreas","Endothelial")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "d", title = "After conditional | Window") +
  guides(fill = "none") +
  ylim(c(-3,105)) +
  scale_fill_manual(values = c("#238b45","#6a51a3","#cb181d","#2171b5","#525252","#d9d9d9"), limits = c("T cell","B cell", "Kidney","Liver","Pancreas","Endothelial")) +
  theme_linedraw() +
  theme(text = element_text(size=36), legend.position = c(0.40,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))

g1 + g2 + g3 + g4

dev.off()
