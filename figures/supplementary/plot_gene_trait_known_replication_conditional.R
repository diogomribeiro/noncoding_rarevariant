
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png/known_replication_conditional.png",1400,700)

data = fread("~/git/noncoding_rvat/revision/known_replication/conditional/full.conditional.summary")
data$pval_cutoff = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[3])))$V1
data$annotation = data.table(unlist(lapply(data$V3, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
data[annotation == "NC_shuffle"]$annotation = "EXPEC."
data[annotation == "JAVIERRE"]$annotation = "HIC"
data[annotation == "PROMOTER"]$annotation = "PROM."
colnames(data) = c("replication","percentage","V3","pval_cutoff","annotation")

# Direct replication
g1 = ggplot(data[replication == "direct_replication"][pval_cutoff == "P1e9"], aes(x=annotation, y=percentage, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  scale_x_discrete(limits = c("NC","CDS","PROM.","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "a") +
  guides(fill = "none") +
  ylim(c(0,102)) +
  scale_fill_manual(values = c("#2171b5","#525252","#02818a","#fed976"), limits = c("NC","CDS","PROM.","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=36), legend.position = c(0.40,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))

g2 = ggplot(data[replication == "window_replication"][pval_cutoff == "P1e9"], aes(x=annotation, y=percentage, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(round(percentage,1),"%")), size = 8, fontface= "bold", vjust = -0.5) +
  scale_x_discrete(limits = c("NC","CDS","PROM.","EXPEC.")) +
  labs(x="Annotation", y = "% known gene-trait", tag = "b") +
  guides(fill = "none") +
  ylim(c(0,102)) +
  scale_fill_manual(values = c("#2171b5","#525252","#02818a","#fed976"), limits = c("NC","CDS","PROM.","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=36), legend.position = c(0.40,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3))

g1 + g2

dev.off()
