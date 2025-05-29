
library(data.table)
library(ggplot2)
library(patchwork)

png("~/git/noncoding_rvat/figures/revision/png/traits_plot_nc_cds.png",2400,2000)

data1 = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/NC_SKATO_process_FDR0.05.tsv")
data2 = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/CDS_SKATO_process_FDR0.05.tsv")

data3 = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/NC_SKATO_process_PVAL1e9.tsv")
data4 = fread("~/git/noncoding_rvat/revision/no_cond_42_traits/CDS_SKATO_process_PVAL1e9.tsv")

data1$annot = "NC"
data2$annot = "CDS"
data3$annot = "NC"
data4$annot = "CDS"

###############
# FDR 5%
###############

data = unique(rbind(data1,data2))

data[TRAITNAME == "Mean corpuscular haemoglobin concentration"]$TRAITNAME = "Mean heamoglobin concentration"
data[TRAITNAME == "High light scatter reticulocyte percentage"]$TRAITNAME = "High light scatter reticulocyte %"

dt = unique(data[,.(GENE,TRAITNAME,annot)])
meltedData = melt(data.table(table(dt$TRAITNAME, dt$annot)), id.vars = c("V1", "V2"), measure.vars = "N")
rank = data.table(table(dt$TRAITNAME))
order = rank[order(N)]$V1

colnames(meltedData) = c("Trait","Annotation","N","genes")
meltedData$Trait <- factor(meltedData$Trait, levels = order)

g1 = ggplot(data = meltedData, aes(y = Trait, x = genes, fill = Annotation) ) +
  geom_bar(stat = "identity", alpha = 0.8, size = 1, color = "black") +
  xlab("Gene associations") +
  scale_fill_manual(values = c("#2171b5","#525252"), limits = c("NC","CDS")) +
  theme_minimal() +
  labs(tag = "a", title = "FDR<5%") +
  theme(legend.position = c(0.70,0.5), legend.background = element_rect(fill = "white"), axis.text.y = element_text( color = "black"), #face="bold",
        text = element_text(size=36), panel.grid.minor.y = element_line(size = 0.03, color = "black"), plot.title = element_text(hjust = 0.5),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),  panel.grid.major.y = element_line(size = 0.03, color = "black"),
        panel.background = element_rect(colour = "black", fill = "white", size = 3)
  )


###############
# P1e9
###############

data = unique(rbind(data3,data4))

data[TRAITNAME == "Mean corpuscular haemoglobin concentration"]$TRAITNAME = "Mean heamoglobin concentration"
data[TRAITNAME == "High light scatter reticulocyte percentage"]$TRAITNAME = "High light scatter reticulocyte %"

dt = unique(data[,.(GENE,TRAITNAME,annot)])
meltedData = melt(data.table(table(dt$TRAITNAME, dt$annot)), id.vars = c("V1", "V2"), measure.vars = "N")
rank = data.table(table(dt$TRAITNAME))
order = rank[order(N)]$V1

colnames(meltedData) = c("Trait","Annotation","N","genes")
meltedData$Trait <- factor(meltedData$Trait, levels = order)

g2 = ggplot(data = meltedData, aes(y = Trait, x = genes, fill = Annotation) ) +
  geom_bar(stat = "identity", alpha = 0.8, size = 1, color = "black") +
  xlab("Gene associations") +
  scale_fill_manual(values = c("#2171b5","#525252"), limits = c("NC","CDS")) +
  labs(tag = "b", title = "P<1e-9") +
  theme_minimal() +
  theme(legend.position = c(0.70,0.5), legend.background = element_rect(fill = "white"), axis.text.y = element_text( color = "black"), #face="bold",
        text = element_text(size=36), panel.grid.minor.y = element_line(size = 0.03, color = "black"), plot.title = element_text(hjust = 0.5),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),  panel.grid.major.y = element_line(size = 0.03, color = "black"),
        panel.background = element_rect(colour = "black", fill = "white", size = 3)
  )

g1 + g2

dev.off()

