#05-Mar-2023 Diogo Ribeiro @ UNIL
# Summary of SAIGE-GENE results across traits

library(data.table)
library(ggplot2)

# data1 = fread("~/burdenNC/data/results/200k/ABC_SKATO_process_FDR0.05.tsv")
# data2 = fread("~/burdenNC/data/results/200k/CRD_SKATO_process_FDR0.05.tsv")
# data3 = fread("~/burdenNC/data/results/200k/JAVIERRE_SKATO_process_FDR0.05.tsv")

# data1 = fread("~/burdenNC/data/results/200k/ABC_SKATO_process_PVAL1e9.tsv")
# data2 = fread("~/burdenNC/data/results/200k/CRD_SKATO_process_PVAL1e9.tsv")
# data3 = fread("~/burdenNC/data/results/200k/JAVIERRE_SKATO_process_PVAL1e9.tsv")
# 
# data1$annot = "ABC"
# data2$annot = "CRD"
# data3$annot = "HIC"

# data1 = fread("~/burdenNC/data/results/200k/NC_SKATO_process_FDR0.05.tsv")
# data2 = fread("~/burdenNC/data/results/200k/CDS_SKATO_process_FDR0.05.tsv")

data1 = fread("~/burdenNC/data/results/200k/NC_SKATO_process_PVAL1e9.tsv")
data2 = fread("~/burdenNC/data/results/200k/CDS_SKATO_process_PVAL1e9.tsv")

data1$annot = "NC"
data2$annot = "CDS"


data = rbind(data1,data2)
# data = rbind(data1,data2,data3)
data = unique(data)

data[TRAITNAME == "Mean corpuscular haemoglobin concentration"]$TRAITNAME = "Mean heamoglobin concentration"
data[TRAITNAME == "High light scatter reticulocyte percentage"]$TRAITNAME = "High light scatter reticulocyte %"


dt = unique(data[,.(GENE,TRAITNAME,annot)])
table(dt$annot)

meltedData = melt(data.table(table(dt$TRAITNAME, dt$annot)), id.vars = c("V1", "V2"), measure.vars = "N")
rank = data.table(table(dt$TRAITNAME))
order = rank[order(N)]$V1

colnames(meltedData) = c("Trait","Annotation","N","genes")
meltedData$Trait <- factor(meltedData$Trait, levels = order)

# png("~/git/burdenNC/src/figures/traits_plot_fdr5_3annot.png",1300,2000)
# png("~/git/burdenNC/src/figures/traits_plot_p1e9_3annot.png",1300,2000)
# png("~/git/burdenNC/src/figures/traits_plot_fdr5_cds.png",1300,2000)
png("~/git/burdenNC/src/figures/traits_plot_p1e9_cds.png",1300,2000)

ggplot(data = meltedData, aes(y = Trait, x = genes, fill = Annotation) ) +
  geom_bar(stat = "identity", alpha = 0.8, size = 1, color = "black") +
  xlab("Gene associations") +
  # geom_text(aes(label = genes)) +
  # scale_fill_manual(values = c("#756bb1","#31a354","#de2d26"), limits = c("CRD","ABC","HIC")) +
  scale_fill_manual(values = c("#2171b5","#525252"), limits = c("NC","CDS")) +
  theme_minimal() +
  theme(legend.position = c(0.70,0.5), legend.background = element_rect(fill = "white"), axis.text.y = element_text( color = "black"), #face="bold",
        text = element_text(size=36), panel.grid.minor.y = element_line(size = 0.03, color = "black"),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),  panel.grid.major.y = element_line(size = 0.03, color = "black"),
        panel.background = element_rect(colour = "black", fill = "white", size = 3)
  )

dev.off()



# # known = fread("~/burdenNC/data/known_replication/genebass_gwas_catalog_75traits.tsv.gz")
# known = fread("~/burdenNC/data/known_replication/genebass_75traits_significant_lof_missense_simple.tsv.gz")
# known = fread("~/burdenNC/data/known_replication/gwas_catalog_all_gene_75trait.tsv.gz")
# known$TRAIT = data.table(unlist(lapply(known$trait_id, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
# knownDT = data.table(table(known$TRAIT))
# p = data.table(table(d[annot == "NC"]$TRAIT))
# m = merge(p,knownDT, by = "V1")
# cor.test(m$N.x, m$N.y)
# plot(m$N.x, m$N.y)

