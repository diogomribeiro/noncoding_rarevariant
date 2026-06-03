
library(data.table)
library(ggplot2)
library(patchwork)

# png("~/git/noncoding_rvat/figures/final/figure5.png",1500,750)
pdf("~/git/noncoding_rvat/figures/final/figure5.pdf", width = 20.83, height = 10.42)

##############################
# Panel A
##############################

# No conditional files
cdsData137k = fread("zcat ~/burdenNC/data/results/137k_30k/CDS_137k_SKATO_process_all.tsv.gz")
cdsData30k = fread("zcat  ~/burdenNC/data/results/137k_30k/CDS_30k_SKATO_process_all.tsv.gz")
ncData137k = fread("zcat ~/burdenNC/data/results/137k_30k/NC_137k_SKATO_process_all.tsv.gz")
ncData30k = fread("zcat  ~/burdenNC/data/results/137k_30k/NC_30k_SKATO_process_all.tsv.gz")
shuffle137k = fread("zcat /home/dribeiro/git/noncoding_rvat/revision/downsample/no_cond_42_traits/137k/NC_shuffle_SKATO_process_all.tsv.gz")

# # Conditional files, note that this was only run for a subset of genes
cdsData137kCond = fread("~/git/noncoding_rvat/revision/downsample/CDS_137k_SKATO_conditional_all.tsv")
cdsData30kCond = fread("~/git/noncoding_rvat/revision/downsample/CDS_30k_SKATO_conditional_all.tsv")
ncData137kCond = fread("~/git/noncoding_rvat/revision/downsample/NC_137k_SKATO_conditional_all.tsv")
ncData30kCond = fread("~/git/noncoding_rvat/revision/downsample/NC_30k_SKATO_conditional_all.tsv")
shuffle137kCond = fread("~/git/noncoding_rvat/revision/downsample/NC_137k_SKATO_conditional_all_shuffle.tsv")

wantedTraits = fread("~/git/noncoding_rvat/supporting_data/42_traits.txt", header = F)

cdsData137k = cdsData137k[TRAIT %in% wantedTraits$V1]
cdsData30k = cdsData30k[TRAIT %in% wantedTraits$V1]
ncData137k = ncData137k[TRAIT %in% wantedTraits$V1]
ncData30k = ncData30k[TRAIT %in% wantedTraits$V1]
shuffle137k = shuffle137k[TRAIT %in% wantedTraits$V1]

###########
## CDS
###########
mergedCDS = merge(cdsData30k,cdsData137k,by=c("GENE","TRAIT"), all.x = T, all.y = T)
paste(nrow(unique(mergedCDS[is.na(PVAL.x)])), "CDS gene-traits not evaluated in 30k")
paste(nrow(unique(mergedCDS[is.na(PVAL.y)])), "CDS gene-traits not evaluated in 137k")
paste(nrow(unique(mergedCDS[!is.na(PVAL.x)][!is.na(PVAL.y)])),"NC gene-traits compared")
mergedCDS = mergedCDS[!is.na(PVAL.x)][!is.na(PVAL.y)]

###########
# NC
###########
mergedNC = merge(ncData30k,ncData137k,by=c("GENE","TRAIT"), all.x = T, all.y = T)
paste(nrow(unique(mergedNC[is.na(PVAL.x)])), "NC gene-traits not evaluated in 30k")
paste(nrow(unique(mergedNC[is.na(PVAL.y)])), "NC gene-traits not evaluated in 137k")
paste(nrow(unique(mergedNC[!is.na(PVAL.x)][!is.na(PVAL.y)])),"NC gene-traits compared")
mergedNC = mergedNC[!is.na(PVAL.x)][!is.na(PVAL.y)]

###########
# Shuffled NC
###########
mergedSH = merge(ncData30k,shuffle137k,by=c("GENE","TRAIT"), all.x = T, all.y = T)
paste(nrow(unique(mergedSH[is.na(PVAL.x)])), "NC gene-traits not evaluated in 30k")
paste(nrow(unique(mergedSH[is.na(PVAL.y)])), "NC gene-traits not evaluated in 137k")
paste(nrow(unique(mergedSH[!is.na(PVAL.x)][!is.na(PVAL.y)])),"NC gene-traits compared")
mergedSH = mergedSH[!is.na(PVAL.x)][!is.na(PVAL.y)]

###########
## CDS COND
###########
mergedCDSCond = merge(cdsData30kCond,cdsData137kCond,by=c("GENE","TRAIT"), all.x = T, all.y = T)
paste(nrow(unique(mergedCDSCond[is.na(PVAL.x)])), "CDS gene-traits not evaluated in 30k")
paste(nrow(unique(mergedCDSCond[is.na(PVAL.y)])), "CDS gene-traits not evaluated in 137k")
paste(nrow(unique(mergedCDSCond[!is.na(PVAL.x)][!is.na(PVAL.y)])),"NC gene-traits compared")
mergedCDSCond = mergedCDSCond[!is.na(PVAL.x)][!is.na(PVAL.y)]

###########
# NC COND
###########
mergedNCCond = merge(ncData30kCond,ncData137kCond,by=c("GENE","TRAIT"), all.x = T, all.y = T)
paste(nrow(unique(mergedNCCond[is.na(PVAL.x)])), "NC gene-traits not evaluated in 30k")
paste(nrow(unique(mergedNCCond[is.na(PVAL.y)])), "NC gene-traits not evaluated in 137k")
paste(nrow(unique(mergedNCCond[!is.na(PVAL.x)][!is.na(PVAL.y)])),"NC gene-traits compared")
mergedNCCond = mergedNCCond[!is.na(PVAL.x)][!is.na(PVAL.y)]

###########
# Shuffled NC
###########
mergedSHCond = merge(ncData30kCond,shuffle137kCond,by=c("GENE","TRAIT"), all.x = T, all.y = T)
paste(nrow(unique(mergedSHCond[is.na(PVAL.x)])), "NC gene-traits not evaluated in 30k")
paste(nrow(unique(mergedSHCond[is.na(PVAL.y)])), "NC gene-traits not evaluated in 137k")
paste(nrow(unique(mergedSHCond[!is.na(PVAL.x)][!is.na(PVAL.y)])),"NC gene-traits compared")
mergedSHCond = mergedSHCond[!is.na(PVAL.x)][!is.na(PVAL.y)]


############
# Replication rate
############
# ncFDR0.05 = nrow(mergedNC[FDR.y < 0.05][PVAL.x < 0.05]) * 100 / nrow(mergedNC[FDR.y < 0.05])
# cdsFDR0.05 = nrow(mergedCDS[FDR.y < 0.05][PVAL.x < 0.05]) * 100 / nrow(mergedCDS[FDR.y < 0.05])
# shFDR0.05 = nrow(mergedSH[FDR.y < 0.05][PVAL.x < 0.05]) * 100 / nrow(mergedSH[FDR.y < 0.05])
ncP1e9 = nrow(mergedNC[PVAL.y < 1e-9][PVAL.x < 0.05]) * 100 / nrow(mergedNC[PVAL.y < 1e-9])
cdsP1e9 = nrow(mergedCDS[PVAL.y < 1e-9][PVAL.x < 0.05]) * 100 / nrow(mergedCDS[PVAL.y < 1e-9])
shP1e9 = nrow(mergedSH[PVAL.y < 1e-9][PVAL.x < 0.05]) * 100 / nrow(mergedSH[PVAL.y < 1e-9])

ncP1e9Cond = nrow(mergedNCCond[PVAL.y < 1e-9][PVAL.x < 0.05]) * 100 / nrow(mergedNCCond[PVAL.y < 1e-9])
cdsP1e9Cond = nrow(mergedCDSCond[PVAL.y < 1e-9][PVAL.x < 0.05]) * 100 / nrow(mergedCDSCond[PVAL.y < 1e-9])
shP1e9Cond = nrow(mergedSHCond[PVAL.y < 1e-9][PVAL.x < 0.05]) * 100 / nrow(mergedSHCond[PVAL.y < 1e-9])

data = data.table(annotation = c("NC", "CDS","EXPECTED","NC","CDS","EXPECTED"), 
                  percentage = c(ncP1e9,cdsP1e9,shP1e9,ncP1e9Cond,cdsP1e9Cond,shP1e9Cond), conditional = c("no","no","no","yes","yes","yes"))

g1 = ggplot(data, aes(x=annotation, y=percentage, alpha = conditional, fill = annotation )) +
  geom_bar( stat = "identity", position = "dodge", size = 1, color = "black") +
  geom_text(aes(label = paste0(round(percentage,1),"%"), group = conditional), position = position_dodge(width = .9), vjust = -0.5, alpha = 1, size = 9, fontface = "bold" ) +
  # geom_text(aes(label = paste0(round(percentage,1),"%"), group = cutoff), position = position_dodge(width = .9), vjust = -0.5, alpha = 1, size = 9, fontface = "bold" ) +
  annotate(geom = "text", x = "NC", y = -3, label = paste0(" N=",nrow(mergedNC[PVAL.y < 1e-9]), "   N=",nrow(mergedNCCond[PVAL.y < 1e-9])), size = 8) +
  annotate(geom = "text", x = "CDS", y = -3, label =  paste0(" N=",nrow(mergedCDS[PVAL.y < 1e-9]), "   N=",nrow(mergedCDSCond[PVAL.y < 1e-9])), size = 8) +
  annotate(geom = "text", x = "EXPECTED", y = -3, label = paste0("N=165   N=16"), size = 8) +
  scale_x_discrete(limits = c("NC","CDS","EXPECTED")) +
  labs(x="Annotation", y = "% replicated gene-trait", tag = "a") +
  guides(fill = "none") +
  ylim(c(-3, 100)) +
  scale_fill_manual(values = c("#2171b5","#525252","grey"), limits = c("NC","CDS","EXPEC.")) +
  theme_linedraw() +
  theme(text = element_text(size=32), legend.position = c(0.86,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)

d1=data

##############################
# Panel B
##############################

## P1e9
data = fread("~/git/noncoding_rvat/revision/downsample/summary_pval1e9.tsv")

data = melt(data)
data[Dataset == "30k"]$Dataset = 26212
data[Dataset == "83k"]$Dataset = 83370
data[Dataset == "137k"]$Dataset = 136740
data[Dataset == "166k"]$Dataset = 166740
data$Dataset = as.numeric(data$Dataset)
colnames(data) = c("sample_size","annotation","results")

cor.test(data[annotation == "NC"]$results, data[annotation == "NC"]$sample_size)
cor.test(data[annotation == "CDS"]$results, data[annotation == "CDS"]$sample_size)

# NC
modelNC = lm(results ~ sample_size, data = data[annotation == "NC"])
data = rbind(data, data.table(sample_size = 400000, annotation = "NC", results = predict(modelNC, data.table(sample_size = 400000))))
# CDS
modelCDS = lm(results ~ sample_size, data = data[annotation == "CDS"])
data = rbind(data, data.table(sample_size = 400000, annotation = "CDS", results = predict(modelCDS, data.table(sample_size = 400000))))
# NC COND
modelNCCond = lm(results ~ sample_size, data = data[annotation == "NC_COND"])
data = rbind(data, data.table(sample_size = 400000, annotation = "NC_COND", results = predict(modelNCCond, data.table(sample_size = 400000))))
# CDS COND
modelCDSCond = lm(results ~ sample_size, data = data[annotation == "CDS_COND"])
data = rbind(data, data.table(sample_size = 400000, annotation = "CDS_COND", results = predict(modelCDSCond, data.table(sample_size = 400000))))

d2=data

options(scipen = 1)
g2 = ggplot() + 
  geom_abline(intercept = coef(modelNC)[1], slope = coef(modelNC)[2], color = "#2171b5", size = 1.5, linetype = "dashed") +
  geom_abline(intercept = coef(modelCDS)[1], slope = coef(modelCDS)[2], color = "#525252", size = 1.5, linetype = "dashed") +
  geom_abline(intercept = coef(modelNCCond)[1], slope = coef(modelNCCond)[2], color = "#2171b5", size = 1.5) +
  geom_abline(intercept = coef(modelCDSCond)[1], slope = coef(modelCDSCond)[2], color = "#525252", size = 1.5) +
  geom_point(data = data[annotation == "CDS" | annotation == "NC"], aes(x = sample_size, y = results, color = annotation), size = 5, shape = 1) +
  geom_point(data = data[annotation == "CDS_COND" | annotation == "NC_COND"], aes(x = sample_size, y = results, color = annotation), size = 5, shape = 17) +
  geom_point(data = data[annotation == "CDS" | annotation == "NC"][sample_size == 400000], aes(x = sample_size, y = results, color = annotation), size = 5, shape = 19, color = "black") +
  geom_point(data = data[annotation == "CDS_COND" | annotation == "NC_COND"][sample_size == 400000], aes(x = sample_size, y = results, color = annotation), size = 5, shape = 17, color = "black") +
  scale_color_manual(name = "", values = c("NC" = "#2171b5","CDS" = "#525252","NC_COND" = "#2171b5","CDS_COND" = "#525252")) +
  labs(x="Sample size", y = "Gene-trait associations", tag = "b") +
  xlim(c(0,420000)) +
  theme_linedraw() +
  theme(text = element_text(size=32), legend.position = c(0.36,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  


g1 + g2

dev.off()


## Source data
library(openxlsx)
write.xlsx(list(fig5a = d1, fig5b = d2), "~git/noncoding_rvat/figures/final/source_data/Ribeiro_SourceData_Fig5.xlsx")

