# Plot number of gene-trait associations depending on sample size

library(data.table)
library(ggplot2)

# png("~/git/noncoding_rvat/figures/revision/png/sample_size.png",900,900)
pdf("~/git/noncoding_rvat/figures/revision/png/sample_size.pdf",10,10)

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


options(scipen = 1)
ggplot() + 
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

dev.off()
