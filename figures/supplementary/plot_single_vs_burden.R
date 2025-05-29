
library(data.table)
library(ggplot2)

fdr5Data = data.table(values = c(1000,406-107,107,410,554-175,175), annotation = c("NC","NC","NC","CDS","CDS","CDS"), association = c("Single","Group","Group-Single","Single","Group","Group-Single"))
p1e9Data = data.table(values = c(290,90-86,86,81,198-162,162), annotation = c("NC","NC","NC","CDS","CDS","CDS"), association = c("Single","Group","Group-Single","Single","Group","Group-Single"))

fdr5NCSum = sum(fdr5Data[1:3,]$values)
fdr5CDSSum = sum(fdr5Data[4:6,]$values)
p1e9NCSum = sum(p1e9Data[1:3,]$values)
p1e9CDSSum = sum(p1e9Data[4:6,]$values)

fdr5Data
p1e9Data

png("fig4a_single_fdr5.png",800,700)
ggplot(data = fdr5Data, aes(x = annotation, y = values, fill = association)) + 
  geom_bar(stat = "identity", color = "black", size = 1.5) +
  annotate(geom = "text", x = "CDS", y = fdr5Data[annotation == "CDS"][association == "Single"]$values, label = paste0(fdr5Data[annotation == "CDS"][association == "Single"]$values," (", round(100 * fdr5Data[annotation == "CDS"][association == "Single"]$values/fdr5CDSSum),"%)"), vjust = 2, size = 8) +
  annotate(geom = "text", x = "CDS", y = fdr5Data[annotation == "CDS"][association == "Single"]$values + 180, label = paste0(fdr5Data[annotation == "CDS"][association == "Group-Single"]$values," (", round(100 * fdr5Data[annotation == "CDS"][association == "Group-Single"]$values/fdr5CDSSum),"%)"), vjust = 2, size = 8) +
  annotate(geom = "text", x = "CDS", y = fdr5Data[annotation == "CDS"][association == "Single"]$values + 500, label = paste0(fdr5Data[annotation == "CDS"][association == "Group"]$values," (", round(100 * fdr5Data[annotation == "CDS"][association == "Group"]$values/fdr5CDSSum),"%)"), vjust = 2, size = 8) +
  annotate(geom = "text", x = "NC", y = fdr5Data[annotation == "NC"][association == "Single"]$values, label = paste0(fdr5Data[annotation == "NC"][association == "Single"]$values," (", round(100 * fdr5Data[annotation == "NC"][association == "Single"]$values/fdr5NCSum),"%)"), vjust = 2, size = 8) +
  annotate(geom = "text", x = "NC", y = fdr5Data[annotation == "NC"][association == "Single"]$values + 100, label = paste0(fdr5Data[annotation == "NC"][association == "Group-Single"]$values," (", round(100 * fdr5Data[annotation == "NC"][association == "Group-Single"]$values/fdr5NCSum),"%)"), vjust = 2, size = 8) +
  annotate(geom = "text", x = "NC", y = fdr5Data[annotation == "NC"][association == "Single"]$values + 350, label = paste0(fdr5Data[annotation == "NC"][association == "Group"]$values," (", round(100 * fdr5Data[annotation == "NC"][association == "Group"]$values/fdr5NCSum),"%)"), vjust = 2, size = 8) +
  # geom_text(aes(label = values)) +
  scale_fill_manual(values = c("#238b45","#74c476","#fc9272")) +
  ylab("gene-trait associations") +
  xlab("Annotation") +
  theme_minimal() + 
  theme(text = element_text(size=32), legend.position = c(0.3,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

dev.off()


png("fig4a_single_p1e9.png",800,700)
ggplot(data = p1e9Data, aes(x = annotation, y = values, fill = association)) + 
  geom_bar(stat = "identity", color = "black", size = 1.5) +
  annotate(geom = "text", x = "CDS", y = p1e9Data[annotation == "CDS"][association == "Single"]$values, label = paste0(p1e9Data[annotation == "CDS"][association == "Single"]$values," (", round(100 * p1e9Data[annotation == "CDS"][association == "Single"]$values/p1e9CDSSum),"%)"), vjust = 3, size = 8) +
  annotate(geom = "text", x = "CDS", y = p1e9Data[annotation == "CDS"][association == "Single"]$values + 150, label = paste0(p1e9Data[annotation == "CDS"][association == "Group-Single"]$values," (", round(100 * p1e9Data[annotation == "CDS"][association == "Group-Single"]$values/p1e9CDSSum),"%)"), vjust = 2, size = 8) +
  annotate(geom = "text", x = "CDS", y = p1e9Data[annotation == "CDS"][association == "Single"]$values + 200, label = paste0(p1e9Data[annotation == "CDS"][association == "Group"]$values," (", round(100 * p1e9Data[annotation == "CDS"][association == "Group"]$values/p1e9CDSSum),"%)"), vjust = 2, size = 8) +
  annotate(geom = "text", x = "NC", y = p1e9Data[annotation == "NC"][association == "Single"]$values, label = paste0(p1e9Data[annotation == "NC"][association == "Single"]$values," (", round(100 * p1e9Data[annotation == "NC"][association == "Single"]$values/p1e9NCSum),"%)"), vjust = 3, size = 8) +
  annotate(geom = "text", x = "NC", y = p1e9Data[annotation == "NC"][association == "Single"]$values + 80, label = paste0(p1e9Data[annotation == "NC"][association == "Group-Single"]$values," (", round(100 * p1e9Data[annotation == "NC"][association == "Group-Single"]$values/p1e9NCSum),"%)"), vjust = 2, size = 8) +
  annotate(geom = "text", x = "NC", y = p1e9Data[annotation == "NC"][association == "Single"]$values + 120, label = paste0(p1e9Data[annotation == "NC"][association == "Group"]$values," (", round(100 * p1e9Data[annotation == "NC"][association == "Group"]$values/p1e9NCSum),"%)"), vjust = 2, size = 8) +
  # geom_text(aes(label = values)) +
  scale_fill_manual(values = c("#238b45","#74c476","#fc9272")) +
  ylab("gene-trait associations") +
  xlab("Annotation") +
  theme_minimal() + 
  theme(text = element_text(size=32), legend.position = c(0.3,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

dev.off()
