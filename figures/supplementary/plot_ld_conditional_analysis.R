
library(data.table)
library(ggplot2)
library(patchwork)

finalDT = fread("~/burdenNC/data/results/200k/conditional/LD/summary.conditional.pvals.tsv", sep="\t")

png("~/git/noncoding_rvat/figures/png/conditional_ld_all.png",1500,800)

# When doing conditional, most p-values get weaker, as seen previously
100*nrow(finalDT[nocondLOG10 > fullLOG10])/nrow(finalDT)
100*nrow(finalDT[nocondLOG10 <= fullLOG10])/nrow(finalDT)

# # When varians in LD are ommited from conditional, p-value tends to become more significant
# 100*nrow(finalDT[subLOG10 > fullLOG10])/nrow(finalDT)
# 100*nrow(finalDT[subLOG10 <= fullLOG10])/nrow(finalDT)

# Indeed, doing conditional with the variants in LD often suffice to decrease p-value significance
100*nrow(finalDT[nocondLOG10 > intLOG10])/nrow(finalDT)
100*nrow(finalDT[nocondLOG10 <= intLOG10])/nrow(finalDT)

# Perc of FDR 5% associations that remain after conditioning for variants in LD
minLOG10 = min(finalDT$nocondLOG10)
100 * nrow(finalDT[nocondLOG10 > minLOG10][intLOG10 > minLOG10]) / nrow(finalDT[nocondLOG10 > minLOG10])
# Perc of 1e-9 associations that remain after conditioning for variants in LD
100 * nrow(finalDT[nocondLOG10 > 9][intLOG10 > 9]) / nrow(finalDT[nocondLOG10 > 9])
nrow(finalDT[intLOG10 > minLOG10])

#################
# No cond vs LD cond
#################
test = cor.test(finalDT$nocondLOG10,finalDT$intLOG10, method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
g1 = ggplot(data = finalDT, aes(nocondLOG10, intLOG10)) + 
  geom_point() +
  # geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 7, fontface = "bold"  ) +
  labs(tag = "a", title = "LD conditional") +
  ylab("LD conditional (-log10)") +
  xlab("No conditional (-log10)") +
  theme_minimal() + 
  theme(text = element_text(size=32), legend.position = c(0.86,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  

##################
# NO cond vs Full cond
##################
test = cor.test(finalDT$nocondLOG10,finalDT$fullLOG10, method = "spearman")
text = paste("Spearman R:",round(test$estimate,2), "P-value",format.pval(test$p.value,2))
g2 = ggplot(data = finalDT, aes(nocondLOG10, fullLOG10)) + 
  geom_point() +
  # geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed", size = 2) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.5, vjust = 1.5, size = 7, fontface = "bold"  ) +
  ylab("Full conditional (-log10)") +
  xlab("No conditional (-log10)") +
  labs(tag = "b", title = "Full conditional") +
  theme_minimal() + 
  theme(text = element_text(size=32), legend.position = c(0.86,0.85), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  


g1 + g2

dev.off()
