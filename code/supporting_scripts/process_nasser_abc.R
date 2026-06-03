
library(data.table)
data = fread("AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
summary(data$ABC.Score)

cellTypes = data.table(table(data$CellType))
cellTypes

idMapping = fread("gencode_id_mapping.tsv")
idMapping = idMapping[,.(V4,V6)]

data = merge(data,idMapping, by.x = "TargetGene", by.y = "V6")

## FILTER FOR ABC SCORE
data = data[ABC.Score > 0.1]

## All blood celltypes
wantedCellTypes = c("OCI-LY7-ENCODE","dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_2_hour-Garber2017","MM.1S-ENCODE","T-cell-ENCODE","CD14-positive_monocyte-ENCODE","CD56-positive_natural_killer_cells-Roadmap","CD8-positive_alpha-beta_T_cell-ENCODE","CD34-positive_mobilized-Roadmap","Karpas-422-ENCODE","CD3-positive_T_cell-Roadmap","HAP1","dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_1_hour-Garber2017","erythroblast-Corces2016","Jurkat_anti-CD3_PMA_4hr-Engreitz","CD14-positive_monocyte_treated_with_RPMI_4h-Novakovic2016","CD14-positive_monocyte_treated_with_RPMI_1h-Novakovic2016","THP_pmaLPS_ATAC_72h","CD14-positive_monocyte_treated_with_LPS_d6-Novakovic2016","THP_pmaLPS_ATAC_96h","CD14-positive_monocyte_treated_with_RPMI_d1-Novakovic2016","CD8-positive_alpha-beta_T_cell-Corces2016","CD14-positive_monocyte_treated_with_RPMI_d6-Novakovic2016","THP_pmaLPS_ATAC_120h","CD14-positive_monocytes-Roadmap","CD19-positive_B_cell-Roadmap","CD4-positive_helper_T_cell-ENCODE","CD4-positive_helper_T_cell-Corces2016","B_cell-ENCODE","BJAB-Engreitz","BJAB_anti-IgM_anti-CD40_4hr-Engreitz","CD14-positive_monocyte-Novakovic2016","CD14-positive_monocyte_treated_with_BG_1h-Novakovic2016","CD14-positive_monocyte_treated_with_BG_4h-Novakovic2016","CD14-positive_monocyte_treated_with_BG_d1-Novakovic2016","CD14-positive_monocyte_treated_with_BG_d6-Novakovic2016","CD14-positive_monocyte_treated_with_LPS_1h-Novakovic2016","CD14-positive_monocyte_treated_with_LPS_4h-Novakovic2016","CD14-positive_monocyte_treated_with_LPS_d1-Novakovic2016","H9-Roadmap","THP1_LPS_4hr-Engreitz","THP1-Engreitz","THP-1_monocyte-VanBortle2017","Jurkat-Engreitz","THP-1_macrophage-VanBortle2017","K562-Roadmap","megakaryocyte-erythroid_progenitor-Corces2016","dendritic_cell_treated_with_Lipopolysaccharide_0_ng-mL_for_0_hour-Garber2017","THP_pmaLPS_ATAC_48h","THP_pmaLPS_ATAC_24h","dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_30_minute-Garber2017","dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_4_hour-Garber2017","THP_pmaLPS_ATAC_12h","dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_6_hour-Garber2017","THP_pmaLPS_ATAC_6h","THP_pmaLPS_ATAC_2h","THP_pmaLPS_ATAC_1h","THP_pmaLPS_ATAC_0h","natural_killer_cell-Corces2016","U937_LPS_4hr-Engreitz","GM12878-Roadmap")
filteredData = data[CellType %in% wantedCellTypes]
write.table(filteredData[,.(chr,start,end,V4)],"abc_bloodtypes_score0.1_hg19.bed",row.names=F,sep="\t",col.names=F,quote=F)
