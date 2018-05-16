library(ggplot2)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)

##################################################################################
######################## Read data ###############################################
##################################################################################

gene_table = read.csv("fibroblast_data/gene_data_features.csv", header = T, row.names = 1)

#### fix NAs based on column-specific values
gene_table[is.na(gene_table$new_div)==1,]$new_div=(-7) 
gene_table[is.na(gene_table$dNdS_29M)==1,]$dNdS_29M=0 
gene_table[is.na(gene_table$gene_age)==1,]$gene_age=-1 
gene_table[is.na(gene_table$DM_pIC4)==1,]$DM_pIC4=-1 

#### Normalize ranges

linMap <- function(x, from, to) {
  (x - min(x)) / max(x - min(x)) * (to - from) + from
}

gene_table$new_div = linMap(gene_table$new_div,0,100)
gene_table$dNdS_29M = linMap(gene_table$dNdS_29M,0,100)
gene_table$gene_age= linMap(gene_table$gene_age,0,100)
gene_table$viral_PPIs = linMap(gene_table$viral_PPIs,0,100)
gene_table$interactions_binding_with_evidence = linMap(gene_table$interactions_binding_with_evidence,0,100)
gene_table$DM_pIC4 = linMap(gene_table$DM_pIC4,0,100)


##################################################################################
######################## Plotting ################################################
##################################################################################

colnames(gene_table) = c(
  "Response\ndivergence",
  "Cell-to-cell\nvariability",
  "Sequence\nevolution",
  "Gene age",
  "Cellular PPI", 
  "Host-virus PPI")

pdf("output_figures/fig_4F.pdf",width=5,height=7, useDingbats = FALSE)
color_range  = colorRampPalette(brewer.pal(9,"Blues"))(100)
pheatmap(gene_table, 
         color = color_range,
         gaps_row = c(18,18,18,32,32,32),
         gaps_col = c(4,4,4,4),
         show_rownames = T, show_colnames = T, cluster_rows = F, cluster_cols = F)

dev.off()

