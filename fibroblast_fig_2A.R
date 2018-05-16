library(reshape2)
library(ggplot2)
library(ggsignif)
source("common_functions.R")

##################################################################################
######################## Read data ###############################################
##################################################################################

# TFBM density in H3K4me3 peaks versus response divergence

#read in gene data
data = read.csv(file = "fibroblast_data/fibroblast_data.csv", row.names = 1, header = T)

gene_table = data[data$human_padj<0.01,]
gene_table = gene_table[is.na(gene_table$divergence)==0,]

divq1 = quantile(gene_table$divergence,0.25)
divq3 = quantile(gene_table$divergence,0.75)

gene_table$class = "b_medium"
gene_table[gene_table$divergence<divq1,]$class = "a_low"
gene_table[gene_table$divergence>divq3,]$class = "c_high"

## peaks
gene_table_with_peaks = gene_table[is.na(gene_table$H3K4_length)==0,]

gene_table_with_peaks$TFBM_density    = gene_table_with_peaks$TFBM/gene_table_with_peaks$H3K4_length

gene_table_with_peaks$class = "b_medium"
gene_table_with_peaks[gene_table_with_peaks$divergence<divq1,]$class = "a_low"
gene_table_with_peaks[gene_table_with_peaks$divergence>divq3,]$class = "c_high"


##############################################################
########################### stats ############################
##############################################################


empirical_p_value <- calc_empirical_p(
  gene_table_with_peaks$TFBM_density[gene_table_with_peaks$class=="a_low"],
  gene_table_with_peaks$TFBM_density[gene_table_with_peaks$class=="c_high"]
)



##############################################################
########################### figure ###########################
##############################################################




p <- ggplot(gene_table_with_peaks[gene_table_with_peaks$class%in%c("a_low","b_medium","c_high"),], 
       aes(x = class, y = TFBM_density,                                                                                                                   
           fill = class
       )) +
  geom_violin()+
  coord_cartesian(ylim = c(0.0, 0.5)) +
  scale_fill_manual(values=c("#6699DD", "#9966CC", "#CC6677")) +
  geom_signif(comparisons=list(c("a_low","c_high")), annotations=formatSignificance(empirical_p_value),
              y_position = 0.45, tip_length = 0, vjust=0.4) +
  labs(x="Response divergence",
       y="TFBM density")  +
  scale_x_discrete(labels=c("Low","Medium","High"),
                   breaks=c("a_low","b_medium","c_high")) +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = (0),colour="black", size = 10, face="bold"))+
  theme(axis.text.y = element_text(angle = (0),colour="black", size = 10)) +
  theme(axis.title.x = element_text(colour="black", size = 10)) +
  theme(axis.title.y = element_text(colour="black", size = 10))

pdf("output_figures/fig_2A.pdf", width = 2.5, height = 3,useDingbats = F)
print(p)
dev.off()

