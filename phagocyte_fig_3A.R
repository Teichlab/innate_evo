library(reshape2)
library(ggplot2)
library(ggsignif)

##################################################################################
######################## Read data ###############################################
##################################################################################

gene_table = read.csv(file = "phagocyte_data/phagocytes_summary_data.csv", row.names = 1, header = T)
gene_table = gene_table[is.na(gene_table$divergence)==0,]
gene_table = gene_table[!is.na(gene_table$mouse_padj) & gene_table$mouse_padj<0.01,]

gene_table = gene_table[is.na(gene_table$dm_mouse2_lps4)==0,]
divq1 = quantile(gene_table$divergence,0.25)
divq3 = quantile(gene_table$divergence,0.75)

gene_table$divergence_group = "b_middle_divergence"
gene_table[gene_table$divergence>divq3,]$divergence_group <- "c_high_divergence"
gene_table[gene_table$divergence<divq1,]$divergence_group <- "a_low_divergence"

##################################################################################
######################## Main Figure DM at 4hrs, mouse ###########################
##################################################################################



p_dm <- calc_empirical_p(
  na.omit(gene_table$dm_mouse2_lps4[gene_table$divergence_group=="a_low_divergence"]),
  na.omit(gene_table$dm_mouse2_lps4[gene_table$divergence_group=="c_high_divergence"])
)


plots_gg = ggplot(gene_table[gene_table$divergence_group%in%c("a_low_divergence","b_middle_divergence","c_high_divergence"),], 
                  aes(
                    x = divergence_group, 
                    y = dm_mouse2_lps4,
                    fill = divergence_group)) +
  geom_boxplot(outlier.colour = NA) + 
  coord_cartesian(ylim = c(-0.2, 0.6)) +
  scale_fill_manual(values=c("#6699DD", "#9966CC", "#CC6677")) +
  geom_signif(comparisons=list(c("a_low_divergence","c_high_divergence")), annotations=formatSignificance(p_dm),
              y_position = 0.5, tip_length = 0, vjust=0.4) +
  labs(x="Response divergence",
       y="Cell-to-cell Variability (DM)") +
  scale_x_discrete(labels=c("Low","Medium","High"),
                   breaks=c("a_low_divergence","b_middle_divergence","c_high_divergence")) +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = (0),colour="black", size = 10, face="bold"))+
  theme(axis.text.y = element_text(angle = (0),colour="black", size = 7)) +
  theme(axis.title.x = element_text(colour="black", size = 10)) +
  theme(axis.title.y = element_text(colour="black", size = 10))



pdf("output_figures/fig_3A_phago.pdf", width = 2.5, height = 3,useDingbats = F)
print(plots_gg)
dev.off()
