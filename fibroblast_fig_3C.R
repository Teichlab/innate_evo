library(ggplot2)
library(gridExtra)
source("common_functions.R")

##################################################################################
######################## Read data ###############################################
##################################################################################

gene_table = read.csv(file = "fibroblast_data/fibroblast_data.csv", row.names = 1, header = T)
gene_table = gene_table[is.na(gene_table$divergence)==0,]
gene_table = gene_table[gene_table$human_padj<0.01,]

#Set categories
gene_table$funct = "x_rest"
gene_table[gene_table$cytokines_and_chemokines==1 | gene_table$cytokines_receptors,]$funct="cyto"
gene_table[gene_table$transcription_factors==1,]$funct="tf"
gene_table[gene_table$kinases_and_phosphatases==1 & gene_table$cytokines_and_chemokines==0,]$funct="kin"

# putting categories in the desired order in the plot
gene_table$funct = factor(gene_table$funct, levels = c("x_rest", "kin", "tf", "cyto"))
gene_table = gene_table[order(gene_table$funct),]

names_gene_table = gene_table[gene_table$name %in% c("IFNB1", "CXCL10", "IKBKE", "CXCL11","CCL5","CCL2","IL6","IL7","IL16","IL1A",
                             "CXCL9","LIF","CCL20","CCL7","IL23A","IL15","CSLP"),]


##################################################################################
######################## SCATTER PLOT OF DM values VS divergence #################
##################################################################################



# top side bar
plot_top <- ggplot(subset(gene_table, funct!="x_rest"), aes(x = divergence, fill = funct)) + geom_density(alpha = 0.75) + 
  scale_fill_manual(values = c("#808000","#008080","#800080"))+
  scale_y_continuous(labels = function(x) sprintf("%.2f", x))+
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        plot.margin = unit(rep(0.1, 4), "cm"))



# scatter 
plot_scatter <- ggplot(gene_table, aes(x = divergence, y = DM_pIC4, colour = funct))+
  guides(colour = guide_legend(nrow = 3))+
  geom_point(size = 2.8, alpha=0.85)+
  geom_text(data = names_gene_table, mapping = aes(x = divergence, y = DM_pIC4, 
                                           label = name, colour = funct, hjust = 1.25, vjust = 1.0),
            angle = 30, size = 3)+
  guides(colour = guide_legend(title = "Function", nrow = 2, ncol = 2))+
  theme_bw()+
  scale_colour_manual(values = c("#c4cccd","#808000","#008080","#800080"))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(rep(0, 4), "cm")) +
  labs(y="Cell-to-cell variability") +
  theme(axis.title.y = element_text(colour="black", size = 15)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  labs(x="Response divergence") +
  theme(axis.text.x = element_text(colour="black", size = 12)) +
  theme(axis.title.x = element_text(colour="black", size = 15)) +
  theme(legend.position = "none")




# right side bar
plot_right <- ggplot(subset(gene_table, funct!="x_rest"), aes(x = DM_pIC4, fill = funct)) + 
  geom_density(alpha = 0.75) + coord_flip() + 
  scale_fill_manual(values = c("#808000","#008080","#800080"))+
  scale_y_continuous(labels = function(x) sprintf("%.2f", x))+
  theme_classic() + 
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        plot.margin = unit(rep(0.1, 4), "cm"))




pdf(file="output_figures/fig_3C.pdf", useDingbats = F)
cowplot::plot_grid(plot_top, NULL,
                   plot_scatter, plot_right, nrow = 2, ncol = 2,
                   rel_heights = c(0.3, 1), rel_widths = c(1, 0.3),
                   align = "hv")
dev.off()
