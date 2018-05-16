library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggplot2)
library(grid)
library(psych)
source("common_functions.R")

##################################################################################
######################## Read data ###############################################
##################################################################################

gene_table = read.csv(file = "fibroblast_data/fibroblast_data.csv", row.names = 1, header = T)
gene_table = gene_table[is.na(gene_table$divergence)==0,]
gene_table = gene_table[gene_table$human_padj<0.01,]

##################################################################################
######################## Statistics ##############################################
##################################################################################

gene_table$viral_class = "a"
gene_table[gene_table$viral_PPIs>0,]$viral_class = "b"
gene_table[gene_table$viral_immunomodulators_PPI>0,]$viral_class = "c"

pval_ab <- calc_empirical_p(
  gene_table$divergence[gene_table$viral_class=="a"],
  gene_table$divergence[gene_table$viral_class=="b"]
)

pval_bc <- calc_empirical_p(
  gene_table$divergence[gene_table$viral_class=="b"],
  gene_table$divergence[gene_table$viral_class=="c"]
)


##################################################################################
######################## Plots ###################################################
##################################################################################


box_gg = ggplot(gene_table[gene_table$viral_class %in% c("a","b","c"),], 
                aes(
                  x = viral_class, 
                  y = divergence,                                                                                                                     #colour = class,
                  fill = viral_class)) +
  geom_violin() +
  coord_cartesian(ylim = c(-6.5, 4.5)) +
  scale_fill_manual(values=c("#2b4070", "#304ACE", "#5B6BBC")) +
  geom_signif(comparisons=list(c("a","b"),c("b","c")), annotations=c(formatSignificance(pval_ab),formatSignificance(pval_bc)),
              y_position = c(4.1,3.7), tip_length = 0, vjust=0.4) +
  labs(x="",
       y="Response divergence") +
  scale_x_discrete(labels=c(
    "Non-viral\nInteractors",
    "Viral\nInteractors",
    "Immunomodulator\nInteractors"),breaks=c("a","b","c")) +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = (0),colour="black", size = 7, face="bold"))+
  theme(axis.text.y = element_text(angle = (0),colour="black", size = 10)) +
  theme(axis.title.x = element_text(colour="black", size = 2)) +
  theme(axis.title.y = element_text(colour="black", size = 10))

pdf("output_figures/fig_4E.pdf", width = 2.8, height = 3.3,useDingbats = F)
print(box_gg)
dev.off()

