library(reshape2)
library(ggplot2)
library(ggsignif)
source("common_functions.R")

##################################################################################
######################## Read data ###############################################
##################################################################################


gene_table = read.csv(file = "fibroblast_data/fibroblast_data.csv", row.names = 1, header = T)
gene_table = gene_table[is.na(gene_table$divergence)==0 & is.na(gene_table$dNdS_29M)==0,]
gene_table = gene_table[gene_table$human_padj<0.01,]

divq1 = quantile(gene_table$divergence,0.25)
divq3 = quantile(gene_table$divergence,0.75)


# prepare for fancy boxplots:
gene_table$divergence_group = "b_middle_divergence"
gene_table[gene_table$divergence<divq1,]$divergence_group <- "a_low_divergence"
gene_table[gene_table$divergence>divq3,]$divergence_group <- "c_high_divergence"

plots_gg <- list()

##################################################################################
######################## Plot dNdS, A ############################################
##################################################################################

p_dnds <- calc_empirical_p(
  gene_table$dNdS_29M[gene_table$divergence_group=="a_low_divergence"],
  gene_table$dNdS_29M[gene_table$divergence_group=="c_high_divergence"]
)


plots_gg[["boxplot_dNdS"]] = 
  ggplot(gene_table[gene_table$divergence_group %in% c("a_low_divergence","b_middle_divergence","c_high_divergence"),], 
                                    aes(x = divergence_group, y = dNdS_29M,fill = divergence_group)) +
  geom_boxplot(outlier.colour = NA) + # option2 - without outliers
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_fill_manual(values=c("#6699DD", "#9966CC", "#CC6677")) +
  geom_signif(comparisons=list(c("a_low_divergence","c_high_divergence")), annotations=formatSignificance(p_dnds),
              y_position = .75, tip_length = 0, vjust=0.4) +
  labs(x="Response divergence",
       y="Sequence evolution (dN/dS)") +
  scale_x_discrete(labels=c("Low","Medium","High"),
                   breaks=c("a_low_divergence","b_middle_divergence","c_high_divergence")) +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = (0),colour="black", size = 10, face="bold"))+
  theme(axis.text.y = element_text(angle = (0),colour="black", size = 7)) +
  theme(axis.title.x = element_text(colour="black", size = 10)) +
  theme(axis.title.y = element_text(colour="black", size = 10))


pdf("output_figures/fig_4A.pdf", width = 2.5, height = 3,useDingbats = F)
print(plots_gg[["boxplot_dNdS"]])
dev.off()


##################################################################################
######################## Plot PPI, D #############################################
##################################################################################


p_ppi <- calc_empirical_p(
  log10(gene_table$cellular_PPIs[gene_table$divergence_group=="a_low_divergence"]),
  log10(gene_table$cellular_PPIs[gene_table$divergence_group=="c_high_divergence"])
)

plots_gg[["boxplot_PPIs"]] = 
  ggplot(gene_table[gene_table$divergence_group %in% c("a_low_divergence","b_middle_divergence","c_high_divergence"),],
         aes(
           x = divergence_group, 
           y = log10(cellular_PPIs),
           fill = divergence_group)) +
  geom_boxplot(outlier.colour = NA) + 
  coord_cartesian(ylim = c(0, 4)) +
  scale_fill_manual(values=c("#6699DD", "#9966CC", "#CC6677")) +
  geom_signif(comparisons=list(c("a_low_divergence","c_high_divergence")), annotations=formatSignificance(p_ppi),
              y_position = 3.6, tip_length = 0, vjust=0.4) +
  labs(x="Response divergence",
       y="log10[Protein-Protein Interactions]") +
  scale_x_discrete(labels=c("Low","Medium","High"),
                   breaks=c("a_low_divergence","b_middle_divergence","c_high_divergence")) +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = (0),colour="black", size = 10, face="bold"))+
  theme(axis.text.y = element_text(angle = (0),colour="black", size = 7)) +
  theme(axis.title.x = element_text(colour="black", size = 10)) +
  theme(axis.title.y = element_text(colour="black", size = 10))


pdf("output_figures/fig_4D.pdf", width = 2.5, height = 3,useDingbats = F)
print(plots_gg[["boxplot_PPIs"]])
dev.off()

##################################################################################
######################## Plot Gain and Loss, B ###################################
##################################################################################


p_gainloss <- calc_empirical_p(
  -log10(0.0001+gene_table$gain_loss_rate_pval[gene_table$divergence_group=="a_low_divergence"]),
  -log10(0.0001+gene_table$gain_loss_rate_pval[gene_table$divergence_group=="c_high_divergence"])
)


plots_gg[["boxplot_gain_loss_rate_pval"]] = 
  ggplot(gene_table[gene_table$divergence_group %in% c("a_low_divergence","b_middle_divergence","c_high_divergence"),], 
         aes(
           x = divergence_group, 
           y = -log10(gain_loss_rate_pval+0.0001),
          fill = divergence_group)) +
  geom_boxplot(outlier.colour = NA) + 
  coord_cartesian(ylim = c(0, 4.5)) +
  scale_fill_manual(values=c("#6699DD", "#9966CC", "#CC6677")) +
  geom_signif(comparisons=list(c("a_low_divergence","c_high_divergence")), annotations=formatSignificance(p_gainloss),
              y_position = 4.2, tip_length = 0, vjust=0.4) +
  labs(x="Response divergence",
       y="Rate of gene gain and loss (-logP)") +
  scale_x_discrete(labels=c("Low","Medium","High"),
                   breaks=c("a_low_divergence","b_middle_divergence","c_high_divergence")) +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = (0),colour="black", size = 10, face="bold"))+
  theme(axis.text.y = element_text(angle = (0),colour="black", size = 7)) +
  theme(axis.title.x = element_text(colour="black", size = 10)) +
  theme(axis.title.y = element_text(colour="black", size = 10))


pdf("output_figures/fig_4B.pdf", width = 2.5, height = 3,useDingbats = F)
print(plots_gg[["boxplot_gain_loss_rate_pval"]])
dev.off()


##################################################################################
######################## Plot Age, C #############################################
##################################################################################
 
p_age <- calc_empirical_p(
  na.omit(gene_table$age_PTHR7_wagner[gene_table$divergence_group=="a_low_divergence"]),
  na.omit(gene_table$age_PTHR7_wagner[gene_table$divergence_group=="c_high_divergence"])
)


plots_gg[["boxplot_Age"]] = 
  ggplot(gene_table[gene_table$divergence_group %in% c("a_low_divergence","b_middle_divergence","c_high_divergence"),], 
         aes(
           x = divergence_group, 
           y = age_PTHR7_wagner,
           fill = divergence_group)) +
  geom_boxplot(outlier.colour = NA) + 
  coord_cartesian(ylim = c(0, 17)) +
  scale_fill_manual(values=c("#6699DD", "#9966CC", "#CC6677")) +
  geom_signif(comparisons=list(c("a_low_divergence","c_high_divergence")), annotations=formatSignificance(p_age),
              y_position = 16, tip_length = 0, vjust=0.4) +
  labs(x="Response divergence",
       y="Gene age") +
  scale_x_discrete(labels=c("Low","Medium","High"),
                   breaks=c("a_low_divergence","b_middle_divergence","c_high_divergence")) +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = (0),colour="black", size = 10, face="bold"))+
  theme(axis.text.y = element_text(angle = (0),colour="black", size = 7)) +
  theme(axis.title.x = element_text(colour="black", size = 10)) +
  theme(axis.title.y = element_text(colour="black", size = 10))


pdf("output_figures/fig_4C.pdf", width = 2.5, height = 3,useDingbats = F)
print(plots_gg[["boxplot_Age"]])
dev.off()
