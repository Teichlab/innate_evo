library(reshape2)
library(ggplot2)

##################################################################################
######################## Read data ###############################################
##################################################################################


# we take only the NAs between -1000 and +500
gene_low_500  = read.table("fibroblast_data/phylopromoters/human_DE_low.txt", header=T, row.names = 1)
gene_med_500  =  read.table("fibroblast_data/phylopromoters/human_DE_mid.txt", header=T, row.names = 1)
gene_high_500 =  read.table("fibroblast_data/phylopromoters/human_DE_high.txt", header=T, row.names = 1)

##################################################################################
######################## Statistics ##############################################
##################################################################################


### calculate means
mean_low_500  = apply(gene_low_500, 2, mean)
mean_med_500  = apply(gene_med_500, 2, mean)
mean_high_500 = apply(gene_high_500, 2, mean)


# stats
ks.test(mean_low_500, mean_high_500, alternative = "greater")
#D^+ = 0.132, p-value = 0.0001646


##################################################################################
######################## Plot figures ############################################
##################################################################################

test1 = data.frame("values" = c(colMeans(gene_low_500), colMeans(gene_med_500), colMeans(gene_high_500)),
                   "position" = rep(-500:-1, 3),
                   "classes" = rep(c("low divergence", "medium divergence", "high divergence"), each = 500))


plot_mean = ggplot(test1, 
                   aes(
                     x = position, 
                     y = values, 
                     group = classes, 
                     colour = classes))+
  #geom_point()
  geom_smooth()+
  # we order alphabetically high, then low, then med
  scale_colour_manual(values =  c( "#CC6688","#1199EE", "#995599") )+
  scale_x_continuous(name = "Sequence position", breaks = c(-500, -250, 0)) + 
  scale_y_continuous(name = "Sequence conservation (phyloP)", breaks = c(0.1, 0.15, 0.2)) + 
  theme_classic() + 
  theme(legend.position = "none")

pdf(file="output_figures/fig_2B.pdf",width=5,height=4, useDingbats = FALSE)
print(plot_mean)
dev.off()


