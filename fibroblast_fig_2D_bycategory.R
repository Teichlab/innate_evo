##################################################################################
######################## Read data ###############################################
##################################################################################

gene_table <- read.table("fibroblast_data/fibroblast_input.csv", 
                         header=T, row.names = NULL, stringsAsFactors = FALSE)
genecatlist <- unique(gene_table$class)
gene_table <- gene_table[!is.na(gene_table$divergence_all_simple),]

genecat_info <- read.csv("fibroblast_data/genecat_info.csv", stringsAsFactors = FALSE)
rownames(genecat_info) <- genecat_info$class


##################################################################################
######################## Statistics ##############################################
##################################################################################


#####################################
## Calculate p-values
catsig <- c()
for(n in genecatlist){
  if(n=="b_DE")
    thep <- 0.5
  else
    thep <- calc_empirical_p(
      gene_table$divergence_all_simple[gene_table$class=="b_DE"],
      gene_table$divergence_all_simple[gene_table$class==n]
    )
  thep<-p.adjust(thep, method = 'BH')
  catsig <- c(catsig, thep)
}
catsigFormatted <- unlist(lapply(catsig,formatSignificance))



##################################################################################
######################## Plotting ################################################
##################################################################################




dotheplot <- function(cats){
  #Which categories to include?
  #includecat <- rep(TRUE, length(genecatlist))
  #includecat <- genecatlist %in% genecatlist[grep("DE",genecatlist)]
  includecat <- genecatlist %in% cats

  #Print P-values
  print(catsig[includecat])
  
  #How many categories?
  num_cat <- sum(includecat)
  
  #Reduce data table to only include these categories
  gtable_df <- gene_table[as.character(gene_table$class) %in% genecatlist[includecat],]
  
  #####################################
  # Calculate quantiles
  qrt75 <- tapply(gene_table$divergence_all_simple, as.character(gene_table$class), function(x) quantile(x, probs=0.75, na.rm=TRUE))
  qrt75 <- qrt75[genecatlist]
  
  qrt95 <- tapply(gene_table$divergence_all_simple, as.character(gene_table$class), function(x) quantile(x, probs=0.95, na.rm=TRUE))
  qrt95 <- qrt95[genecatlist]
  
  
  #Extract titles and colors to use
  col_fill <- genecat_info[genecatlist,]$color
  
  barplot_gg <- ggplot(gtable_df, 
                       aes(
                         x = class, 
                         y = divergence_all_simple, 
                         fill = class)) +
    geom_boxplot(width=.7, outlier.colour = NA) + 
    #coord_cartesian(ylim = c(-9, 1.0)) +
    scale_fill_manual(values=col_fill[includecat]) +
    
    geom_signif(
      annotations=catsigFormatted[includecat],
      color=rep(col_fill[includecat],3),
      y_position = max(qrt95[includecat]), tip_length = 0, vjust=0.4,
      xmin=c(1:num_cat),xmax=c(1:num_cat)) +
    
    labs(y="Response divergence", x=NULL) +
    scale_x_discrete(
      labels=genecat_info[genecatlist[includecat],]$prettyname,
      breaks=genecatlist[includecat]) +
    theme_classic() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(
      angle = (-45),colour=col_fill[includecat],
      size = 7, face="bold",hjust=0, vjust=1))+
    theme(axis.text.y = element_text(angle = (0),colour="black", size = 7)) +
    theme(axis.title.x = element_text(colour="black", size = 10)) +
    theme(axis.title.y = element_text(colour="black", size = 10))
  return(barplot_gg)
}


barplot_gg <- dotheplot(c("a_all","b_DE","d_viral_sensors_DE","g_cytokine_receptors_DE",
                          "j_transcription_factors_DE","l_remodeling_DE","n_kinases_DE","p_ligases_DE",
                          "r_enzymes_non_regulatory_DE"))
pdf(sprintf("output_figures/fig_2D.pdf"),  width=6.0, height=3.5, useDingbats = FALSE)
print(barplot_gg)
dev.off()

barplot_gg <- dotheplot(c("a_all","t_self_defense_DE","v_inflammation_DE",
                          "x_apoptosis_DE","z_regulation_union_DE"))
pdf(sprintf("output_figures/fig_S4A.pdf"),  width=6.0, height=3.5, useDingbats = FALSE)
print(barplot_gg)
dev.off()

