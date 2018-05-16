source("common_functions.R")

##################################################################################
######################## Read data ###############################################
##################################################################################

gene_table <- read.table("phagocyte_data/phagocyte_functions.csv", 
                         header=T, row.names = 1, stringsAsFactors = FALSE)
genecatlist <- unique(gene_table$class)
gene_table <- gene_table[!is.na(gene_table$divergence),]

genecat_info <- read.csv("phagocyte_data/genecat_info.csv", stringsAsFactors = FALSE)
rownames(genecat_info) <- genecat_info$class


##################################################################################
######################## Statistics ##############################################
##################################################################################


#####################################
## Calculate p-values
catsig <- c()
for(n in genecatlist){
  thep <- calc_empirical_p(
    gene_table$divergence[gene_table$class=="b_DE"],
    gene_table$divergence[gene_table$class==n]
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
  qrt75 <- tapply(gene_table$divergence, as.character(gene_table$class), function(x) quantile(x, probs=0.75, na.rm=TRUE))
  qrt75 <- qrt75[genecatlist]
  
  qrt95 <- tapply(gene_table$divergence, as.character(gene_table$class), function(x) quantile(x, probs=0.95, na.rm=TRUE))
  qrt95 <- qrt95[genecatlist]
  
  
  #Extract titles and colors to use
  col_fill <- genecat_info[genecatlist,]$color
  
  barplot_gg <- ggplot(gtable_df, 
                       aes(
                         x = class, 
                         y = divergence, 
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


barplot_gg <- dotheplot(c("a_all","b_DE","c_viral_sensors","d_cytokines_and_receptors",
                          "g_transcription_factors","h_chromatin_modulators_wide",
                          "i_kinases_and_phosphatases","f_ligases_and_deubiquitinases",
                          "j_enzymes"))
pdf(sprintf("output_figures/fig_2D_phago.pdf"),  width=6.0, height=3.5, useDingbats = FALSE)
print(barplot_gg)
dev.off()


