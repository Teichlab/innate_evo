
# This script was originally written by JK Kim & JC Marioni, for:
# Single Cell RNA-Sequencing of Pluripotent States Unlocks Modular Transcriptional Variation, Cell Stem Cell 17, 471-485, 2015.
# Aleksandra A Kolodziejczyk, Jong Kyoung Kim, Jason CH Tsang, Tomislav Ilicic, Johan Henriksson, Kedar N Natarajan, Alex C Tuck, Xuefei Gao, Marc Bühler, Pentao Liu, John C Marioni, Sarah A Teichmann



library(biomaRt)
library(GenomicRanges)
library(plyr)
library(ggplot2)
library(fields)
library(DESeq)
library(zoo)


##################################################################################
######################## An example for DM analysis ##############################
######### rat pahgocytes (10X data) with 4 hours of LPS stimulation ##############
##################################################################################




##################################################################################
######################## Helper functions ########################################
##################################################################################


Quantile.loess<- function(Y, X = NULL, 
                          number.of.splits = NULL,
                          window.size = 20,# this is the default - but later we change to 50
                          percent.of.overlap.between.two.windows = NULL,
                          the.distance.between.each.window = NULL,
                          the.quant = .95,
                          window.alignment = c("center"), 
                          window.function = function(x) {quantile(x, the.quant)}, ...)
{
  if(!is.null(number.of.splits)) {window.size <- ceiling(length(Y)/number.of.splits)}
  if(is.null(the.distance.between.each.window)) {the.distance.between.each.window <- window.size}
  if(!is.null(percent.of.overlap.between.two.windows)) 
  {
    the.distance.between.each.window <- window.size * (1-percent.of.overlap.between.two.windows)
  }
  if(!require(zoo))   
  {
    print("zoo is not installed - please install it.")
    install.packages("zoo")
  }
  if(is.null(X)) {X <- index(Y)} 
  zoo.Y <- zoo(x = Y, order.by = X)
  # calculates rolling medians
  new.Y <- rollapply(zoo.Y, width = window.size, 
                     FUN = window.function,
                     by = the.distance.between.each.window,
                     align = "left")
  new.X <- attributes(new.Y)$index  
  # regression with loess function - fitting to rolling medians
  new.Y.loess <- loess(new.Y~new.X, family = "sym",...)$fitted 
  return(list(y = new.Y, x = new.X, y.loess = new.Y.loess))
}


# for genes with 10X data - no normalization by gene length is done  
# DM un-adjusted by length for 10X - used here
DM <- function(meanGenes, CV2Genes, windowSize=50) {
  meanGenesExpressed = meanGenes[meanGenes > 0 & !is.na(CV2Genes) & CV2Genes > 0]
  CV2GenesExpressed = CV2Genes[meanGenes > 0 & !is.na(CV2Genes) & CV2Genes > 0]
  qloess <- Quantile.loess(log10(CV2GenesExpressed), log10(meanGenesExpressed),
                           the.quant=.5, window.size=windowSize, percent.of.overlap.between.two.windows=.5)
  DM <- aaply(seq(1, length(meanGenesExpressed), 1), 1, function(x) {
    mindex<-which.min(abs(qloess[[2]]-log10(meanGenesExpressed[x])))
    as.vector(log10(CV2GenesExpressed[x]) - qloess[[1]][mindex])}, .expand=FALSE, .progress="text"
  )
  names(DM) <- names(meanGenesExpressed)
  DM
}

  

##################################################################################
######################## Calculate DM ############################################
##################################################################################


# The input includes a list of cleaned cells to be used in cell-to-cell variability analysis.
# The cells were cleaned in previous stages (see Methods).
# Similar files for all individuals, conditions and species are in ArrayExpress (Accession: E-MTAB-6754)
gene_expression_matrix =  read.table(file="phagocyte_data/rat2_lps4_filtered_by_cell_cluster0.txt", header = T, row.names = 1)
  
    
# Remove genes that are not expressed in at least 0.5% of the cells
gene_expression_matrix_filtered = gene_expression_matrix[rowSums(gene_expression_matrix>0)>ncol(gene_expression_matrix)*0.005,]
  

# Compute DM
meanGenes = rowMeans(gene_expression_matrix_filtered)
CV2Genes = apply(gene_expression_matrix_filtered , 1, var) / meanGenes^2
names(meanGenes) = rownames(gene_expression_matrix_filtered )
names(CV2Genes) = rownames(gene_expression_matrix_filtered )
DMlevels = DM(meanGenes, CV2Genes) 

# Store results 
out_file =  "output_figures/rat2_lps4_DM.txt"
write.table(DMlevels, file = out_file, quote = F)
