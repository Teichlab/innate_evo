
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
######################## Helper functions ########################################
##################################################################################


# or USE source("http://www.r-statistics.com/wp-content/uploads/2010/04/Quantile.loess_.r.txt")
Quantile.loess<- function(Y, X = NULL, 
                          number.of.splits = NULL,
                          window.size = 20,# this is the default - but later we change to 50
                          percent.of.overlap.between.two.windows = NULL,
                          the.distance.between.each.window = NULL,
                          the.quant = .95,
                          window.alignment = c("center"), 
                          window.function = function(x) {quantile(x, the.quant)}, ...){
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
  new.Y <- rollapply(zoo.Y, width = window.size, 
                     FUN = window.function,
                     by = the.distance.between.each.window,
                    align = "left")
                     #CHANGED FROM ORIGINAL!!:         align = window.alignment)
  new.X <- attributes(new.Y)$index  
  new.Y.loess <- loess(new.Y~new.X, family = "sym",...)$fitted 
  return(list(y = new.Y, x = new.X, y.loess = new.Y.loess))
}


DMLengthAdjusted <- function(ncountAll, cellType, geneLength, windowSize=50, minCount=10) { 
  colIndex = colnames(ncountAll) == cellType
  meanGenes <- rowMeans(ncountAll[,colIndex])
  CV2Genes = apply(as.matrix(ncountAll[,colIndex]), 1, var) / meanGenes^2
  meanGenesExpressed = meanGenes[meanGenes >= minCount & !is.na(CV2Genes) & CV2Genes > 0]
  CV2GenesExpressed = CV2Genes[meanGenes >= minCount & !is.na(CV2Genes) & CV2Genes > 0]
  qloess <- Quantile.loess(log10(CV2GenesExpressed), log10(meanGenesExpressed),
                           the.quant=.5, window.size=windowSize, percent.of.overlap.between.two.windows=.5)
  DM <- aaply(seq(1, length(meanGenesExpressed), 1), 1, function(x) {
    mindex<-which.min(abs(qloess[[2]]-log10(meanGenesExpressed[x])))
    as.vector(log10(CV2GenesExpressed[x]) - qloess[[1]][mindex])}, .expand=FALSE, .progress="text"
  )
  names(DM) <- names(meanGenesExpressed)
  
  DMLength = merge(data.frame(gene=geneLength[,1], legnth=geneLength[,2]), data.frame(gene=names(DM), DM=DM), by.x="gene", by.y="gene")
  qloess <- Quantile.loess(DMLength[,3], log10(DMLength[,2]),
                           the.quant=.5, window.size=windowSize, percent.of.overlap.between.two.windows=.5)
  DM2 <- aaply(seq(1, nrow(DMLength), 1), 1, function(x) {
    mindex<-which.min(abs(qloess[[2]]-log10(DMLength[x,2])))
    as.vector(DMLength[x,3] - qloess[[1]][mindex])}, .expand=FALSE, .progress="text"
  )
  names(DM2) = DMLength[,1]
  DM2
}


##################################################################################
######################## Process data ############################################
##################################################################################

list_fibro_cond <- c("unst","PIC2","PIC4","PIC8")
list_fibro_species <- c("human","mouse","rat","rhesus")

for(curspecies in list_fibro_species){
  for(curcond in list_fibro_cond){
    print(sprintf("%s - %s",curspecies, curcond))

    ##### Prepare the count table
    countGenes = read.table(file=sprintf("fibroblast_data/counts_%s_%s.csv",curspecies,curcond), 
                                 stringsAsFactor=FALSE, header=FALSE, check.names=FALSE, row.names = 1)
    if(rownames(countGenes)[1]=="gene"){
      #Human is a different format
      countGenes = read.table(file=sprintf("fibroblast_data/counts_%s_%s.csv",curspecies,curcond), 
                              stringsAsFactor=FALSE, header=TRUE, check.names=FALSE, row.names = 1)
    }
    
    ##### Include only genes that are not completely zeros
    countGenes = countGenes[rowSums(countGenes)>0,] 
    
    ##### Remove ERCCs
    countGenes = countGenes[!grepl("ERCC", rownames(countGenes)),]
    
    ##### Normalize read counts
    countGenes = round(countGenes)
    ncountGenes = t(t(countGenes) / estimateSizeFactorsForMatrix(countGenes))
    
    ##### Calculate DM
    geneLength = read.table(sprintf("fibroblast_data/genelength_%s.csv",curspecies), sep="\t", stringsAsFactor=FALSE)
    
    meanCutoff = 10
    colnames(ncountGenes) = rep("cell", ncol(ncountGenes))
    DMOS25 = DMLengthAdjusted(ncountGenes, "cell", geneLength, minCount=meanCutoff)
    
    write.table(DMOS25, sprintf("fibroblast_out/dm_%s_%s.csv",curspecies,curcond))
  }
}
