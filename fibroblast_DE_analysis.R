# load libararies
library(limma)
library(edgeR)
library(zoo)
library(lmtest)

##################################################################################
######################## An example for DE analysis ##############################
######### rat fibroblast with poly I:C (dsRNA) stimulation  vs control############
##################################################################################



#########################################
######## get data #######################
#########################################

# count matrix
countsMatrix<-read.table("fibroblast_data/rat_counts_pIC_ctrl.txt", header=T, row.names=1)
countsMatrix<-round(countsMatrix) 

# samples description
info_table<-read.table("fibroblast_data/rat_samples_info.txt", header=T, row.names=1)

info_table$treatment = factor(info_table$treatment, levels = c("LF", "PIC")) # reorganized to get UNST as the intercept; not alphabetiaclly

info_table$treatment_and_time<-factor(paste0(info_table$treatment, info_table$timepoint)) #we add treatment_and_time as a factor combining the two old factors 


#########################################
######## Use edgeR ######################
#########################################

edgeR_DE_obj <- DGEList(counts=countsMatrix, group=info_table$treatment_and_time)
edgeR_DE_obj <- calcNormFactors(edgeR_DE_obj, method="TMM") 
edgeR_DE_obj <- estimateDisp(edgeR_DE_obj)
et <- exactTest(edgeR_DE_obj, pair=c("LF4","PIC4"))
et$table$fdr <- p.adjust(et$table$PValue, "BH")


write.table(et$table, file="output_figures/rat_edgeR_PIC_ctrl.txt")
