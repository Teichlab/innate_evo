library(ggplot2)
library(gridExtra)
source("common_functions.R")


##################################################################################
######################## Read data ###############################################
##################################################################################

gene_table = read.csv(file = "fibroblast_data/fibroblast_data.csv", row.names = 1, header = T)
gene_table = gene_table[is.na(gene_table$divergence)==0,]
gene_table = gene_table[gene_table$human_padj<0.01,]

gene_table$is_cyto   <- gene_table$cytokines_and_chemokines==1 | gene_table$cytokines_receptors==1
gene_table$is_tf     <- gene_table$transcription_factors==1
gene_table$is_kinase <- gene_table$kinases_and_phosphatases==1 & gene_table$cytokines_and_chemokines==0

gene_table$funct = "x_rest"
gene_table[gene_table$is_cyto,]$funct   = "cyto"
gene_table[gene_table$is_tf,]$funct     = "tf"
gene_table[gene_table$is_kinase,]$funct = "kin"
gene_table$funct <- factor(gene_table$funct, levels=c("x_rest","kin","tf","cyto"))


##################################################################################
######################## Statistics ##############################################
##################################################################################

list_timepoints <- c("pIC2","pIC4","pIC8")

#for pIC2, pIC4, pIC8, loop and calculate p-values
pout <- NULL
for(curtp in list_timepoints){
  curcond <- sprintf("DM_%s", curtp)
  thispval <- data.frame(
    tp=curtp,
    p_cyto_tf=calc_empirical_p(
      na.omit(gene_table[gene_table$is_cyto,curcond]),
      na.omit(gene_table[gene_table$is_tf,curcond])
    ),
    p_cyto_kinases=calc_empirical_p(
      na.omit(gene_table[gene_table$is_cyto,curcond]),
      na.omit(gene_table[gene_table$is_kinase,curcond])
    ))

  pout <- rbind(pout, thispval)
}
write.table(pout, "fibroblast_out/pval_3D.csv")



##################################################################################
######################## Plot ####################################################
##################################################################################



violin_df <- gene_table[
  gene_table$funct %!in% "x_rest",
  c("funct",sprintf("DM_%s",list_timepoints))]
colnames(violin_df) = c("funct", "dsRNA 2h", "dsRNA 4h", "dsRNA 8h")
violin_df = reshape2::melt(violin_df)
violin_df = violin_df[complete.cases(violin_df),]

g1 = ggplot(violin_df, 
            aes(
              y = value, 
              x = variable, 
              fill = funct))+
  geom_violin(position=position_dodge(width=0.85))+
  scale_fill_manual(values = c("#808000","#008080","#800080"))+
  scale_y_continuous(name = "", breaks = c(0, 1, 2))+
  scale_x_discrete(name = "")+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  theme(legend.position = "none")

pdf(file="output_figures/fig_3D.pdf", useDingbats = F)
print(g1)
dev.off()

