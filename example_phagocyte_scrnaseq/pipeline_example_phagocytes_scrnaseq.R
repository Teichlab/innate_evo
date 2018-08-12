# A pipeline for phagocyte single-cell QC, clustering and DE analysis.
# We here give an example of the rabbit unstimulated cells. Raw input files can be found in ArrayExpress

library(Seurat)
library(dplyr)
library(Matrix)
library(sqldf)




# Function to read 10X data and create an R-object with ENSEMBL IDs 
Read10XbyENSEMBLid <- function(data.dir = NULL){
  full.data <- list()
  for (i in seq_along(data.dir)) {
    run <- data.dir[i]
    if (! dir.exists(run)){
      stop("Directory provided does not exist")
    }
    if(!grepl("\\/$", run)){
      run <- paste(run, "/", sep = "")
    }
    barcode.loc <- paste0(run, "barcodes.tsv")
    gene.loc <- paste0(run, "genes.tsv")
    matrix.loc <- paste0(run, "matrix.mtx")
    if (!file.exists(barcode.loc)){
      stop("Barcode file missing")
    }
    if (! file.exists(gene.loc)){
      stop("Gene name file missing")
    }
    if (! file.exists(matrix.loc)){
      stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(
        x = as.character(
          x = sapply(
            X = cell.names,
            FUN = ExtractField, field = 1,
            delim = "-"
          )
        )
      )
    }
    rownames(x = data) <- make.unique(
      names = as.character(
        x = sapply(
          X = gene.names,
          FUN = ExtractField,
          field = 1,
          delim = "\\t"
        )
      )
    )
    if (is.null(x = names(x = data.dir))) {
      if(i < 2){
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_",
                                   cell.names)
    }
    full.data <- append(x = full.data, values = data)
  }
  full.data <- do.call(cbind, full.data)
  return(full.data)
}


# tables of coding and mitochondrial genes (can be found in this directory)
coding_genes = read.table("Oryctolagus_cuniculus.OryCun2.0.84_coding_with_names_and_transcripts.txt", header = F, row.names = 1)
MT_genes = read.table("Oryctolagus_cuniculus.OryCun2.0.84_MT_genes.txt", header = F, row.names = 1)

unst.data = Read10XbyENSEMBLid(data.dir = "rabbit1_unst_25220_8/")



unst <- CreateSeuratObject(raw.data = unst.data, min.cells = 0, min.genes = 200,
                           project = "unst")
unst_matrix = as.matrix(unst@raw.data) # convert to sparse matrix

# get MT data
unst_matrix_MT = unst_matrix[rownames(MT_genes),]

# subset to include only coding genes  
unst_matrix_coding = unst_matrix[rownames(coding_genes),]
unst_matrix_coding_filtered = unst_matrix_coding


################### QC for each individual library


# MT fraction
unst_matrix_MT = unst_matrix_coding_filtered[rownames(MT_genes),]

sum_MT1 = colSums(unst_matrix_MT)
sum_all = colSums(unst_matrix_coding_filtered) 

percent.mito = sum_MT1/sum_all*100

### Metadata table
metadata_unst = data.frame(row.names = colnames(unst_matrix_coding_filtered),
                           "p_mt" = percent.mito,
                           "n_genes" = colSums(unst_matrix_coding_filtered>0),
                           "total_reads" = colSums(unst_matrix_coding_filtered))

# (2) remove cells by high fraction of MT or by low number of UMI

good_cells = rownames(metadata_unst[metadata_unst$p_mt<=10 &
                                      metadata_unst$n_genes<=quantile(metadata_unst$n_genes, 9/10) &
                                      metadata_unst$n_genes>=200,
                                    ])

metadata_unst_filtered = metadata_unst[good_cells,]

unst_matrix_filtered = unst_matrix_coding_filtered[,good_cells]

write.table(x = unst_matrix_filtered, file = "rabbit1_unst_filtered_by_cell.txt", quote = F)


###############################################
############################################### analyze and cluster
###############################################


### merge 3 rabbit UNST - these are the inputs from the previous stage, and QC them
unst_matrix_filtered_cleaned_rabbit1 = read.table(file = "rabbit1_unst_filtered_by_cell.txt", row.names = 1, header = 1)
unst_matrix_filtered_cleaned_rabbit2 = read.table(file = "rabbit2_unst_filtered_by_cell.txt", row.names = 1, header = 1)
unst_matrix_filtered_cleaned_rabbit3 = read.table(file = "rabbit3_unst_filtered_by_cell.txt", row.names = 1, header = 1)

merged_table_rabbit_raw = Reduce(mergeAndFormat, list(unst_matrix_filtered_cleaned_rabbit1, 
                                                      unst_matrix_filtered_cleaned_rabbit2,
                                                      unst_matrix_filtered_cleaned_rabbit3))

# we remove genes that are not expressed in at least 0.5% of the ENTIRE cell population (across 3 batches)
merged_table_rabbit = merged_table_rabbit_raw[rowSums(merged_table_rabbit_raw>0)>ncol(merged_table_rabbit_raw)*0.005,]

merged_matrix_MT = merged_table_rabbit[rownames(MT_genes),]

sum_MT1 = colSums(merged_matrix_MT)
sum_all = colSums(merged_table_rabbit) 
percent.mito = sum_MT1/sum_all*100

metadata_merged = data.frame(row.names = colnames(merged_table_rabbit),
                             "p_mt" = percent.mito,
                             "n_genes" = colSums(merged_table_rabbit>0),
                             "total_reads" = colSums(merged_table_rabbit),
                             "individual" = c(rep("rabbit1", ncol(unst_matrix_filtered_cleaned_rabbit1)),
                                              rep("rabbit2",ncol(unst_matrix_filtered_cleaned_rabbit2)),
                                              rep("rabbit3",ncol(unst_matrix_filtered_cleaned_rabbit3)))
)

# remove cells that don't have a minimum number of genes or that their MT fraction is too high (>10%)
quant10 = quantile(metadata_merged$n_genes, c(1,9)/10)
quant10 = as.numeric(quant10[1])

good_cells1 = rownames(metadata_merged[metadata_merged$p_mt<=10 & metadata_merged$n_genes>=quant10,]) ## lower 10%

metadata_merged_filtered1 = metadata_merged[good_cells1,]
rabbit_merged_filtered1 = merged_table_rabbit[,good_cells1]

# normalize data
merged_table_norm1 = limma::removeBatchEffect(log2(rabbit_merged_filtered1+1), 
                                              batch = metadata_merged_filtered1$individual,
                                              covariates = metadata_merged_filtered1$total_reads)

## identify clusters with Seurat's functions
merged_table_obj = CreateSeuratObject(merged_table_norm1, project = "merged")
merged_table_obj@scale.data = merged_table_obj@data
merged_table_obj = FindVariableGenes(merged_table_obj)

# run PCA
merged_table_obj = RunPCA(merged_table_obj, pc.genes = row.names(merged_table_obj@data))

# find clusters, based on 20 PCs and with a resolution of 0.1
seurat_clusters01 = FindClusters(merged_table_obj, dims.use = 1:20, resolution = 0.1)
clusters_data01 = as.matrix(seurat_clusters01@ident)
clusters_data_df01 = as.data.frame(seurat_clusters01@ident, rownames = 1)

write.table(clusters_data_df03, file = "rabbit_all3_unst_cells_to_clusters_01.txt", quote = F, row.names = T)



###############################################
############################################### 
###############################################


# Subset original input based on cleaned cells from the first cluster (code not given)


###############################################
############################################### 
###############################################


# merge files
mergeAndFormat = function(...){ 
  x = merge(..., all=T, by = 0)
  rownames(x) = x[,1]
  x = x[,-1]
  return(x)
}


#read condition 1
rabbit1_unst = read.table(file = "rabbit1_unst_filtered_by_cell_cluster0.txt", row.names = 1, header = T)
rabbit2_unst = read.table(file = "rabbit2_unst_filtered_by_cell_cluster0.txt", row.names = 1, header = T)
rabbit3_unst = read.table(file = "rabbit3_unst_filtered_by_cell_cluster0.txt", row.names = 1, header = T)

# read condition 2
rabbit1_lps4 = read.table(file = "rabbit1_lps4_filtered_by_cell_cluster0.txt", row.names = 1, header = T)
rabbit2_lps4 = read.table(file = "rabbit2_lps4_filtered_by_cell_cluster0.txt", row.names = 1, header = T)
rabbit3_lps4 = read.table(file = "rabbit3_lps4_filtered_by_cell_cluster0.txt", row.names = 1, header = T)


# these input files should include only high-quality cells - from specific clusters that are analogous across conditions

# merge data and create metadata 
all_merged = Reduce(mergeAndFormat, list(rabbit1_unst, rabbit2_unst,rabbit3_unst, 
                                         rabbit1_lps4, rabbit2_lps4, rabbit3_lps4))

treatments <- factor(c(rep("UNST",ncol(rabbit1_unst)),
                       rep("UNST",ncol(rabbit2_unst)),
                       rep("UNST",ncol(rabbit3_unst)),
                       rep("lps4",ncol(rabbit1_lps4)),
                       rep("lps4",ncol(rabbit2_lps4)),
                       rep("lps4",ncol(rabbit3_lps4))
              ),levels=c("UNST","lps4"))

metadata_merged = data.frame(row.names = colnames(all_merged),
                             "n_genes" = colSums(all_merged>0),
                             "total_reads" = colSums(all_merged),
                             "treatment" = treatments,
                             "individual" = c(
                               rep("rabbit1",ncol(rabbit1_unst)),
                               rep("rabbit2",ncol(rabbit2_unst)),
                               rep("rabbit3",ncol(rabbit3_unst)),
                               rep("rabbit1",ncol(rabbit1_lps4)),
                               rep("rabbit2",ncol(rabbit2_lps4)),
                               rep("rabbit3",ncol(rabbit3_lps4))
                             )
)

metadata_merged$experiment <- sprintf("%s_%s", metadata_merged$treatment,metadata_merged$individual)

# Create "in vitro" bulk objects by combining cells
listexp <- sort(unique(metadata_merged$experiment)) # list of individual experiments
countbulk <- matrix(nrow=nrow(all_merged), ncol=length(listexp)) # per experiment, get a count matrix
rownames(countbulk) <- rownames(all_merged)
for(i in 1:length(listexp)){ 
  countbulk[,i] <- apply(all_merged[,metadata_merged$experiment==listexp[i]],1,sum) # sum up all cells per experiment
}
condbulk <- sqldf("select distinct experiment, treatment, individual from metadata_merged order by experiment") # design matrix for DESeq

####### Perform DE test using DEseq2 (wald test)
dds_bulk = DESeq2::DESeqDataSetFromMatrix(countData = countbulk, 
                                          colData = condbulk,
                                          design = ~individual+treatment) 
dds_bulk = DESeq2::DESeq(dds_bulk) 
result_bulk <- DESeq2::results(dds_bulk, contrast=c("treatment","UNST","lps4")) 
result_bulk <- result_bulk[order(result_bulk$padj),]

write.table(result_bulk, file = "rabbit1_2_3_bulk_UNST_vs_lps4_DE_Deseq2.txt",quote = F)
