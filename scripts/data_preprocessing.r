library(dplyr)
library(NMF)
library(vcd)
library(survival)
library(ggstatsplot)
library(ggplot2)
library(DESeq2)
library(ggfortify)
library(glmnet) # package for regularization in GLMs
library(repr)


## input data
dir_nam = '../data/'
files_nam = c('TCGA-UCEC_expression_data.csv', 'TCGA-UCEC_clinical_data.csv')
res_dir = '../outputs/'
genemapping_file <- paste(dir_nam, 'TCGA-UCEC_geneSymbol_ensembleID_mapping.csv', sep='')
hmr2_genes_mapping_file <- paste(dir_nam, 'HMRdatabase2_00_GENES.csv', sep='')

## output directory
output_dir = '../outputs/'
vst_hmr2_filtered_data_file <- paste(output_dir, 'vst_hmr2_filtered_data.csv', sep='')

## Load data
read_exp_file <- function(file){
    exp_data <- read.csv(file)
    row.names(exp_data) <- exp_data$barcode
    exp_data$barcode <- NULL
    exp_data <- t(exp_data)
    print(dim(exp_data))
    return(exp_data)
}

read_clinical_file <- function(file){
    clinical_data <- read.csv(file)
    row.names(clinical_data) <- clinical_data$barcode
    clinical_data$barcode <- NULL
    print(dim(clinical_data))
    return(clinical_data)    
}

ucec_exp_data <- read_exp_file(paste(dir_nam, files_nam[1], sep=''))
#head(ucec_exp_data, 3)

ucec_clinical_data <- read_clinical_file(paste(dir_nam, files_nam[2], sep=''))
#head(ucec_clinical_data, 1)

gene_mapping <- read.csv(genemapping_file)
rownames(gene_mapping) <- gene_mapping$X
gene_mapping$X <- NULL
#print(dim(gene_mapping))
#head(gene_mapping, 2)

hmr2_genes_mapping =read.csv(hmr2_genes_mapping_file)
rownames(hmr2_genes_mapping) <- hmr2_genes_mapping$GENE.NAME
hmr2_genes_mapping$X. <- NULL
#print(dim(hmr2_genes_mapping))
#head(hmr2_genes_mapping, 2)

ucec_prim_exp_data <- ucec_exp_data[, rownames(ucec_clinical_data[ucec_clinical_data['sample_type'] == 'Primary Tumor',])]
ucec_prim_clinical_data <-ucec_clinical_data[rownames(ucec_clinical_data[ucec_clinical_data['sample_type']=='Primary Tumor',]),]
#print(dim(ucec_prim_exp_data))
#print(dim(ucec_prim_clinical_data))

vst_norm <- function(count_data){
    ## Pre-filtering the dataset 
    keep = rowSums(count_data >= 10) >= 5
    counts_keep <- count_data[keep,]
    print(nrow(counts_keep))

    ## Variance stabilizing transformation [Sample- level Quality Control]
    vst_prim <- vst(counts_keep, blind = FALSE)
    print(dim(vst_prim))
    return(vst_prim)
}

# Normalization
ucec_vst_data = vst_norm(ucec_prim_exp_data)
#head(ucec_vst_data,2)

# Genes filtering from HMR2
hmr2_genes = read.csv(hmr2_genes_mapping_file)
hmr2_genes <- hmr2_genes$GENE.NAME
ucec_vst_data = as.data.frame(ucec_vst_data[intersect(rownames(ucec_vst_data), hmr2_genes), ])
cat("VST data dimension:", dim(ucec_vst_data), sep=" ")

# Save the output
write.csv(ucec_vst_data, file=vst_hmr2_filtered_data_file, quote=F)
