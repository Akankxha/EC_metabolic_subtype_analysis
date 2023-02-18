# Import Important Libraries
library("pheatmap")
library("RColorBrewer")
library('apeglm')
library("tidyverse")
library('ggrepel')
library('DEGreport')
library("EnhancedVolcano")

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

## Output Dir
output_dir = '../outputs/'

## input data
dir_nam = '../data/'
files_nam = c('_expression_data.csv', '_clinical_data.csv', '_vst_data.csv')
res_dir = '../outputs/'
genemapping_file <- paste(dir_nam, 'TCGA-UCEC_geneSymbol_ensembleID_mapping.csv')
hmr2_genes_mapping <- paste(dir_nam, 'HMRdatabase2_00.xlsx - GENES.csv')
hmr2_genes = hmr2_genes_mapping
vst_hmr2_filtered_data_file <- paste(output_dir, 'vst_hmr2_filtered_data.csv')

# load existing functions
source('../R-Scripts/data_preprocessing.r')


## Load Data
ucec_vst_data <- read.csv(vst_hmr2_filtered_data_file)
ucec_clinical_data <- read_clinical_file(paste(dir_nam, files_nam[2], sep=''))
ucec_exp_data <- read_exp_file(paste(dir_nam, files_nam[1], sep=''))

ucec_count_data <- ucec_prim_exp_data[rownames(ucec_vst_data), ]
print(dim(ucec_count_data))

# Check that sample names match in both count data and stage data
all(colnames(ucec_count_data) %in% rownames(ucec_prim_clinical_data['nmf_cluster']))
all(colnames(ucec_count_data) == rownames(ucec_prim_clinical_data['nmf_cluster']))

# construct the DESeqDataSet object from the matrix of counts and the sample information

ddsMat <- DESeqDataSetFromMatrix(countData = ucec_count_data,
                                 colData = ucec_prim_clinical_data['nmf_cluster'],
                                 design = ~nmf_cluster)
ddsMat

ddsMat <- DESeq(ddsMat)
ddsMat

# To know the name of coef in function lfcShrink
resultsNames(ddsMat)

## Extract results 
dds_res <- results(ddsMat, alpha = 0.05)
head(dds_res)

# Apply fold change shrinkage
dds_res <- lfcShrink(ddsMat, coef="nmf_cluster_Cluster_2_vs_Cluster_1", type="apeglm")
dds_res

dds_res %>% data.frame() %>% head()

summary(dds_res, alpha = 0.05)

summary(dds_res, alpha = 0.01)

# Create a tibble of results
dds_res_tb <- dds_res %>% data.frame()
print(dim(dds_res_tb))
head(dds_res_tb)

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- c(1, 0.5, 0.6)

# Subset the tibble to keep only significant genes
sig_res1 <- dds_res_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff[1])
print(dim(sig_res1))
#head(sig_res1, 4)

sig_res2 <- dds_res_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff[2])
print(dim(sig_res2))
#head(sig_res2, 4)

sig_res3 <- dds_res_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff[3])
print(dim(sig_res3))
#head(sig_res3, 4)

table(sig_res1['log2FoldChange'] > 0)

table(sig_res2['log2FoldChange'] > 0)

table(sig_res3['log2FoldChange'] > 0)

up_regulated_genes1 <- rownames(sig_res1[sig_res1['log2FoldChange']>0,])
print(length(up_regulated_genes1))

down_regulated_genes1 <- rownames(sig_res1[sig_res1['log2FoldChange']<0,])
print(length(down_regulated_genes1))

up_regulated_genes2 <- rownames(sig_res2[sig_res2['log2FoldChange']>0,])
print(length(up_regulated_genes2))

down_regulated_genes2 <- rownames(sig_res2[sig_res2['log2FoldChange']<0,])
print(length(down_regulated_genes2))

up_regulated_genes3 <- rownames(sig_res3[sig_res3['log2FoldChange']>0,])
print(length(up_regulated_genes3))

down_regulated_genes3 <- rownames(sig_res3[sig_res3['log2FoldChange']<0,])
print(length(down_regulated_genes3))

#print_genes(gene_mapping[up_regulated_genes3,'external_gene_name'])
#print_genes(gene_mapping[down_regulated_genes3,'external_gene_name'])

gene_mapping[rownames(head(sig_res1[order(-sig_res1[,2]),], 10)), "external_gene_name"]

head(sig_res1[order(-sig_res1[,2]),], 10)

gene_mapping[rownames(head(sig_res1[order(sig_res1[,2]),], 10)), "external_gene_name"]

head(sig_res1[order(sig_res1[,2]),], 10)

library(EnhancedVolcano)

sig_res1$external_gene_name <- gene_mapping[rownames(sig_res1),'external_gene_name']
sig_res1$entrez_id <- hmr2_genes_mapping[rownames(sig_res1),'GENE.ID.2']
sig_res1$hmr2_gene_names <- hmr2_genes_mapping[rownames(sig_res1),'SHORT.NAME']
#write.csv(sig_res1, file=paste(res_dir,'volcano_plots/clust2_DEGs_table.csv', sep=''),quote=F)

dds_res_tb <- dds_res_tb[order(dds_res_tb[,5], -dds_res_tb[,2]),]

dds_res_tb$external_gene_name <- gene_mapping[rownames(dds_res_tb),'external_gene_name']
dds_res_tb$entrez_id <- hmr2_genes_mapping[rownames(dds_res_tb),'GENE.ID.2']
dds_res_tb$hmr2_gene_names <- hmr2_genes_mapping[rownames(dds_res_tb),'SHORT.NAME']
#write.csv(dds_res_tb, file=paste(res_dir,'volcano_plots/clust2_deseq2_result.csv', sep=''),quote=F)
head(dds_res_tb)

#filter(dds_res_tb, log2FoldChange < -1.5)[order(filter(dds_res_tb, log2FoldChange < -1.5)$log2FoldChange), ]

options(repr.plot.width = 5, repr.plot.height = 5)
vp1 <- EnhancedVolcano(dds_res_tb,
                       lab = dds_res_tb$external_gene_name,
                       selectLab = c( "L1CAM","FOLR3", 'TKTL1', 'GCK', "LPCAT2", "CBS", "SULT1E1",
                                     "ELOVL2",  "GLDC", "ERBB2", "SRD5A3","PLPP2", "PLA2G2C", 
                                     "DPEP1"  ),
                       x = 'log2FoldChange',
                       y = 'padj', 
                       xlab = "log2 fold change",
                       ylab = "-log10 adjusted p-value",
                       ylim = c(0,160),
                       xlim = c(-6,6),
                       colAlpha = 1,
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       title = "Cluster2 / Cluster1 (logFC > 1)",
                       pointSize = 1,
                       labSize = 3,
                       drawConnectors = TRUE,widthConnectors = 0.5,
                       boxedLabels = TRUE
                       )

print(vp1)

#volcano_plot_file1 <- file.path(paste(res_dir, 'DEGs_results/clust21_logFC1.tiff', sep =''))
volcano_plot_file1 <- file.path(paste(res_dir, 'revision/clust21_logFC1.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(volcano_plot_file1, height = 6, width = 6, units = 'in', res= 300)

# Print your heatmap
print(vp1)

# Close the PNG file:
dev.off()

vp2 <- EnhancedVolcano(dds_res_tb,
                       lab = dds_res_tb$external_gene_name,
                       selectLab = c( "L1CAM","FOLR3", 'TKTL1', 'GCK', "LPCAT2", "ACSL5", "CA6", 'ZDHHC19', 
                                     "ELOVL2",  "GLDC", "ERBB2", "SRD5A3", "SMPD3","GSTA3", "ADH7", "SULT1E1"),
                       x = 'log2FoldChange',
                       y = 'padj', 
                       xlab = "log2 fold change",
                       ylab = "-log10 adjusted p-value",
                       ylim = c(0,160),
                       xlim = c(-6,6),
                       colAlpha = 1,
                       pCutoff = 0.05,
                       FCcutoff = 0.5,
                       title = "Cluster2 / Cluster1 (logFC > 0.5)",
                       pointSize = 1,
                       labSize = 3,
                       drawConnectors = TRUE,widthConnectors = 0.5,
                       boxedLabels = TRUE
                       )

print(vp2)

volcano_plot_file1b <- file.path(paste(res_dir, 'volcano_plots/clust21_logFC_0.5.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(volcano_plot_file1b, height = 6, width = 6, units = 'in', res=600)


# Print your heatmap
print(vp2)

# Close the PNG file:
dev.off()

vp3 <- EnhancedVolcano(dds_res_tb,
                       lab = dds_res_tb$external_gene_name,
                       selectLab = c( "L1CAM","FOLR3", 'TKTL1', 'GCK', "LPCAT2", "ACSL5", "", 'ZDHHC19', 
                                     "ELOVL2",  "GLDC", "ERBB2", "SRD5A3", "SMPD3","GSTA3", "", "SULT1E1"),
                       x = 'log2FoldChange',
                       y = 'padj', 
                       xlab = "log2 fold change",
                       ylab = "-log10 adjusted p-value",
                       ylim = c(0,160),
                       xlim = c(-6,6),
                       colAlpha = 1,
                       pCutoff = 0.05,
                       FCcutoff = 0.6,
                       title = "Cluster2 / Cluster1 (logFC > 0.6)",
                       pointSize = 1,
                       labSize = 3,
                       drawConnectors = TRUE,widthConnectors = 0.5,
                       boxedLabels = TRUE
                       )

print(vp3)

volcano_plot_file1c <- file.path(paste(res_dir, 'volcano_plots/clust21_logFC_0.6.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(volcano_plot_file1c, height = 6, width = 6, units = 'in', res=600)


# Print your heatmap
print(vp3)

# Close the PNG file:
dev.off()

dea_ucec_exp = ucec_exp_data[rownames(ucec_vst_data), ]
print(dim(dea_ucec_exp))
#head(dea_ucec_exp)

ucec_clinical_data$sample_type <- sub(" ", "_", ucec_clinical_data$sample_type)
ucec_clinical_data$sample_type <- sub(" ", "_", ucec_clinical_data$sample_type)
table(ucec_clinical_data["sample_type"])

# Check that sample names match in both count data and stage data
all(colnames(dea_ucec_exp) %in% rownames(ucec_clinical_data["sample_type"]))
all(colnames(dea_ucec_exp) == rownames(ucec_clinical_data["sample_type"]))

# construct the DESeqDataSet object from the matrix of counts and the sample information

tn_ddsMat <- DESeqDataSetFromMatrix(countData = dea_ucec_exp,
                                 colData = ucec_clinical_data["sample_type"],
                                 design = ~sample_type)
tn_ddsMat

tn_ddsMat$sample_type <- relevel(tn_ddsMat$sample_type, ref = "Solid_Tissue_Normal")
tn_ddsMat <- DESeq(tn_ddsMat)
tn_ddsMat

# To know the name of coef in function lfcShrink
resultsNames(tn_ddsMat)

## Extract results 
tn_dds_res <- results(tn_ddsMat, alpha = 0.05)
head(tn_dds_res)

# Apply fold change shrinkage
tn_dds_res <- lfcShrink(tn_ddsMat, coef="sample_type_Primary_Tumor_vs_Solid_Tissue_Normal", type="apeglm")
tn_dds_res

tn_dds_res %>% data.frame() %>% head()

summary(tn_dds_res, alpha = 0.05)

# Create a tibble of results
tn_dds_res_tb <- tn_dds_res %>% data.frame()
print(dim(tn_dds_res_tb))
head(tn_dds_res_tb)

# Subset the tibble to keep only significant genes
tn_sig_res1 <- tn_dds_res_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff[1])
print(dim(tn_sig_res1))
#head(tn_sig_res1, 4)

table(tn_sig_res1['log2FoldChange'] > 0)

tn_up_regulated_genes1 <- rownames(tn_sig_res1[tn_sig_res1['log2FoldChange']>0,])
print(length(tn_up_regulated_genes1))

tn_down_regulated_genes1 <- rownames(tn_sig_res1[tn_sig_res1['log2FoldChange']<0,])
print(length(tn_down_regulated_genes1))

tn_dds_res_tb <- tn_dds_res_tb[order(tn_dds_res_tb[,5], -tn_dds_res_tb[,2]),]

tn_vp1 <- EnhancedVolcano(tn_dds_res_tb,
                       lab = gene_mapping[rownames(tn_dds_res_tb), 'external_gene_name'],
                       selectLab = gene_mapping[rownames(tn_dds_res_tb)[1:20],'external_gene_name'],
                       x = 'log2FoldChange',
                       y = 'padj', 
                       xlab = "log2 fold change",
                       ylab = "-log10 adjusted p-value",
                       ylim = c(0,150),
                       xlim = c(-6,6),
                       colAlpha = 1,
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       title = "Tumor / Solid_tissue_normal (logFC > 1)",
                       pointSize = 1,
                       labSize = 3,
                       drawConnectors = TRUE,widthConnectors = 0.5,
                       boxedLabels = TRUE
                       )

print(tn_vp1)

#print_genes(gene_mapping[tn_down_regulated_genes1, "external_gene_name"])

table(ucec_prim_clinical_data[ucec_prim_clinical_data['nmf_cluster']== 'Cluster_1',]['detailed_sample_type'])

table(ucec_prim_clinical_data[ucec_prim_clinical_data['nmf_cluster']== 'Cluster_2',]['detailed_sample_type'])

clust1_nt_samples <- c()
clust1_index <- filter(ucec_prim_clinical_data, detailed_sample_type == "normal_matched_tumor", 
                        nmf_cluster == 'Cluster_1')["index"]

for (i in clust1_index$index){
    clust1_nt_samples <- c(clust1_nt_samples,rownames(ucec_clinical_data)[grep(i, rownames(ucec_clinical_data))])
}
print(length(clust1_nt_samples))

clust2_nt_samples <- c()
clust2_index <- filter(ucec_prim_clinical_data, detailed_sample_type == "normal_matched_tumor", 
                        nmf_cluster == 'Cluster_2')["index"]

for (i in clust2_index$index){
    clust2_nt_samples <- c(clust2_nt_samples,rownames(ucec_clinical_data)[grep(i, rownames(ucec_clinical_data))])
}
print(length(clust2_nt_samples))

length(c(clust1_nt_samples, clust2_nt_samples))

clust1_count_data <- ucec_exp_data[rownames(ucec_vst_data), clust1_nt_samples]
clust2_count_data <- ucec_exp_data[rownames(ucec_vst_data), clust2_nt_samples]

ucec_clinical_data$sample_type <- gsub(" ", "_", ucec_clinical_data$sample_type)
clust1_sample_type <- ucec_clinical_data[clust1_nt_samples, ]['sample_type']
clust2_sample_type <- ucec_clinical_data[clust2_nt_samples, ]['sample_type']
print(dim(clust1_count_data))
print(dim(clust2_count_data))
print(dim(clust1_sample_type))
print(dim(clust2_sample_type))

# Check that sample names match in both count data and stage data
all(colnames(clust1_count_data) %in% rownames(clust1_sample_type['sample_type']))
all(colnames(clust1_count_data) == rownames(clust1_sample_type['sample_type']))

# construct the DESeqDataSet object from the matrix of counts and the sample information
clust1_ddsMat <- DESeqDataSetFromMatrix(countData = clust1_count_data,
                                 colData = clust1_sample_type['sample_type'],
                                 design = ~sample_type)
clust1_ddsMat

clust1_ddsMat$sample_type <- relevel(clust1_ddsMat$sample_type, ref = "Solid_Tissue_Normal")
clust1_ddsMat <- DESeq(clust1_ddsMat)
clust1_ddsMat

# To know the name of coef in function lfcShrink
resultsNames(clust1_ddsMat)

## Extract results 
clust1_dds_res <- results(clust1_ddsMat, alpha = 0.05)
head(clust1_dds_res)

# Apply fold change shrinkage
clust1_dds_res <- lfcShrink(clust1_ddsMat, coef="sample_type_Primary_Tumor_vs_Solid_Tissue_Normal", type="apeglm")
clust1_dds_res

clust1_dds_res %>% data.frame() %>% head()

summary(clust1_dds_res, alpha = 0.05)

# Create a tibble of results
clust1_dds_res_tb <- clust1_dds_res %>% data.frame()
print(dim(clust1_dds_res_tb))
head(clust1_dds_res_tb)

# Subset the tibble to keep only significant genes
clust1_sig_res1 <- clust1_dds_res_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff[1])
print(dim(clust1_sig_res1))
#head(clust1_sig_res1, 4)

table(clust1_sig_res1['log2FoldChange'] > 0)

clust1_up_regulated_genes1 <- rownames(clust1_sig_res1[clust1_sig_res1['log2FoldChange']>0,])
print(length(clust1_up_regulated_genes1))

clust1_down_regulated_genes1 <- rownames(clust1_sig_res1[clust1_sig_res1['log2FoldChange']<0,])
print(length(clust1_down_regulated_genes1))

#print_genes(gene_mapping[clust1_up_regulated_genes1,'external_gene_name'])
#print_genes(gene_mapping[clust1_down_regulated_genes1,'external_gene_name'])

clust1_dds_res_tb <- clust1_dds_res_tb[order(clust1_dds_res_tb[,5], -clust1_dds_res_tb[,2]),]

clust1_vp1 <- EnhancedVolcano(clust1_dds_res_tb,
                       lab = gene_mapping[rownames(clust1_dds_res_tb), 'external_gene_name'],
                       selectLab = gene_mapping[rownames(clust1_dds_res_tb)[1:15],'external_gene_name'],
                       x = 'log2FoldChange',
                       y = 'padj', 
                       xlab = "log2 fold change",
                       ylab = "-log10 adjusted p-value",
                       ylim = c(0,60),
                       xlim = c(-9,9),
                       colAlpha = 1,
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       title = "Cluster1 [Tumor / Normal] (|logFC| > 1)",
                       pointSize = 1,
                       labSize = 3,
                       drawConnectors = TRUE,widthConnectors = 0.5,
                       boxedLabels = TRUE
                       )

print(clust1_vp1)

volcano_plot_file2 <- file.path(paste(res_dir, 'clust1_volcano_plot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(volcano_plot_file2, height = 6, width = 6, units = 'in', res=300)

# Print your heatmap
print(clust1_vp1)

# Close the PNG file:
dev.off()

# Check that sample names match in both count data and stage data
all(colnames(clust2_count_data) %in% rownames(clust2_sample_type['sample_type']))
all(colnames(clust2_count_data) == rownames(clust2_sample_type['sample_type']))

# construct the DESeqDataSet object from the matrix of counts and the sample information
clust2_ddsMat <- DESeqDataSetFromMatrix(countData = clust2_count_data,
                                 colData = clust2_sample_type['sample_type'],
                                 design = ~sample_type)
clust2_ddsMat

clust2_ddsMat$sample_type <- relevel(clust2_ddsMat$sample_type, ref = "Solid_Tissue_Normal")
clust2_ddsMat <- DESeq(clust2_ddsMat)
clust2_ddsMat

# To know the name of coef in function lfcShrink
resultsNames(clust2_ddsMat)

## Extract results 
clust2_dds_res <- results(clust2_ddsMat, alpha = 0.05)
head(clust2_dds_res)

# Apply fold change shrinkage
clust2_dds_res <- lfcShrink(clust2_ddsMat, coef="sample_type_Primary_Tumor_vs_Solid_Tissue_Normal", type="apeglm")
clust2_dds_res

clust2_dds_res %>% data.frame() %>% head()

summary(clust2_dds_res, alpha = 0.05)

# Create a tibble of results
clust2_dds_res_tb <- clust2_dds_res %>% data.frame()
print(dim(clust2_dds_res_tb))
head(clust2_dds_res_tb)

# Subset the tibble to keep only significant genes
clust2_sig_res1 <- clust2_dds_res_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff[1])
print(dim(clust2_sig_res1))
#head(clust2_sig_res1, 4)

print(table(clust2_sig_res1['log2FoldChange'] > 0))

clust2_up_regulated_genes1 <- rownames(clust2_sig_res1[clust2_sig_res1['log2FoldChange']>0,])
print(length(clust2_up_regulated_genes1))

clust2_down_regulated_genes1 <- rownames(clust2_sig_res1[clust2_sig_res1['log2FoldChange']<0,])
print(length(clust2_down_regulated_genes1))

#print_genes(gene_mapping[clust2_up_regulated_genes1,'external_gene_name'])
#print_genes(gene_mapping[clust2_down_regulated_genes1,'external_gene_name'])

clust2_dds_res_tb <- clust2_dds_res_tb[order(clust2_dds_res_tb[,5], -clust2_dds_res_tb[,2]),]

clust2_vp1 <- EnhancedVolcano(clust2_dds_res_tb,
                       lab = gene_mapping[rownames(clust2_dds_res_tb), 'external_gene_name'],
                       selectLab = gene_mapping[rownames(clust2_dds_res_tb)[1:15],'external_gene_name'],
                       x = 'log2FoldChange',
                       y = 'padj', 
                       xlab = "log2 fold change",
                       ylab = "-log10 adjusted p-value",
                       ylim = c(0,60),
                       xlim = c(-9,9),
                       colAlpha = 1,
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       title = "Cluster2 [Tumor / Normal] (|logFC| > 1)",
                       pointSize = 1,
                       labSize = 3,
                       drawConnectors = TRUE,widthConnectors = 0.5,
                       boxedLabels = TRUE
                       )

print(clust2_vp1)

volcano_plot_file3 <- file.path(paste(res_dir, 'clust2_volcano_plot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(volcano_plot_file3, height = 6, width = 6, units = 'in', res=300)


# Print your heatmap
print(clust2_vp1)

# Close the PNG file:
dev.off()

ucec_clinical_data$sample_type <- gsub(" ", "_", ucec_clinical_data$sample_type)
mtn_sample_type <- filter(ucec_clinical_data, detailed_sample_type== "normal_matched_tumor" | 
                          detailed_sample_type== "tumor_matched_normal_sample")["sample_type"]
mtn_count_data <- ucec_exp_data[rownames(ucec_vst_data), rownames(mtn_sample_type)]
print(dim(mtn_count_data))
print(dim(mtn_sample_type))
table(mtn_sample_type$sample_type)

# Check that sample names match in both count data and stage data
all(colnames(mtn_count_data) %in% rownames(mtn_sample_type["sample_type"]))
all(colnames(mtn_count_data) == rownames(mtn_sample_type["sample_type"]))

# construct the DESeqDataSet object from the matrix of counts and the sample information
mtn_ddsMat <- DESeqDataSetFromMatrix(countData = mtn_count_data,
                                 colData = mtn_sample_type['sample_type'],
                                 design = ~sample_type)
mtn_ddsMat

mtn_ddsMat$sample_type <- relevel(mtn_ddsMat$sample_type, ref = "Solid_Tissue_Normal")
mtn_ddsMat <- DESeq(mtn_ddsMat)
mtn_ddsMat

# To know the name of coef in function lfcShrink
resultsNames(mtn_ddsMat)

## Extract results 
mtn_dds_res <- results(mtn_ddsMat, alpha = 0.05)

# Apply fold change shrinkage
mtn_dds_res <- lfcShrink(mtn_ddsMat, coef="sample_type_Primary_Tumor_vs_Solid_Tissue_Normal", type="apeglm")

mtn_dds_res %>% data.frame() %>% head()

summary(mtn_dds_res, alpha = 0.05)

# Create a tibble of results
mtn_dds_res_tb <- mtn_dds_res %>% data.frame()
print(dim(mtn_dds_res_tb))
head(mtn_dds_res_tb, 3)

# Subset the tibble to keep only significant genes
mtn_sig_res1 <- mtn_dds_res_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff[1])
print(dim(mtn_sig_res1))
#head(mtn_sig_res1, 4)

print(table(mtn_sig_res1['log2FoldChange'] > 0))

mtn_up_regulated_genes1 <- rownames(mtn_sig_res1[mtn_sig_res1['log2FoldChange']>0,])
print(length(mtn_up_regulated_genes1))

mtn_down_regulated_genes1 <- rownames(mtn_sig_res1[mtn_sig_res1['log2FoldChange']<0,])
print(length(mtn_down_regulated_genes1))

gene_mapping[rownames(head(mtn_sig_res1[order(-mtn_sig_res1[2] ),], 10)), "external_gene_name"]

mtn_sig_res1$external_gene_name <- gene_mapping[rownames(mtn_sig_res1),'external_gene_name']
mtn_sig_res1$entrez_id <- hmr2_genes_mapping[rownames(mtn_sig_res1),'GENE.ID.2']
mtn_sig_res1$hmr2_gene_names <- hmr2_genes_mapping[rownames(mtn_sig_res1),'SHORT.NAME']
#write.csv(mtn_sig_res1, file=paste(res_dir,'DEGs_results/tumor_DEGs_table.csv', sep=''),quote=F)

mtn_dds_res_tb <- mtn_dds_res_tb[order(mtn_dds_res_tb[,5], -mtn_dds_res_tb[,2]),]
mtn_dds_res_tb$external_gene_name <- gene_mapping[rownames(mtn_dds_res_tb),'external_gene_name']
mtn_dds_res_tb$entrez_id <- hmr2_genes_mapping[rownames(mtn_dds_res_tb),'GENE.ID.2']
mtn_dds_res_tb$hmr2_gene_names <- hmr2_genes_mapping[rownames(mtn_dds_res_tb),'SHORT.NAME']
#write.csv(mtn_dds_res_tb, file=paste(res_dir,'DEGs_results/tumor_deseq2_result.csv', sep=''),quote=F)

#mtn_dds_res_tb[abs(mtn_dds_res_tb$log2FoldChange)>5,]

mtn_vp1 <- EnhancedVolcano(mtn_dds_res_tb,
                       lab = mtn_dds_res_tb$external_gene_name,
                       selectLab = mtn_dds_res_tb$external_gene_name[1:15],
                       x = 'log2FoldChange',
                       y = 'padj', 
                       xlab = "log2 fold change",
                       ylab = "-log10 adjusted p-value",
                       ylim = c(0, 80),
                       xlim = c(-9,9),
                       colAlpha = 1,
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       title = " Tumor / Normal (|logFC| > 1)",
                       pointSize = 1,
                       labSize = 3,
                       drawConnectors = TRUE,widthConnectors = 0.5,
                       boxedLabels = TRUE
                       )

print(mtn_vp1)

mtn_volcano_plot_file <- file.path(paste(res_dir, 'volcano_plots/mTN_logFC1.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(mtn_volcano_plot_file, height = 6, width = 6, units = 'in', res=600)


# Print your heatmap
print(mtn_vp1)

# Close the PNG file:
dev.off()

