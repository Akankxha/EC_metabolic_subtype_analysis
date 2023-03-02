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
#head(dds_res)

# Apply fold change shrinkage
dds_res <- lfcShrink(ddsMat, coef="nmf_cluster_Cluster_2_vs_Cluster_1", type="apeglm")
#dds_res

#dds_res %>% data.frame() %>% head()

summary(dds_res, alpha = 0.05)


# Create a tibble of results
dds_res_tb <- dds_res %>% data.frame()
print(dim(dds_res_tb))
#head(dds_res_tb)

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- c(1)

# Subset the tibble to keep only significant genes
sig_res1 <- dds_res_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff[1])
print(dim(sig_res1))
#head(sig_res1, 4)


table(sig_res1['log2FoldChange'] > 0)


up_regulated_genes1 <- rownames(sig_res1[sig_res1['log2FoldChange']>0,])
print(length(up_regulated_genes1))

down_regulated_genes1 <- rownames(sig_res1[sig_res1['log2FoldChange']<0,])
print(length(down_regulated_genes1))

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
volcano_plot_file1 <- file.path(paste(res_dir, 'clust21_logFC1.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(volcano_plot_file1, height = 6, width = 6, units = 'in', res= 300)

# Print your heatmap
print(vp1)

# Close the PNG file:
dev.off()


##############################################################

# DEA - normal vs tumor
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

