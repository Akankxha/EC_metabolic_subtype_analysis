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

## PCA 
pca_df = merge(t(ucec_exp_data[rownames(ucec_vst_data),]) , ucec_clinical_data["sample_type"], by="row.names") 
row.names(pca_df) = pca_df$Row.names
pca_df = pca_df[, -1]
print(dim(pca_df))
#pca_df
pca_res <- prcomp( pca_df[, - 3585], center = TRUE, scale. = TRUE)

## plot
autoplot(pca_res, data = pca_df, colour = "sample_type", frame = TRUE)

pca_plot_file <- file.path(paste(res_dir, 'pca_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(pca_plot_file, height = 5, width = 6, units = 'in', res=300)

# Print your heatmap
autoplot(pca_res, data = pca_df, colour = 'sample_type')

# Close the PNG file:
dev.off()

pca_plot_file2 <- file.path(paste(res_dir, 'pca_plot_boundary.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(pca_plot_file2, height = 5, width = 6, units = 'in', res=300)

# Print your heatmap
autoplot(pca_res, data = pca_df, colour = 'sample_type', frame = TRUE)

# Close the PNG file:
dev.off()

autoplot(pca_res, data = pca_df, colour = 'sample_type', frame = TRUE)

pca_df2 = merge(t(ucec_exp_data[rownames(ucec_vst_data),]) , ucec_clinical_data["detailed_sample_type"], 
                by="row.names") 
row.names(pca_df2) = pca_df2$Row.names
pca_df2 = pca_df2[, -1]
print(dim(pca_df2))
#pca_df
pca_res2 <- prcomp( pca_df2[, - 3585], center = TRUE, scale. = TRUE)

pca_plot_file3 <- file.path(paste(res_dir, 'pca_3grp_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(pca_plot_file3, height = 5, width = 5, units = 'in', res=600)

# Print your heatmap
autoplot(pca_res2, data = pca_df2, colour = 'detailed_sample_type')

# Close the PNG file:
dev.off()

pca_plot_file4 <- file.path(paste(res_dir, 'pca_3grp_plot_boundary.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(pca_plot_file4, height = 5, width = 5, units = 'in', res=600)

# Print your heatmap
autoplot(pca_res2, data = pca_df2, colour = 'detailed_sample_type', frame = TRUE)

# Close the PNG file:
dev.off()

autoplot(pca_res2, data = pca_df2, colour = 'detailed_sample_type', frame = TRUE)