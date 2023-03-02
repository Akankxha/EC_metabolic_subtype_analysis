# Import the libraries
library(dplyr)
library(NMF)
library(vcd)
library(survival)
library(ggstatsplot)
library(ggplot2)
library(DESeq2)
library(ggfortify)
library(ggthemes)

## Output Dir
output_dir = '../outputs/'

## input data
dir_nam = '../data/'
files_nam = c('TCGA-UCEC_expression_data.csv', 'TCGA-UCEC_clinical_data.csv')
res_dir = '../outputs/'
genemapping_file <- paste(dir_nam, 'TCGA-UCEC_geneSymbol_ensembleID_mapping.csv', sep='')
hmr2_genes_mapping <- paste(dir_nam, 'HMRdatabase2_00_GENES.csv', sep='')
vst_hmr2_filtered_data_file <- paste(output_dir, 'vst_hmr2_filtered_data.csv', sep='')
ucec_prim_clinical_data_file <- paste(output_dir, 'ucec_prim_clinical_data.csv', sep='')

# parameter
num_run_nmf = 50
num_run_nmf = 2


## Load Data
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

ucec_vst_data <- read.csv(vst_hmr2_filtered_data_file)
rownames(ucec_vst_data) <- ucec_vst_data$X
ucec_vst_data$X <- NULL

ucec_clinical_data <- read_clinical_file(paste(dir_nam, files_nam[2], sep=''))
ucec_exp_data <- read_exp_file(paste(dir_nam, files_nam[1], sep=''))

## MAD Score computations
compute_sort_mad <- function(vst_data, top = 1500){
    mad_df = as.data.frame(sapply(as.data.frame(t(vst_data)), mad))
    colnames(mad_df) <- 'median_abs_dev'
    mad_df = arrange(mad_df, median_abs_dev)
    top_mad_genes <- rownames(tail(mad_df, top))
    return(list(top_mad_genes, mad_df))  
}

ucec_mad <- compute_sort_mad(ucec_vst_data, top = 1000)
mad_df <-  ucec_mad[[2]]
ucec_top_mad_genes <- ucec_mad[[1]]
print(dim(mad_df))
print(summary(mad_df))
length(ucec_top_mad_genes)
#write.csv(mad_df, file=paste(res_dir,'genefeatures_mad_data.csv', sep=''),quote=F)
# write.csv(gene_mapping[ucec_top_mad_genes, ],file=paste(res_dir,'mad_genes.csv', sep=''),
#          quote=F,row.names = FALSE)
#print_genes(gene_mapping[ucec_top_mad_genes,'external_gene_name'])

hist(mad_df$median_abs_dev)

# list all available algorithms
print(nmfAlgorithm()[1:6])

nmf_run <- function(data, rank , method = 'brunet', seed_meth = "random", n_run = 10){
    res <- nmf(data, rank , method, seed= seed_meth, nrun = n_run, .options='t')
    return(res)
}

# Rank Estimation
ucec_rank_estimation_res = nmf_run(ucec_vst_data[ucec_top_mad_genes,], 2:7, seed_meth =123456, n_run = num_run_nmf)

ucec_rank_estimation_res_summary <- t(as.data.frame(summary(ucec_rank_estimation_res)))
write.csv(ucec_rank_estimation_res_summary,file=paste(res_dir,'rank_estimation_summary.csv', sep=''),quote=F)
#ucec_rank_estimation_res_summary


ucec_rank_summ_plot <-  plot(ucec_rank_estimation_res) + theme_base(base_size = 20) +
                             #theme_classic(base_size = 20) + 
                             theme(axis.text = element_text(face="bold"), plot.title = element_text(face="bold"),
                                  plot.caption = element_text(face = "bold")) + 
                                  labs(title = "(a). NMF Rank Estimation")
                    
#ucec_rank_summ_plot

ucec_rank_summ_plot_file <- file.path(paste(res_dir, 'rank_summ_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(ucec_rank_summ_plot_file, height = 5, width = 5, units = 'in', res=300)

# Print your heatmap
ucec_rank_summ_plot

# Close the PNG file:
dev.off()

consensusmap(ucec_rank_estimation_res)

ucec_rank_consensus_plot_file <- file.path(paste(res_dir, 'rank_consensus_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(ucec_rank_consensus_plot_file, height = 10, width = 10, units = 'in', res=300)

# Print your heatmap
consensusmap(ucec_rank_estimation_res)

# Close the PNG file:
dev.off()

# Experimentation with NMF
ucec_method_res = nmf_run(ucec_vst_data[ucec_top_mad_genes,], 2, method = nmfAlgorithm()[1:6], seed_meth =123456, n_run = num_run_nmf)
#t(compare(ucec_method_res))
write.csv(t(compare(ucec_method_res)),file=paste(res_dir, 'algo_compare_summary.csv', sep=''),quote=F)

plot(ucec_method_res)

ucec_residuals_plot_file <- file.path(paste(res_dir, 'algo_compare_residual_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(ucec_residuals_plot_file, height = 5, width = 5, units = 'in', res=300)

plot(ucec_method_res)

# Close the PNG file:
dev.off()

consensusmap(ucec_method_res)

ucec_consensus_plot_file <- file.path(paste(res_dir, 'algo_compare_consensus_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(ucec_consensus_plot_file, height = 9, width = 10, units = 'in', res=300)

# Print your heatmap
consensusmap(ucec_method_res)

# Close the PNG file:
dev.off()


# Run the nmf with optimal number of cluster and best method
ucec_best_res = nmf_run(ucec_vst_data[ucec_top_mad_genes,], 2, method = 'offset', seed_meth ='random', n_run = num_run_nmf)

ucec_best_res_summary <- t(as.data.frame(summary(ucec_best_res)))
#write.csv(t(ucec_best_res_summary),file=paste(res_dir,'best_res_summary.csv', sep=''),quote=F)
#ucec_best_res_summary
#compare(ucec_best_res)

plot(ucec_best_res)

ucec_prim_clinical_data <-ucec_clinical_data[rownames(ucec_clinical_data[ucec_clinical_data['sample_type']=='Primary Tumor',]),]
ucec_prim_clinical_data$nmf_cluster <- NA
ucec_prim_clinical_data$Subtypes <- NA

ucec_prim_exp_data <- ucec_exp_data[, rownames(ucec_clinical_data[ucec_clinical_data['sample_type'] == 'Primary Tumor',])]

ucec_cp <- consensusmap(ucec_best_res, , annCol= ucec_clinical_data[colnames(ucec_prim_exp_data),
            ][c('clinical_stage', 'histological_type', 'histological_grade')], 
            main='UCEC Consensus Matrix (Rank = 2)', tracks =c(), Rowv = FALSE)

ucec_prim_clinical_data[colnames(ucec_vst_data[lapply(cut(ucec_cp$Rowv,0.5)$lower,
                                    function(l)rapply(l,function(i) i))[[1]]]),'nmf_cluster']<-'Cluster_1'
ucec_prim_clinical_data[colnames(ucec_vst_data[lapply(cut(ucec_cp$Rowv,0.5)$lower,
                                    function(l)rapply(l,function(i) i))[[2]]]),'nmf_cluster'] <- 'Cluster_2'
                                                      
ucec_prim_clinical_data[colnames(ucec_vst_data[lapply(cut(ucec_cp$Rowv,0.5)$lower,
                                    function(l)rapply(l,function(i) i))[[1]]]),'Subtypes']<-'Metabolic_subtype-1'
ucec_prim_clinical_data[colnames(ucec_vst_data[lapply(cut(ucec_cp$Rowv,0.5)$lower,
                                    function(l)rapply(l,function(i) i))[[2]]]),'Subtypes'] <-'Metabolic_subtype-2'
                                                                                                    

# write.csv(ucec_prim_clinical_data,file=paste(res_dir,'prim_clinical_cluster_assign_data.csv', sep=''),quote=F)

print(table(ucec_prim_clinical_data['nmf_cluster']))                                                           

print(table(ucec_prim_clinical_data['Subtypes']))          

# Generate the consensus plot with clinical information                                                      
old_hist_type = c("Endometrioid endometrial adenocarcinoma", "Serous endometrial adenocarcinoma", 
                  "Mixed serous and endometrioid" )
new_hist_type = c("Endometrioid", "Serous", "Mixed")
hist_map = setNames(new_hist_type, old_hist_type)
ucec_clinical_data$histological_type <- hist_map[ucec_clinical_data$histological_type]
print(dim(ucec_clinical_data))

consensusmap(ucec_best_res,
        annCol=ucec_clinical_data[colnames(ucec_prim_exp_data),][c('clinical_stage', 'histological_type', 'histological_grade')],
        annRow=ucec_clinical_data[colnames(ucec_prim_exp_data),][c('clinical_stage', 'histological_type', 'histological_grade')],
        main='(b). Consensus Matrix', 
        tracks=c(),
        info = FALSE,
        Rowv = FALSE,
        labRow = NA,
        fontsize = 14
        )

ucec_best_consensus_plot_file <- file.path(paste(res_dir, 'best_consensus_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(ucec_best_consensus_plot_file, height = 6, width = 6, units = 'in', res= 300)

# Print your heatmap
consensusmap(ucec_best_res,
        annCol=ucec_clinical_data[colnames(ucec_prim_exp_data),][c('clinical_stage', 'histological_type', 'histological_grade')],
        annRow=ucec_clinical_data[colnames(ucec_prim_exp_data),][c('clinical_stage', 'histological_type', 'histological_grade')],
        main='(b). Consensus Matrix', 
        tracks=c(),
        info = FALSE,
        Rowv = FALSE,
        labRow = NA,
        fontsize = 14
        )

# Close the PNG file:
dev.off()

# Save the output
write.csv(ucec_prim_clinical_data, file=ucec_prim_clinical_data_file, quote=F)
