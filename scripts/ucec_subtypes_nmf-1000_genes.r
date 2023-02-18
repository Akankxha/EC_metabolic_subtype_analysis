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

dir_nam = '/home/akansha/downloads/Pan_Cancer_dataset/gene_exp/prepared_datasets/TCGA-UCEC'
files_nam = c('_expression_data.csv', '_clinical_data.csv', '_vst_data.csv')
res_dir = '/home/akansha/downloads/metabolic_characterization/met_subtypes_nmf/ucec_nmf_results_542sample_1000genes/'

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
head(ucec_exp_data, 3)

ucec_clinical_data <- read_clinical_file(paste(dir_nam, files_nam[2], sep=''))
head(ucec_clinical_data, 1)

genemapping_file <- '/home/akansha/downloads/Pan_Cancer_dataset/gene_exp/TCGA-UCEC_geneSymbol_ensembleID_mapping.csv'
gene_mapping <- read.csv(genemapping_file)
rownames(gene_mapping) <- gene_mapping$X
gene_mapping$X <- NULL
print(dim(gene_mapping))
head(gene_mapping, 2)

hmr2_genes_mapping =read.csv('/home/akansha/downloads/metabolic_characterization/HMRdatabase2_00.xlsx - GENES.csv')
rownames(hmr2_genes_mapping) <- hmr2_genes_mapping$GENE.NAME
hmr2_genes_mapping$X. <- NULL
print(dim(hmr2_genes_mapping))
head(hmr2_genes_mapping, 2)

ucec_prim_exp_data <- ucec_exp_data[, rownames(ucec_clinical_data[ucec_clinical_data['sample_type'] == 'Primary Tumor',])]
ucec_prim_clinical_data <-ucec_clinical_data[rownames(ucec_clinical_data[ucec_clinical_data['sample_type']=='Primary Tumor',]),]
print(dim(ucec_prim_exp_data))
print(dim(ucec_prim_clinical_data))

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

ucec_vst_data = vst_norm(ucec_prim_exp_data)
head(ucec_vst_data,2)

hmr2_gene_file = '/home/akansha/downloads/metabolic_characterization/HMRdatabase2_00.xlsx - GENES.csv'
hmr2_genes = read.csv(hmr2_gene_file)
hmr2_genes <- hmr2_genes$GENE.NAME
length(hmr2_genes)

length(intersect(rownames(ucec_prim_exp_data), hmr2_genes))

ucec_vst_data = as.data.frame(ucec_vst_data[intersect(rownames(ucec_vst_data), hmr2_genes), ])
print(dim(ucec_vst_data))

pca_df = merge(t(ucec_exp_data[rownames(ucec_vst_data),]) , ucec_clinical_data["sample_type"], by="row.names") 
row.names(pca_df) = pca_df$Row.names
pca_df = pca_df[, -1]
print(dim(pca_df))
#pca_df
pca_res <- prcomp( pca_df[, - 3585], center = TRUE, scale. = TRUE)

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

ucec_rank_estimation_res = nmf_run(ucec_vst_data[ucec_top_mad_genes,], 2:7, seed_meth =123456, n_run = 50)

ucec_rank_estimation_res_summary <- t(as.data.frame(summary(ucec_rank_estimation_res)))
write.csv(ucec_rank_estimation_res_summary,file=paste(res_dir,'rank_estimation_summary.csv', sep=''),quote=F)
ucec_rank_estimation_res_summary

library(ggthemes)

options(repr.plot.width = 10, repr.plot.height = 8)
ucec_rank_summ_plot <-  plot(ucec_rank_estimation_res) + theme_base(base_size = 20) +
                             #theme_classic(base_size = 20) + 
                             theme(axis.text = element_text(face="bold"), plot.title = element_text(face="bold"),
                                  plot.caption = element_text(face = "bold")) + 
                                  labs(title = "(a). NMF Rank Estimation")
                    
ucec_rank_summ_plot

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

ucec_method_res = nmf_run(ucec_vst_data[ucec_top_mad_genes,], 2, method = nmfAlgorithm()[1:6], seed_meth =123456, n_run = 50)

t(compare(ucec_method_res))
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

ucec_best_res = nmf_run(ucec_vst_data[ucec_top_mad_genes,], 2, method = 'offset', seed_meth ='random', n_run = 50)

ucec_best_res_summary <- t(as.data.frame(summary(ucec_best_res)))
#write.csv(t(ucec_best_res_summary),file=paste(res_dir,'best_res_summary.csv', sep=''),quote=F)
ucec_best_res_summary
#compare(ucec_best_res)

plot(ucec_best_res)

ucec_prim_clinical_data$nmf_cluster <- NA
ucec_prim_clinical_data$Subtypes <- NA
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

m = matrix(rnorm(100), nrow = 10)
ha = rowAnnotation(foo = 1:10, show_legend = FALSE)
ht = Heatmap(m) + ha
lgd = color_mapping_legend(ha@anno_list[[1]]@color_mapping)
# you can also construct a Legend object by
lgd = Legend(title = "foo", 
    col_fun = circlize::colorRamp2(c(0, 10), c("white", "red")),
    at = seq(0, 10, by = 2))
draw(ht, annotation_legend_list = lgd, annotation_legend_side = "left")

old_hist_type = c("Endometrioid endometrial adenocarcinoma", "Serous endometrial adenocarcinoma", 
                  "Mixed serous and endometrioid" )
new_hist_type = c("Endometrioid", "Serous", "Mixed")
hist_map = setNames(new_hist_type, old_hist_type)
ucec_clinical_data$histological_type <- hist_map[ucec_clinical_data$histological_type]
print(dim(ucec_clinical_data))

options(repr.plot.width = 7, repr.plot.height = 7)
consensusmap(ucec_best_res,
        annCol=ucec_clinical_data[colnames(ucec_prim_exp_data),][c('clinical_stage', 'histological_type', 'histological_grade')],
        annRow=ucec_clinical_data[colnames(ucec_prim_exp_data),][c('clinical_stage', 'histological_type', 'histological_grade')],
        main='(b). Consensus Matrix', 
        tracks=c(),
        info = FALSE,
        Rowv = FALSE,
        labRow = NA,
        color = "Blues",
        fontsize = 14
        )

ucec_best_consensus_plot_file <- file.path(paste(res_dir, 'best_consensus_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(ucec_best_consensus_plot_file, height = 6, width = 6, units = 'in', res= 300)

# Print your heatmap
consensusmap(ucec_best_res, annCol= ucec_clinical_data[colnames(ucec_prim_exp_data),][c('clinical_stage',
            'histological_type', 'histological_grade')], main='UCEC Consensus Plot', tracks =c(),
            Rowv = FALSE)

# Close the PNG file:
dev.off()

# basis components
basismap(ucec_best_res, subsetRow=TRUE)
# mixture coefficients
coefmap(ucec_best_res)

# basis components
basismap(ucec_best_res)

ucec_best_basis_plot_file <- file.path(paste(res_dir, 'best_basis_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(ucec_best_basis_plot_file, height = 6, width = 6, units = 'in', res=300)

# Print your heatmap
basismap(ucec_best_res, subsetRow=TRUE)

# Close the PNG file:
dev.off()

ucec_best_coef_plot_file <- file.path(paste(res_dir, 'best_coef_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(ucec_best_coef_plot_file, height = 6, width = 6, units = 'in', res=300)

# Print your heatmap
coefmap(ucec_best_res)

# Close the PNG file:
dev.off()

ucec_best_full_basis_plot_file <- file.path(paste(res_dir, 'best_full_basis_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(ucec_best_full_basis_plot_file, height = 6, width = 6, units = 'in', res=300)

# Print your heatmap
basismap(ucec_best_res)

# Close the PNG file:
dev.off()

# Fitting the survival model
ucec_surv_fun = survfit(Surv(ucec_prim_clinical_data$OS.time/365,
                             ucec_prim_clinical_data$OS)~ ucec_prim_clinical_data$nmf_cluster)
ucec_surv_fun

#summary(ucec_surv_fun)
survdiff(Surv(ucec_prim_clinical_data$OS.time/365,ucec_prim_clinical_data$OS)~ ucec_prim_clinical_data$nmf_cluster)

options(repr.plot.width = 8, repr.plot.height = 6)
par (cex.lab = "1.2", cex.main = "1.5")
plot(ucec_surv_fun, 
     main = "(d). KM Plot for Metabolic Subtypes", 
     xlim = c(1,15), 
     xlab = substitute(paste(bold('Time (in years)'))) , 
     ylab = substitute(paste(bold('Survival Probability'))), 
     col = c("#F8766D", "#00BFC4"), 
     las= 1, 
     lwd = 2,  
     mark.time = TRUE, 
     cex=1.2)
legend(1, 0.2, legend = c('Cluster_1 (n=309)', 'Cluster_2 (n=233)'),lwd = 2, col = c("#F8766D", "#00BFC4"), bty = '', cex = 1.2)
text(10, 0.9,  'p-value = 3e-06', cex=1.3) 

#ucec_survival_plot_file <- file.path(paste(res_dir, 'survival_plot.tiff', sep =''))
ucec_survival_plot_file <- file.path(paste(res_dir, 'revision/survival_plot.tiff', sep =''))


# Open a PNG file - width and height arguments control the size of the output
tiff(ucec_survival_plot_file, height = 6, width = 6, units = 'in', res= 300)

# Print your heatmap
plot(ucec_surv_fun, 
     main = "KM Plot for UCEC Metabolic Subtypes", 
     xlim = c(1,15), 
     xlab = substitute(paste(bold('Time (in years)'))) , 
     ylab = substitute(paste(bold('Survival Probability'))), 
     col = c("#F8766D", "#00BFC4"), 
     las= 1, 
     lwd = 2,  
     mark.time = TRUE, 
     cex=1.2)
legend(1, 0.2, legend = c('Cluster_1 (n=309)', 'Cluster_2 (n=233)'),lwd = 2, col = c("#F8766D", "#00BFC4"), bty = '', cex = 1.2)
text(10, 0.9,  'p-value = 3e-06', cex=1.3) 

# Close the PNG file:
dev.off()

#hist(ucec_prim_clinical_data$OS.time)
#table(ucec_prim_clinical_data$OS.time)
#hist(ucec_prim_clinical_data$age_at_initial_pathologic_diagnosis)




#stage_ct <- as.matrix(table(ucec_prim_clinical_data$clinical_stage, ucec_prim_clinical_data$nmf_cluster))
stage_ct <- as.matrix(table(ucec_prim_clinical_data$clinical_stage, ucec_prim_clinical_data$Subtypes))
stage_ct[,1] <- (stage_ct[,1]*100)/sum(stage_ct[,1])
stage_ct[,2] <- (stage_ct[,2]*100)/sum(stage_ct[,2])
stage_ct

barplot(t(stage_ct), beside = TRUE, xlab = "Clinical Stage", ylab = "Percentage of Samples (%)", 
        main = "Association between Clinical stage and Clusters", col = c("#8470FF", "#F08080")) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#8470FF", "#F08080")) 

stage_barplot_file <- file.path(paste(res_dir, 'barplots/stage_barplot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(stage_barplot_file, height = 6, width = 6, units = 'in', res=600)


# Print your heatmap
barplot(t(stage_ct), beside = TRUE, xlab = "Clinical Stage", ylab = "Percentage of Samples (%)", 
        main = "Association between Clinical stage and Clusters", col = c("#8470FF", "#F08080")) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#8470FF", "#F08080")) 

# Close the PNG file:
dev.off()

ucec_stage_f_test <-fisher.test(table(ucec_prim_clinical_data$clinical_stage, ucec_prim_clinical_data$nmf_cluster))
print(ucec_stage_f_test)

stage_cramer <- assocstats(table(ucec_prim_clinical_data$clinical_stage, ucec_prim_clinical_data$nmf_cluster))
stage_cramer

ucec_stage_test <- chisq.test(table(ucec_prim_clinical_data$clinical_stage, ucec_prim_clinical_data$nmf_cluster))
ucec_stage_test

# plot
#ggbarstats(data = ucec_prim_clinical_data,
#           x = clinical_stage,
#           y = nmf_cluster, 
#           title = 'Pearsons chi-square test between NMF Clusters & Clinical Stage\n',
#           cex = 2.5) + labs(caption = NULL) 

# Save the plot
#ucec_chisq_stage_plot_file <- file.path(paste(res_dir, 'chi_sq_test_stage_plot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
#tiff(ucec_chisq_stage_plot_file, height = 6, width = 6, units = 'in', res=300)


# Print your heatmap
#ggbarstats(data = ucec_prim_clinical_data,
#           x = clinical_stage,
#           y = nmf_cluster, 
#           title = 'Pearsons chi-square test between NMF Clusters & Clinical Stage\n',
#           cex = 2.5, 
#          ) + labs(caption = NULL) # remove caption

# Close the PNG file:
#dev.off()

# Other functions for pearson chi-square test
# 1. summary(table(ucec_prim_clinical_data$clinical_stage, ucec_prim_clinical_data$nmf_cluster))
# 2. assocstats(table(ucec_prim_clinical_data$clinical_stage, ucec_prim_clinical_data$nmf_cluster)) - gives cramer's V

# Mosaic Plot
#mosaic(~ clinical_stage + nmf_cluster,
#  direction = c("v", "h"),
#  data = ucec_prim_clinical_data,
#  shade = TRUE
#)

#type_ct <- as.matrix(table(ucec_prim_clinical_data$histological_type, 
#                           ucec_prim_clinical_data$nmf_cluster))
type_ct <- as.matrix(table(ucec_prim_clinical_data$histological_type, 
                           ucec_prim_clinical_data$Subtypes))
type_ct[,1] <- (type_ct[,1]*100)/sum(type_ct[,1])
type_ct[,2] <- (type_ct[,2]*100)/sum(type_ct[,2])
rownames(type_ct) <- c("Endometrioid", "Mixed", "Serous")
type_ct

barplot(t(type_ct), beside = TRUE, xlab = "Histological Type", 
        ylab = "Percentage of Samples (%)", 
        main = "Association between Histological Type and Clusters", 
        col = c("#8470FF", "#F08080")) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#8470FF", "#F08080")) 

type_barplot_file <- file.path(paste(res_dir, 'barplots/type_barplot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(type_barplot_file, height = 6, width = 6, units = 'in', res=600)

# Print your heatmap
barplot(t(type_ct), beside = TRUE, xlab = "Histological Type", ylab = "Percentage of Samples (%)", 
        main = "Association between Histological Type and Clusters", col = c("#8470FF", "#F08080")) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#8470FF", "#F08080"))  

# Close the PNG file:
dev.off()

ucec_hist_type_f_test <- fisher.test(table(ucec_prim_clinical_data$histological_type, 
                                           ucec_prim_clinical_data$nmf_cluster))
print(ucec_hist_type_f_test)

assocstats(table(ucec_prim_clinical_data$histological_type, ucec_prim_clinical_data$nmf_cluster))

ucec_hist_type_test <- chisq.test(table(ucec_prim_clinical_data$histological_type, 
                                        ucec_prim_clinical_data$nmf_cluster))
ucec_hist_type_test

# plot
#ggbarstats(data = ucec_prim_clinical_data,
#           x = histological_type,
#           y = nmf_cluster, 
#           title = 'Pearsons chi-square test between NMF Clusters & Histological Type\n',
#           cex = 2.5, 
#          ) + labs(caption = NULL) # remove caption

#ucec_chisq_hist_type_plot_file <- file.path(paste(res_dir, 'chi_sq_test_hist_type_plot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
#tiff(ucec_chisq_hist_type_plot_file, height = 6, width = 6, units = 'in', res=300)


# Print your heatmap
#ggbarstats(data = ucec_prim_clinical_data,
#           x = histological_type,
#           y = nmf_cluster, 
#           title = 'Pearsons chi-square test between NMF Clusters & Histological Type\n',
#           cex = 2.5, 
#          ) + labs(caption = NULL) # remove caption

# Close the PNG file:
#dev.off()

#grade_ct <- as.matrix(table(ucec_prim_clinical_data$histological_grade, ucec_prim_clinical_data$nmf_cluster))
grade_ct <- as.matrix(table(ucec_prim_clinical_data$histological_grade, ucec_prim_clinical_data$Subtypes))
grade_ct[,1] <- (grade_ct[,1]*100)/sum(grade_ct[,1])
grade_ct[,2] <- (grade_ct[,2]*100)/sum(grade_ct[,2])
grade_ct

barplot(t(grade_ct), beside = TRUE, xlab = "Histological Grade", ylab = "Percentage of Samples (%)", 
        main = "Association between Histological Grade and Clusters", col = c("#8470FF", "#F08080")) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#8470FF", "#F08080")) 

grade_barplot_file <- file.path(paste(res_dir, 'barplots/grade_barplot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(grade_barplot_file, height = 6, width = 6, units = 'in', res=600)

# Print your heatmap
barplot(t(grade_ct), beside = TRUE, xlab = "Histological Grade", ylab = "Percentage of Samples (%)", 
        main = "Association between Histological Grade and Clusters", col = c("#8470FF", "#F08080")) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#8470FF", "#F08080")) 

# Close the PNG file:
dev.off()

ucec_grade_f_test <- fisher.test(table(ucec_prim_clinical_data$histological_grade, 
                                       ucec_prim_clinical_data$nmf_cluster))
print(ucec_grade_f_test)

assocstats(table(ucec_prim_clinical_data$histological_grade, ucec_prim_clinical_data$nmf_cluster))

ucec_grade_test <- chisq.test(table(ucec_prim_clinical_data$histological_grade, ucec_prim_clinical_data$nmf_cluster))
ucec_grade_test

# combine plot and statistical test with ggbarstats
#ucec_cv <- assocstats(table(ucec_prim_clinical_data$histological_grade, ucec_prim_clinical_data$nmf_cluster))
#ggbarstats(
#    ucec_prim_clinical_data, histological_grade, nmf_cluster,
#    results.subtitle = FALSE,
#    title = "Fisher's exact test\n",
#    subtitle = paste0("p-value = ", ucec_grade_f_test$p.value, "\n\nCramer's V = ", ucec_cv$cramer),
#    cex = 2.5)

#ucec_chisq_grade_plot_file <- file.path(paste(res_dir, 'chi_sq_test_grade_plot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
#png(ucec_chisq_grade_plot_file, width = 600, height = 600)
#tiff(ucec_chisq_grade_plot_file, height = 6, width = 6, units = 'in', res=300)


# Print your heatmap
#ggbarstats(
#    ucec_prim_clinical_data, histological_grade, nmf_cluster,
#    results.subtitle = FALSE,
#   title = "Fisher's exact test\n",
#    subtitle = paste0("p-value = ", ucec_grade_f_test$p.value, "\n\nCramer's V = ", ucec_cv$cramer),
#    cex = 2.5)

# Close the PNG file:
#dev.off()

ucec_prim_clinical_data$age <- ifelse(ucec_prim_clinical_data$age_at_initial_pathologic_diagnosis < 50,
  "<50", ">=50"
)
print(table(ucec_prim_clinical_data$age))
#print(table(ucec_prim_clinical_data$age, ucec_prim_clinical_data$nmf_cluster))
#age_ct <- as.matrix(table(ucec_prim_clinical_data$age, ucec_prim_clinical_data$nmf_cluster))
age_ct <- as.matrix(table(ucec_prim_clinical_data$age, ucec_prim_clinical_data$Subtypes))
age_ct[,1] <- (age_ct[,1]*100)/sum(age_ct[,1])
age_ct[,2] <- (age_ct[,2]*100)/sum(age_ct[,2])
age_ct

barplot(t(age_ct), beside = TRUE, xlab = "Age", ylab = "Percentage of Samples (%)", 
        main = "Association between Age and Clusters", col = c("#8470FF", "#F08080"),
        width = c(0.08), space = c(0.01, 1), xlim = c(0, 0.7)) 
legend("topleft", c("Cluster-1", "Cluster-2"), fill = c("#8470FF", "#F08080")) 

age_barplot_file <- file.path(paste(res_dir, 'barplots/age_barplot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(age_barplot_file, height = 6, width = 6, units = 'in', res=600)

# Print your heatmap
barplot(t(age_ct), beside = TRUE, xlab = "Age", ylab = "Percentage of Samples (%)", 
        main = "Association between Age and Clusters", col = c("#8470FF", "#F08080"),
        width = c(0.08), space = c(0.01, 1), xlim = c(0, 0.7)) 
legend("topleft", c("Cluster-1", "Cluster-2"), fill = c("#8470FF", "#F08080")) 

# Close the PNG file:
dev.off()

# (bottom, left, top, and right.)
# par(mar = c(4, 4, 0.1, 0.5))  


stage_ct_df <- as.data.frame(stage_ct)
colnames(stage_ct_df) <- c("clinical_stages", "Clusters", "Percentage")
stage_bar <- ggbarplot(stage_ct_df, x = "clinical_stages", y = "Percentage", fill = "Clusters", label = FALSE,  
          position = position_dodge(0.8),
          ylab = "Percentage of Samples (%)", xlab = "", title = "Clinical stage", 
          ggtheme = theme_classic(base_size = 18), legend = "top") +
          theme(plot.title = element_text(hjust = 0.5, face="bold"), 
                axis.text = element_text(face="bold"), 
                legend.title=element_blank())  +
          scale_y_continuous(breaks=seq(0, 100, 10)) 

type_ct_df <- as.data.frame(type_ct)
colnames(type_ct_df) <- c("histological_types", "Clusters", "Percentage")
type_bar  <-  ggbarplot(type_ct_df, x = "histological_types", y = "Percentage", fill = "Clusters", label = FALSE,  
          position = position_dodge(0.8),
          ylab = "Percentage of Samples (%)", xlab = "", title = "Histological Type", 
          ggtheme = theme_classic(base_size = 18), legend = "top") +
          theme(plot.title = element_text(hjust = 0.5, face="bold"), 
                axis.text = element_text(face="bold"), 
                legend.title=element_blank())  +
          scale_y_continuous(breaks=seq(0, 100, 10))

grade_ct_df <- as.data.frame(grade_ct)
colnames(grade_ct_df) <- c("histological_grades", "Clusters", "Percentage")
grade_bar  <-  ggbarplot(grade_ct_df, x = "histological_grades", y = "Percentage", fill = "Clusters", label = FALSE,  
          position = position_dodge(0.8),
          ylab = "Percentage of Samples (%)", xlab = "", title = "Histological Grade", 
          ggtheme = theme_classic(base_size = 18), legend = "top") +
          theme(plot.title = element_text(hjust = 0.5, face="bold"), 
                axis.text = element_text(face="bold"), 
                legend.title=element_blank())  +
          scale_y_continuous(breaks=seq(0, 100, 10))

age_ct_df <- as.data.frame(age_ct)
colnames(age_ct_df) <- c("Ages", "Clusters", "Percentage")
age_bar  <-  ggbarplot(age_ct_df, x = "Ages", y = "Percentage", fill = "Clusters", label = FALSE,  
          position = position_dodge(0.8),
          ylab = "Percentage of Samples (%)", xlab = "", title = "Age", 
          ggtheme = theme_classic(base_size = 18), legend = "top") +
          theme(plot.title = element_text(hjust = 0.5, face="bold"), 
                axis.text = element_text(face="bold"), 
                legend.title=element_blank())  +
          scale_y_continuous(breaks=seq(0, 100, 10))

options(repr.plot.width = 12, repr.plot.height = 10)
clinical_barplot <- ggarrange(stage_bar, type_bar,  grade_bar, age_bar,
          labels = c("a", "b", "c", "d"), font.label = list(size = 22),
          ncol = 2, nrow = 2) 
clinical_barplot

library(ggpubr)

stage_ct_df <- as.data.frame(stage_ct)
colnames(stage_ct_df) <- c("clinical_stages", "Clusters", "Percentage")
stage_bar <- ggbarplot(stage_ct_df, x = "clinical_stages", y = "Percentage", fill = "Clusters", label = FALSE,  
          position = position_dodge(0.8),
          ylab = "Percentage of Samples (%)", xlab = "", title = "Clinical stage", 
          ggtheme = theme_classic(base_size = 18), legend = "top") +
          theme(plot.title = element_text(hjust = 0.5, face="bold"), 
                axis.text = element_text(face="bold"), 
                legend.title=element_blank())  +
          scale_y_continuous(breaks=seq(0, 100, 10)) 

type_ct_df <- as.data.frame(type_ct)
colnames(type_ct_df) <- c("histological_types", "Clusters", "Percentage")
type_bar  <-  ggbarplot(type_ct_df, x = "histological_types", y = "Percentage", fill = "Clusters", label = FALSE,  
          position = position_dodge(0.8),
          ylab = "Percentage of Samples (%)", xlab = "", title = "Histological Type", 
          ggtheme = theme_classic(base_size = 18), legend = "top") +
          theme(plot.title = element_text(hjust = 0.5, face="bold"), 
                axis.text = element_text(face="bold"), 
                legend.title=element_blank())  +
          scale_y_continuous(breaks=seq(0, 100, 10))

grade_ct_df <- as.data.frame(grade_ct)
colnames(grade_ct_df) <- c("histological_grades", "Clusters", "Percentage")
grade_bar  <-  ggbarplot(grade_ct_df, x = "histological_grades", y = "Percentage", fill = "Clusters", label = FALSE,  
          position = position_dodge(0.8),
          ylab = "Percentage of Samples (%)", xlab = "", title = "Histological Grade", 
          ggtheme = theme_classic(base_size = 18), legend = "top") +
          theme(plot.title = element_text(hjust = 0.5, face="bold"), 
                axis.text = element_text(face="bold"), 
                legend.title=element_blank())  +
          scale_y_continuous(breaks=seq(0, 100, 10))

age_ct_df <- as.data.frame(age_ct)
colnames(age_ct_df) <- c("Ages", "Clusters", "Percentage")
age_bar  <-  ggbarplot(age_ct_df, x = "Ages", y = "Percentage", fill = "Clusters", label = FALSE,  
          position = position_dodge(0.8),
          ylab = "Percentage of Samples (%)", xlab = "", title = "Age", 
          ggtheme = theme_classic(base_size = 18), legend = "top") +
          theme(plot.title = element_text(hjust = 0.5, face="bold"), 
                axis.text = element_text(face="bold"), 
                legend.title=element_blank())  +
          scale_y_continuous(breaks=seq(0, 100, 10))

options(repr.plot.width = 12, repr.plot.height = 10)
clinical_barplot <- ggarrange(stage_bar, type_bar,  grade_bar, age_bar,
          labels = c("a.", "b.", "c.", "d."), font.label = list(size = 24),
          ncol = 2, nrow = 2) 
clinical_barplot

clinical_barplot_file <- file.path(paste(res_dir, 'revision/clinical_barplot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(clinical_barplot_file, height = 10, width = 12, units = 'in', res=300)

# Print your heatmap
clinical_barplot

# Close the PNG file:
dev.off()



options(repr.plot.width = 12, repr.plot.height = 8)
par(mar = c(3, 5, 1.5, 0.1)) 
par(mfrow=c(2,2), cex.main= 1.5 )
barplot(t(stage_ct), beside = TRUE, ylab = "Percentage of Samples (%)", 
        main = "Clinical stage", col = c("#F8766D", "#00BFC4")) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#F8766D", "#00BFC4")) 

barplot(t(type_ct), beside = TRUE, ylab = "Percentage of Samples (%)", 
        main = "Histological Type", col = c("#F8766D", "#00BFC4")) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#F8766D", "#00BFC4"))  


barplot(t(grade_ct), beside = TRUE, ylab = "Percentage of Samples (%)", 
        main = "Histological Grade", col = c("#F8766D", "#00BFC4")) 
legend("topleft", c("Cluster-1", "Cluster-2"), fill = c("#F8766D", "#00BFC4")) 


barplot(t(age_ct), beside = TRUE,  ylab = "Percentage of Samples (%)", 
        main = "Age", col = c("#F8766D", "#00BFC4"),
        width = c(0.08), space = c(0.01, 1), xlim = c(0, 0.7)) 
legend("topleft", c("Cluster-1", "Cluster-2"), fill = c("#F8766D", "#00BFC4")) 

combined_clinical_barplot_file <- file.path(paste(res_dir, 'barplots/combined_clinical_barplot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(combined_clinical_barplot_file, height = 6, width = 9, units = 'in', res= 300)

# Print your heatmap
par(mar = c(3, 5, 1.5, 0.1)) 
par(mfrow=c(2,2), cex.main= 1.5 )
barplot(t(stage_ct), beside = TRUE, ylab = "Percentage of Samples (%)", 
        main = "Clinical stage", col = c("#F8766D", "#00BFC4")) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#F8766D", "#00BFC4")) 


barplot(t(type_ct), beside = TRUE, ylab = "Percentage of Samples (%)", 
        main = "Histological Type", col = c("#F8766D", "#00BFC4")) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#F8766D", "#00BFC4"))  


barplot(t(grade_ct), beside = TRUE, ylab = "Percentage of Samples (%)", 
        main = "Histological Grade", col = c("#F8766D", "#00BFC4")) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#F8766D", "#00BFC4")) 

barplot(t(age_ct), beside = TRUE,  ylab = "Percentage of Samples (%)", 
        main = "Age", col = c("#F8766D", "#00BFC4"),
        width = c(0.08), space = c(0.01, 1), xlim = c(0, 0.7)) 
legend("topright", c("Cluster-1", "Cluster-2"), fill = c("#F8766D", "#00BFC4")) 

# Close the PNG file:
dev.off()

ucec_age_f_test <- fisher.test(table(ucec_prim_clinical_data$age, ucec_prim_clinical_data$nmf_cluster))
print(ucec_age_f_test)

assocstats(table(ucec_prim_clinical_data$age, ucec_prim_clinical_data$nmf_cluster))

ucec_age_test <- chisq.test(table(ucec_prim_clinical_data$age, ucec_prim_clinical_data$nmf_cluster))
ucec_age_test

# plot
#ggbarstats(data = ucec_prim_clinical_data,
#           x = age,
#           y = nmf_cluster, 
#           title = 'Pearsons chi-square test between NMF Clusters & Age\n',
#           cex = 2.5, 
#          ) + labs(caption = NULL) # remove caption

#ucec_chisq_age_plot_file <- file.path(paste(res_dir, 'chi_sq_test_age_plot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
#tiff(ucec_chisq_age_plot_file, height = 6, width = 6, units = 'in', res=300)


# Print your heatmap
#ggbarstats(data = ucec_prim_clinical_data,
#           x = age,
#           y = nmf_cluster, 
#           title = 'Pearsons chi-square test between NMF Clusters & Age\n',
#           cex = 2.5, 
#          ) + labs(caption = NULL) 

# Close the PNG file:
#dev.off()

W = ucec_best_res@fit@W
H = ucec_best_res@fit@H
print(dim(W))
print(dim(H))

# only compute the scores
ucec_best_res_feat_score <- featureScore(ucec_best_res)
print(summary(ucec_best_res_feat_score))

# compute the scores and characterize each metagene
ucec_best_res_feat <- extractFeatures(ucec_best_res)
print(str(ucec_best_res_feat))

clust1_feat <- W[ucec_best_res_feat[[1]],]
clust2_feat <- W[ucec_best_res_feat[[2]],]
print(dim(clust1_feat))
print(dim(clust2_feat))

clust1_genes <- rownames(clust1_feat)
clust2_genes <- rownames(clust2_feat)
print(length(clust1_genes))
print(length(clust2_genes))

print_genes <- function(lst_genes){
    for (i in lst_genes){
        cat(noquote(i), sep='\n')
    }
}

#print_genes(gene_mapping[clust1_genes,'external_gene_name'])
#print_genes(gene_mapping[clust2_genes,'external_gene_name'])

# Import Important Libraries
library("pheatmap")
library("RColorBrewer")
library('apeglm')
library("tidyverse")
library('ggrepel')
library('DEGreport')
library("EnhancedVolcano")

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

# print_genes(gene_mapping[mtn_down_regulated_genes1,'external_gene_name'])



summary(t(ucec_vst_data["ENSG00000130707", ]))

ass1 = as.data.frame(t(ucec_vst_data["ENSG00000130707", ]))
ass1$binary<- ifelse(ass1$ENSG00000130707> 11.5,"high","low")
head(ass1,3)

# Fitting the survival model
ass1_surv_fun = survfit(Surv(ucec_prim_clinical_data$OS.time/365,
                             ucec_prim_clinical_data$OS)~ ass1$binary)
ass1_surv_fun

survdiff(Surv(ucec_prim_clinical_data$OS.time/365,ucec_prim_clinical_data$OS)~ ass1$binary)

plot(ass1_surv_fun, main = "KM Plot for ASS1 expression", xlim = c(1,15), xlab = 'Time (in years)', ylab = 'Survival Probability', 
     col = c('red', 'blue'), las= 1, lwd = 2,  mark.time = TRUE, cex=1)
legend(1, 0.15, legend = c('high (n=272)', 'low (n=270)'),lwd = 2, col = c('red', 'blue'), bty = '', cex = 1)
text(10, 0.9,  'p-value = 2e-06', cex=1.3) 

summary(t(ucec_vst_data["ENSG00000141736", ]))  # ERBB2

erbb2 = as.data.frame(t(ucec_vst_data["ENSG00000141736", ]))
erbb2$binary<- ifelse(erbb2$ENSG00000141736>= 12.685,"high","low")
head(erbb2,3)

# Fitting the survival model
surv_fun_erbb2 = survfit(Surv(ucec_prim_clinical_data$OS.time/365,
                             ucec_prim_clinical_data$OS)~ erbb2$binary)
surv_fun_erbb2

survdiff(Surv(ucec_prim_clinical_data$OS.time/365,ucec_prim_clinical_data$OS)~ erbb2$binary)

plot(surv_fun_erbb2, main = "KM Plot for ERBB2 expression", xlim = c(1,15), xlab = 'Time (in years)', ylab = 'Survival Probability', 
     col = c('red', 'blue'), las= 1, lwd = 2,  mark.time = TRUE, cex=1)
legend(1, 0.15, legend = c('high (n=270)', 'low (n=272)'),lwd = 2, col = c('red', 'blue'), bty = '', cex = 1)
text(10, 0.9,  'p-value = 0.009', cex=1.3) 

summary(t(ucec_vst_data["ENSG00000135069", ]))  # PSAT1

psat1 = as.data.frame(t(ucec_vst_data["ENSG00000135069", ]))
psat1$binary<- ifelse(psat1$ENSG00000135069> 10.980,"high","low")
head(psat1,3)

# Fitting the survival model
surv_fun_psat1 = survfit(Surv(ucec_prim_clinical_data$OS.time/365,
                             ucec_prim_clinical_data$OS)~ psat1$binary)
surv_fun_psat1

survdiff(Surv(ucec_prim_clinical_data$OS.time/365,ucec_prim_clinical_data$OS)~ psat1$binary)

plot(surv_fun_psat1, main = "KM Plot for PSAT1 expression", xlim = c(1,15), xlab = 'Time (in years)', ylab = 'Survival Probability', 
     col = c('red', 'blue'), las= 1, lwd = 2,  mark.time = TRUE, cex=1)
legend(1, 0.15, legend = c('high (n=271)', 'low (n=271)'),lwd = 2, col = c('red', 'blue'), bty = '', cex = 1)
text(10, 0.9,  'p-value = 0.02', cex=1.3) 

summary(t(ucec_vst_data["ENSG00000106633", ]))  # GCK

gck = as.data.frame(t(ucec_vst_data["ENSG00000106633", ]))
gck$binary<- ifelse(gck$ENSG00000106633> 4.3,"high","low")
head(gck,3)

# Fitting the survival model
ucec_surv_fun_gck = survfit(Surv(ucec_prim_clinical_data$OS.time/365,
                             ucec_prim_clinical_data$OS)~ gck$binary)
ucec_surv_fun_gck

survdiff(Surv(ucec_prim_clinical_data$OS.time/365,ucec_prim_clinical_data$OS)~ gck$binary)

plot(ucec_surv_fun_gck, main = "KM Plot for GCK expression", xlim = c(1,15), xlab = 'Time (in years)', ylab = 'Survival Probability', 
     col = c('red', 'blue'), las= 1, lwd = 2,  mark.time = TRUE, cex=1)
legend(1, 0.15, legend = c('high (n=272)', 'low (n=270)'),lwd = 2, col = c('red', 'blue'), bty = '', cex = 1)
text(10, 0.9,  'p-value = 9e-04', cex=1.3) 

gck_survival_plot_file <- file.path(paste(res_dir, 'gck_survival_plot.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(gck_survival_plot_file, height = 6, width = 6, units = 'in', res=300)

# Print your heatmap
plot(ucec_surv_fun_gck, main = "KM Plot for GCK expression", xlim = c(1,15), xlab = 'Time (in years)', ylab = 'Survival Probability', 
     col = c('red', 'blue'), las= 1, lwd = 2,  mark.time = TRUE, cex=1)
legend(1, 0.15, legend = c('high (n=272)', 'low (n=270)'),lwd = 2, col = c('red', 'blue'), bty = '', cex = 1)
text(10, 0.9,  'p-value = 9e-04', cex=1.3) 

# Close the PNG file:
dev.off()

gck_vst_exp <- merge(ucec_prim_clinical_data["nmf_cluster"], as.data.frame(t(ucec_vst_data["ENSG00000106633", ])),
                     by = 'row.names', all = TRUE)
head(gck_vst_exp,3)

boxplot(ENSG00000106633~nmf_cluster, data = gck_vst_exp, xlab = "Clusters", ylab = "normalised expression",
        main = "GCK",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = FALSE)

mtn_vst_data = vst_norm(mtn_count_data)
head(mtn_vst_data,2)

mtn_gck_vst_exp <-merge(mtn_sample_type, as.data.frame(t(mtn_vst_data))["ENSG00000106633"],
                     by = 'row.names', all = TRUE)
head(mtn_gck_vst_exp,3)

ucec_all_vst_data <- vst_norm(ucec_exp_data)
ucec_all_vst_data <- ucec_all_vst_data[rownames(ucec_vst_data), ]
print(dim(ucec_all_vst_data))
head(ucec_all_vst_data,2)

merged_ucec_all_vst_data <- merge(as.data.frame(t(ucec_all_vst_data)), 
           ucec_prim_clinical_data["nmf_cluster"], 
           by = 'row.names', all.x = TRUE)
merged_ucec_all_vst_data$nmf_cluster[is.na(merged_ucec_all_vst_data$nmf_cluster)]<- "Solid_Tissue_Normal"
rownames(merged_ucec_all_vst_data) <- merged_ucec_all_vst_data$Row.names
merged_ucec_all_vst_data$Row.names <- NULL
head(merged_ucec_all_vst_data, 2)

table(merged_ucec_all_vst_data$nmf_cluster)

boxplot(ENSG00000130707~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "ASS1",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000106633~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "GCK",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000007350~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "TKTL1",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000036473~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "OTC",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000109193~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "SULT1E1",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000145692~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "BHMT",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000135069~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "PSAT1",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000178445~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "GLDC",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000021826~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "CPS1",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000118137~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "APOA1",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000142319~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "SLC6A3",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000134014~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "ELP3",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000172508~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "CARNS1",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

boxplot(ENSG00000141744~nmf_cluster, 
        data = merged_ucec_all_vst_data, xlab = "sample_type",
        ylab = "normalised expression",
        main = "PNMT",  col = c("#8470FF", "#F08080"),
        notch = TRUE, outline = TRUE)

library("readxl")
library("KEGGREST")
library(EnrichmentBrowser)
library(stringr)

hmr2_db_df <- read_excel('/home/akansha/downloads/metabolic_characterization/HMRdatabase2_00.xlsx')[, 
                                                                c("GENE ASSOCIATION", "SUBSYSTEM", "COMPARTMENT")]
print(dim(hmr2_db_df))
#head(hmr2_db_df,3)

print(table(is.na(hmr2_db_df['SUBSYSTEM'])))

print(table(is.na(hmr2_db_df['GENE ASSOCIATION'])))

hmr2_db_df <- hmr2_db_df[complete.cases(hmr2_db_df),]
hmr2_db_df$compart_pathway <- str_c(hmr2_db_df$COMPARTMENT, "_", hmr2_db_df$SUBSYSTEM)
print(dim(hmr2_db_df))
#head(hmr2_db_df, 3)

cat("\nTotal Number of pathways in HMR2: ", dim(unique(hmr2_db_df['SUBSYSTEM']))[1])
cat("\nTotal Number of compartment-wise pathways in HMR2: ", dim(unique(hmr2_db_df['compart_pathway']))[1])


hmr2_db <- list()
for (sys in unique(hmr2_db_df['compart_pathway'])$compart_pathway){
    hmr2_db[[sys]] <- c(unique(unlist(strsplit(c(filter(hmr2_db_df, SUBSYSTEM == substring(sys, 3,) & COMPARTMENT == substring(sys, 1,1)))[[1]], ';'))))

}
print(length(names(hmr2_db)))
hmr2_db[2]

hmr2_db_df <- read_excel('/home/akansha/downloads/metabolic_characterization/HMRdatabase2_00.xlsx')[, 
                                                                c("GENE ASSOCIATION", "SUBSYSTEM")]
print(dim(hmr2_db_df))
#head(hmr2_db_df,3)

print(table(is.na(hmr2_db_df['SUBSYSTEM'])))

print(table(is.na(hmr2_db_df['GENE ASSOCIATION'])))

hmr2_db_df <- hmr2_db_df[complete.cases(hmr2_db_df),]
print(dim(hmr2_db_df))
#head(hmr2_db_df, 3)

cat("\nTotal Number of pathways in HMR2: ", dim(unique(hmr2_db_df['SUBSYSTEM']))[1])

hmr2_db <- list()
for (sys in unique(hmr2_db_df['SUBSYSTEM'])$SUBSYSTEM){
    hmr2_db[[sys]] <- c(unique(unlist(strsplit(c(hmr2_db_df[hmr2_db_df['SUBSYSTEM'] == sys,
                                                            'GENE ASSOCIATION'])[[1]], ';'))))

}
print(length(names(hmr2_db)))
hmr2_db[2]

kegg_pathway_db <- getGenesets(org = "hsa", db = "kegg", cache = TRUE, return.type="list")
cat("Total Number of pathways in KEGG: ", length(kegg_pathway_db))

web_kegg_met_db <- read.csv("/home/akansha/downloads/metabolic_characterization/met_subtypes_nmf/kegg_metabolic_pathways.csv")
web_kegg_met_db$X <-NULL
print(dim(web_kegg_met_db))
#head(web_kegg_met_db, 8)

lst_common_met_pathways <- c()
for (pathway in web_kegg_met_db$kegg_Metabolic_pathways){
    com_path <- names(kegg_pathway_db)[grep(strsplit(pathway, '_')[[1]][1] , names(kegg_pathway_db))]
    lst_common_met_pathways <- c(lst_common_met_pathways, com_path)
}
print(length(lst_common_met_pathways))

kegg_met_pathway_db <-  kegg_pathway_db[lst_common_met_pathways]

cat("Total Number of metabolic pathways in KEGG: ", length(kegg_met_pathway_db))
kegg_met_pathway_db[2]

create_contigency_table <- function(n_genes_db, n_genes_pathway, 
                                    n_genes_interest_pathdb, n_genes_interest_pathway){
    
    ## Prepare a two-dimensional contingency table
    contingency.table <- data.frame(matrix(nrow=2, ncol=2))
    colnames(contingency.table) <- c("interesting_genes", "not_interesting_genes")
    rownames(contingency.table) <- c("enriched_pathway", "not_enriched_pathway")
    
    contingency.table["enriched_pathway", "interesting_genes"] <- n_genes_interest_pathway ## Number of marked genes in the selection
    contingency.table["not_enriched_pathway", "interesting_genes"] <- n_genes_interest_pathdb - n_genes_interest_pathway ## Number of non-marked genes in the selection
    contingency.table["enriched_pathway", "not_interesting_genes"] <- n_genes_pathway - n_genes_interest_pathway ## Number of marked genes outside of the selection
    contingency.table["not_enriched_pathway", "not_interesting_genes"] <- n_genes_db - n_genes_interest_pathdb - n_genes_pathway + n_genes_interest_pathway ## Number of non-marked genes in the selection
    #print(contingency.table)
    return (contingency.table)
}

print_contigency_table_marginal_sums <- function(contingency.table){
    ## Create a contingency table with marginal sums
    contingency.row.sum <- apply(contingency.table, 1, sum)
    contingency.col.sum <- apply(contingency.table, 2, sum)
    contingency.table.margins <- cbind(contingency.table, contingency.row.sum)
    contingency.table.margins <- rbind(contingency.table.margins, apply(contingency.table.margins, 2, sum))
    names(contingency.table.margins) <- c(names(contingency.table), "total")
    rownames(contingency.table.margins) <- c(rownames(contingency.table), "total")
    print(contingency.table.margins)
    
}

fischer_exact_test <- function(contigency_table){
    f_test_result <- fisher.test(x=contigency_table, alternative="greater")
    #print(f_test_result)
    return (f_test_result$p.value)
}

#fischer_exact_test(rbind(c(19, 592), c(40, 12937)))

print(fischer_exact_test(rbind(c(9, 52), c(255, 3449))))  # Glycine, serine
fischer_exact_test(rbind(c(12, 49), c(486, 3218)))  # Glycine, serine

gsea <- function(lst_genes_submittted, pathwayDB){
    lst_p_values <- c()
    lst_interesting_genes <- list()
    pathwayDB_genes <- unique(unlist(pathwayDB))
    #print(length(pathwayDB_genes))
    
    lst_genes_interest <- intersect(pathwayDB_genes, lst_genes_submittted)
    for (path in seq(1:length(names(pathwayDB)))){
        #print(path)
        pathway_genes <- pathwayDB[[path]]
        genes_of_interest_in_pathway <- intersect(pathway_genes, lst_genes_interest)
        #print(genes_of_interest_in_pathway)
        
        contegency_tb <- create_contigency_table(length(pathwayDB_genes), length(pathway_genes), 
                            length(lst_genes_interest), length(genes_of_interest_in_pathway))
        
        #print(names(pathwayDB)[path])
        #print_contigency_table_marginal_sums(contegency_tb)
        
        lst_p_values <- c(lst_p_values, fischer_exact_test(contegency_tb)) 
        lst_interesting_genes[[path]] <- gene_mapping[genes_of_interest_in_pathway,'external_gene_name']
    }
    adjp_values <- p.adjust(lst_p_values, method = p.adjust.methods[5], n = length(lst_p_values))
    print(length(lst_p_values))
    return (list(lst_p_values, adjp_values, lst_interesting_genes))
}

gsea_result_df <- function(db, p_values, adjp_values, gene_list){ 
    result_df <- data.frame("pathways" = names(db), "p_value" = p_values, 
                             "adj_p_value" = adjp_values)
    result_df$gene_lst <- gene_list
    result_df <- result_df[order(result_df$p_value),]
    result_df <- filter(result_df, p_value != 1)
    #result_df <- result_df[result_df["p_value"]< 0.05,]
    print(dim(result_df))
    return (result_df)
    
}

#print(length(names(hmr2_db)))
#hmr2_db <- hmr2_db[names(hmr2_db) %in% c('Isolated', 'Transport, nuclear',
#                            'Transport, endoplasmic reticular','Transport, mitochondrial',
#                              'Transport, peroxisomal', 'Transport, extracellular',
#                              'Miscellaneous', 'isolated', 'Transport, lysosomal', 
#                              'Transport, Golgi apparatus') == FALSE]
# Remove multiple list elements
#length(names(hmr2_db))

#clust1_hmr2_gsea <- gsea(clust1_genes, hmr2_db)
#clust1_hmr2_gsea_result_df <- gsea_result_df(hmr2_db, clust1_hmr2_gsea[[1]], clust1_hmr2_gsea[[2]])
#clust1_hmr2_gsea_result_df
#write.table(clust1_hmr2_gsea_result_df, file=paste(res_dir,'pathways_clust1_genes_hmr2.tsv', sep=''),
#            quote=F, row.names = FALSE, sep = "\t")

#clust2_hmr2_gsea <- gsea(clust2_genes, hmr2_db)
#clust2_hmr2_gsea_result_df <- gsea_result_df(hmr2_db, clust2_hmr2_gsea[[1]], clust2_hmr2_gsea[[2]])
#clust2_hmr2_gsea_result_df
#write.table(clust2_hmr2_gsea_result_df, file=paste(res_dir,'pathways_clust2_genes_hmr2.tsv', sep=''),
#            quote=F, row.names = FALSE, sep = "\t")

hmr2_gsea_up_clust2_1 <- gsea(up_regulated_genes1, hmr2_db)
hmr2_gsea_up_clust2_1_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_up_clust2_1[[1]], 
                                 hmr2_gsea_up_clust2_1[[2]], hmr2_gsea_up_clust2_1[[3]])

hmr2_gsea_up_clust2_1_result_df$gene_lst = as.character(hmr2_gsea_up_clust2_1_result_df$gene_lst)
head(hmr2_gsea_up_clust2_1_result_df, 3)

write.table(hmr2_gsea_up_clust2_1_result_df, file=paste(res_dir,
            'HMR2_pathway_enrichment_results/clust2_up_logFC_1.tsv', 
             sep=''),quote=F, row.names = FALSE, sep = "\t")

hmr2_gsea_down_clust2_1 <- gsea(down_regulated_genes1, hmr2_db)
hmr2_gsea_down_clust2_1_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_down_clust2_1[[1]], 
                                                    hmr2_gsea_down_clust2_1[[2]], hmr2_gsea_down_clust2_1[[3]])

hmr2_gsea_down_clust2_1_result_df$gene_lst = as.character(hmr2_gsea_down_clust2_1_result_df$gene_lst)
head(hmr2_gsea_down_clust2_1_result_df, 3)

write.table(hmr2_gsea_down_clust2_1_result_df, file=paste(res_dir,
            'HMR2_pathway_enrichment_results/clust2_down_logFC_1.tsv', 
             sep=''),quote=F, row.names = FALSE, sep = "\t")

hmr2_gsea_all_clust2_1 <- gsea(c(up_regulated_genes1, down_regulated_genes1), hmr2_db)
hmr2_gsea_all_clust2_1_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_all_clust2_1[[1]], 
                                                    hmr2_gsea_all_clust2_1[[2]], hmr2_gsea_all_clust2_1[[3]])

hmr2_gsea_all_clust2_1_result_df$gene_lst = as.character(hmr2_gsea_all_clust2_1_result_df$gene_lst)
head(hmr2_gsea_all_clust2_1_result_df, 3)

write.table(hmr2_gsea_all_clust2_1_result_df, file=paste(res_dir,
            'HMR2_pathway_enrichment_results/clust2_all_logFC_1.tsv', 
             sep=''),quote=F, row.names = FALSE, sep = "\t")

hmr2_gsea_up_mtn <- gsea(mtn_up_regulated_genes1, hmr2_db)
hmr2_gsea_up_mtn_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_up_mtn[[1]], 
                                 hmr2_gsea_up_mtn[[2]], hmr2_gsea_up_mtn[[3]])

hmr2_gsea_up_mtn_result_df$gene_lst = as.character(hmr2_gsea_up_mtn_result_df$gene_lst)
head(hmr2_gsea_up_mtn_result_df, 3)

#write.table(hmr2_gsea_up_mtn_result_df, file=paste(res_dir,
#            'HMR2_pathway_enrichment_results/mTN_up_logFC_1.tsv',sep=''),quote=F,row.names = FALSE, sep = "\t")

hmr2_gsea_down_mtn <- gsea(mtn_down_regulated_genes1, hmr2_db)
hmr2_gsea_down_mtn_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_down_mtn[[1]], 
                                 hmr2_gsea_down_mtn[[2]], hmr2_gsea_down_mtn[[3]])

hmr2_gsea_down_mtn_result_df$gene_lst = as.character(hmr2_gsea_down_mtn_result_df$gene_lst)
head(hmr2_gsea_down_mtn_result_df, 3)

#write.table(hmr2_gsea_down_mtn_result_df, file=paste(res_dir,
#            'HMR2_pathway_enrichment_results/mTN_down_logFC_1.tsv', 
#            sep=''),quote=F, row.names = FALSE, sep = "\t")

hmr2_gsea_all_mtn <- gsea(c(mtn_up_regulated_genes1, mtn_down_regulated_genes1), hmr2_db)
hmr2_gsea_all_mtn_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_all_mtn[[1]], 
                                 hmr2_gsea_all_mtn[[2]], hmr2_gsea_all_mtn[[3]])

hmr2_gsea_all_mtn_result_df$gene_lst = as.character(hmr2_gsea_all_mtn_result_df$gene_lst)
head(hmr2_gsea_all_mtn_result_df, 3)

#write.table(hmr2_gsea_all_mtn_result_df, file=paste(res_dir,
#            'HMR2_pathway_enrichment_results/mTN_all_logFC_1.tsv', 
#             sep=''),quote=F, row.names = FALSE, sep = "\t")

hmr2_gsea_up_clust1 <- gsea(clust1_up_regulated_genes1, hmr2_db)
hmr2_gsea_up_clust1_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_up_clust1[[1]], hmr2_gsea_up_clust1[[2]])
hmr2_gsea_up_clust1_result_df
#write.table(hmr2_gsea_up_clust1_result_df, file=paste(res_dir,'pathways_up_clust1_hmr2.tsv', sep=''),
#            quote=F, row.names = FALSE, sep = "\t")

hmr2_gsea_down_clust1 <- gsea(clust1_down_regulated_genes1, hmr2_db)
hmr2_gsea_down_clust1_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_down_clust1[[1]], hmr2_gsea_down_clust1[[2]])
hmr2_gsea_down_clust1_result_df
#write.table(hmr2_gsea_down_clust1_result_df, file=paste(res_dir,'pathways_down_clust1_hmr2.tsv', sep=''),
#            quote=F, row.names = FALSE, sep = "\t")

print_genes(hmr2_gsea_up_clust1_result_df$pathways)

hmr2_gsea_up_clust2 <- gsea(clust2_up_regulated_genes1, hmr2_db)
hmr2_gsea_up_clust2_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_up_clust2[[1]], hmr2_gsea_up_clust2[[2]])
hmr2_gsea_up_clust2_result_df
#write.table(hmr2_gsea_up_clust2_result_df, file=paste(res_dir,'pathways_up_clust2_hmr2.tsv', sep=''),
#            quote=F, row.names = FALSE, sep ="\t")

hmr2_gsea_down_clust2 <- gsea(clust2_down_regulated_genes1, hmr2_db)
hmr2_gsea_down_clust2_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_down_clust2[[1]], hmr2_gsea_down_clust2[[2]])
hmr2_gsea_down_clust2_result_df
#write.table(hmr2_gsea_down_clust2_result_df, file=paste(res_dir,'pathways_down_clust2_hmr2.tsv', sep=''),
#            quote=F, row.names = FALSE, sep = '\t' )

#print_genes(hmr2_gsea_down_clust2_result_df$pathways)

kegg_gsea_up_clust2_1 <- gsea(hmr2_genes_mapping[up_regulated_genes1,'GENE.ID.2'], kegg_met_pathway_db)
kegg_gsea_up_clust2_1_result_df <-gsea_result_df(kegg_met_pathway_db, kegg_gsea_up_clust2_1[[1]], kegg_gsea_up_clust2_1[[2]])
head(kegg_gsea_up_clust2_1_result_df, 5)

kegg_gsea_up_clust2_1 <- gsea(hmr2_genes_mapping[up_regulated_genes1,'GENE.ID.2'], kegg_pathway_db)
kegg_gsea_up_clust2_1_result_df <-gsea_result_df(kegg_pathway_db, kegg_gsea_up_clust2_1[[1]], kegg_gsea_up_clust2_1[[2]])
head(kegg_gsea_up_clust2_1_result_df, 10)

hmr2_gsea_uni_down_clust2 <- gsea(uniq_c2_down_genes, hmr2_db)
hmr2_gsea_uni_down_clust2_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_uni_down_clust2[[1]], 
                                                  hmr2_gsea_uni_down_clust2[[2]])
head(hmr2_gsea_uni_down_clust2_result_df, 4)

hmr2_gsea_uni_up_clust2 <- gsea(uniq_c2_up_genes, hmr2_db)
hmr2_gsea_uni_up_clust2_result_df <- gsea_result_df(hmr2_db, hmr2_gsea_uni_up_clust2[[1]], 
                                                  hmr2_gsea_uni_up_clust2[[2]])
head(hmr2_gsea_uni_up_clust2_result_df, 1)

tn_dds_res_tb$external_gene_name <- gene_mapping[rownames(tn_dds_res_tb),'external_gene_name']
tn_dds_res_tb$entrez_id <- hmr2_genes_mapping[rownames(tn_dds_res_tb),'GENE.ID.2']
tn_dds_res_tb$hmr2_gene_names <- hmr2_genes_mapping[rownames(tn_dds_res_tb),'SHORT.NAME']
#write.csv(tn_dds_res_tb, file=paste(res_dir,'TN_deseq2_result.csv', sep=''),quote=F)

clust1_dds_res_tb$external_gene_name <- gene_mapping[rownames(clust1_dds_res_tb),'external_gene_name']
clust1_dds_res_tb$entrez_id <- hmr2_genes_mapping[rownames(clust1_dds_res_tb),'GENE.ID.2']
clust1_dds_res_tb$hmr2_gene_names <- hmr2_genes_mapping[rownames(clust1_dds_res_tb),'SHORT.NAME']
#write.csv(clust1_dds_res_tb, file=paste(res_dir,'clust1_deseq2_result.csv', sep=''),quote=F)

clust2_dds_res_tb$external_gene_name <- gene_mapping[rownames(clust2_dds_res_tb),'external_gene_name']
clust2_dds_res_tb$entrez_id <- hmr2_genes_mapping[rownames(clust2_dds_res_tb),'GENE.ID.2']
clust2_dds_res_tb$hmr2_gene_names <- hmr2_genes_mapping[rownames(clust2_dds_res_tb),'SHORT.NAME']
#write.csv(clust2_dds_res_tb, file=paste(res_dir,'clust2_deseq2_result.csv', sep=''),quote=F)

prim_clinical_clusters <- read.csv(paste(res_dir,'prim_clinical_cluster_assign_data.csv', sep=''))
rownames(prim_clinical_clusters) <- prim_clinical_clusters$X
prim_clinical_clusters <- prim_clinical_clusters[,-1]
print(dim(prim_clinical_clusters))
table(prim_clinical_clusters$nmf_cluster)

pca_df3 = merge(t(ucec_vst_data) , prim_clinical_clusters["nmf_cluster"], by="row.names") 
row.names(pca_df3) = pca_df3$Row.names
pca_df3 = pca_df3[, -1]
print(dim(pca_df3))
#pca_df
pca_res3 <- prcomp( pca_df3[, - 3585], center = TRUE, scale. = TRUE)

pca_plot_file5 <- file.path(paste(res_dir, 'pca_2clust_all_met_plot_boundary.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(pca_plot_file5, height = 5, width = 5, units = 'in', res=600)

# Print your heatmap
autoplot(pca_res3, data = pca_df3, colour = 'nmf_cluster', frame = TRUE)

# Close the PNG file:
dev.off()

pca_plot_file6 <- file.path(paste(res_dir, 'pca_2clust_all_met_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(pca_plot_file6, height = 5, width = 5, units = 'in', res=600)

# Print your heatmap
autoplot(pca_res3, data = pca_df3, colour = 'nmf_cluster')

# Close the PNG file:
dev.off()

autoplot(pca_res3, data = pca_df3, colour = 'nmf_cluster', frame = TRUE)

pca_df4 = merge(t(ucec_vst_data[ucec_top_mad_genes,]) , prim_clinical_clusters["nmf_cluster"], by="row.names") 
row.names(pca_df4) = pca_df4$Row.names
pca_df4 = pca_df4[, -1]
print(dim(pca_df4))
pca_res4 <- prcomp( pca_df4[, -1001], center = TRUE, scale. = TRUE)

pca_plot_file7 <- file.path(paste(res_dir, 'pca_2clust_MAD_met_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(pca_plot_file7, height = 5, width = 5, units = 'in', res=600)

# Print your heatmap
autoplot(pca_res4, data = pca_df4, colour = 'nmf_cluster')

# Close the PNG file:
dev.off()

pca_plot_file8 <- file.path(paste(res_dir, 'pca_2clust_MAD_met_plot_boundary.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(pca_plot_file8, height = 5, width = 5, units = 'in', res=600)

# Print your heatmap
autoplot(pca_res4, data = pca_df4, colour = 'nmf_cluster', frame = TRUE)

# Close the PNG file:
dev.off()

autoplot(pca_res4, data = pca_df4, colour = 'nmf_cluster', frame = TRUE)

vst_exp_surv <- as.data.frame(t(ucec_vst_data[union(up_regulated_genes1, ucec_top_mad_genes),]))
print(dim(vst_exp_surv))
#vst_exp_surv

vst_exp_surv1 <- as.data.frame(t(ucec_vst_data[rownames(sig_res1),]))
print(dim(vst_exp_surv1))
#vst_exp_surv

res.cox <- coxph(Surv(ucec_prim_clinical_data$OS.time, 
                      ucec_prim_clinical_data$OS) ~ vst_exp_surv[,"ENSG00000106633"])
res.cox

summary(res.cox)

row = 1
mat = matrix(ncol = 6, nrow = 0)  
univariate_res  = data.frame(mat)
colnames(univariate_res) <- c("Gene", "coeficient_beta", "Hazard_Ratio", "wald.test", "p.value", "Ensembl_id")
for (col in colnames(vst_exp_surv)){
    res.cox <- coxph(Surv(ucec_prim_clinical_data$OS.time, 
                      ucec_prim_clinical_data$OS) ~ vst_exp_surv[,col])
    res.cox <- summary(res.cox)
    p.value<-  signif(res.cox$wald["pvalue"], digits=2)
    wald.test<-signif(res.cox$wald["test"], digits=2)
    beta<-signif(res.cox$coef[1], digits=2);#coeficient beta
    HR <-signif(res.cox$coef[2], digits=2);#exp(beta)
    res<-c(gene_mapping[col ,'external_gene_name'], beta, HR, wald.test, p.value, col)
    univariate_res[row, ] <- res
    row <- row + 1
    
}

univariate_res[, c("coeficient_beta", "Hazard_Ratio", "wald.test", "p.value")]  <- lapply(univariate_res[, 
                                c("coeficient_beta", "Hazard_Ratio", "wald.test", "p.value")] ,as.numeric)
sig_univariate_res <- filter(univariate_res, p.value < 0.05)
print(dim(sig_univariate_res))
head(sig_univariate_res,4)

row = 1
mat = matrix(ncol = 6, nrow = 0)  
univariate_res  = data.frame(mat)
colnames(univariate_res) <- c("Gene", "coeficient_beta", "Hazard_Ratio", "wald.test", "p.value", "Ensembl_id")
for (col in colnames(vst_exp_surv1)){
    res.cox <- coxph(Surv(ucec_prim_clinical_data$OS.time, 
                      ucec_prim_clinical_data$OS) ~ vst_exp_surv1[,col])
    res.cox <- summary(res.cox)
    p.value<-  signif(res.cox$wald["pvalue"], digits=2)
    wald.test<-signif(res.cox$wald["test"], digits=2)
    beta<-signif(res.cox$coef[1], digits=2);#coeficient beta
    HR <-signif(res.cox$coef[2], digits=2);#exp(beta)
    res<-c(gene_mapping[col ,'external_gene_name'], beta, HR, wald.test, p.value, col)
    univariate_res[row, ] <- res
    row <- row + 1
    
}

univariate_res[, c("coeficient_beta", "Hazard_Ratio", "wald.test", "p.value")]  <- lapply(univariate_res[, 
                                c("coeficient_beta", "Hazard_Ratio", "wald.test", "p.value")] ,as.numeric)
sig_univariate_res <- filter(univariate_res, p.value < 0.05)
print(dim(sig_univariate_res))
head(sig_univariate_res,4)

#write.csv(univariate_res, file=paste(res_dir,'DEGs_results/univariate_results.csv', sep=''),quote=F, row.names = FALSE)
#write.csv(sig_univariate_res, file=paste(res_dir,'DEGs_results/prognostic_DEGs_table.csv', sep=''),quote=F, row.names = FALSE)

table(sig_univariate_res$Hazard_Ratio > 1)

"LPCAT2" %in% print_genes(high_exp_better_surv_genes)

table(sig_univariate_res$coeficient_beta > 0)

high_exp_better_surv_genes <- filter(sig_univariate_res, coeficient_beta < 0)$Gene
high_exp_poor_surv_genes <- filter(sig_univariate_res, coeficient_beta > 0)$Gene
print(length(high_exp_better_surv_genes))
print(length(high_exp_poor_surv_genes))

high_exp_better_surv_genes

high_exp_poor_surv_genes

chdh = as.data.frame(t(ucec_vst_data["ENSG00000016391", ]))
summary(chdh)

chdh$binary<- ifelse(chdh$ENSG00000016391 > 9.171,"high","low")
#head(chdh,3)

# Fitting the survival model
ucec_surv_fun_chdh = survfit(Surv(ucec_prim_clinical_data$OS.time/365,
                             ucec_prim_clinical_data$OS)~ chdh$binary)
ucec_surv_fun_chdh

survdiff(Surv(ucec_prim_clinical_data$OS.time/365,ucec_prim_clinical_data$OS)~ chdh$binary)

plot(ucec_surv_fun_chdh, main = "KM Plot for CHDH expression", xlim = c(1,15), xlab = 'Time (in years)', ylab = 'Survival Probability', 
     col = c('red', 'blue'), las= 1, lwd = 2,  mark.time = TRUE, cex=1)
legend(1, 0.15, legend = c('high (n=271)', 'low (n=271)'),lwd = 2, col = c('red', 'blue'), bty = '', cex = 1)
text(10, 0.9,  'p-value = 0.04', cex=1.3) 

chat = as.data.frame(t(ucec_vst_data["ENSG00000070748", ]))
summary(chat)

chat$binary<- ifelse(chat$ENSG00000070748 > 3.354,"high","low")
#head(chdh,3)

# Fitting the survival model
surv_fun_chat = survfit(Surv(ucec_prim_clinical_data$OS.time/365,
                             ucec_prim_clinical_data$OS)~ chat$binary)
surv_fun_chat

survdiff(Surv(ucec_prim_clinical_data$OS.time/365,ucec_prim_clinical_data$OS)~ chat$binary)

plot(surv_fun_chat, main = "KM Plot for CHAT expression", xlim = c(1,15), xlab = 'Time (in years)', ylab = 'Survival Probability', 
     col = c('red', 'blue'), las= 1, lwd = 2,  mark.time = TRUE, cex=1)
legend(1, 0.15, legend = c('high (n=190)', 'low (n=352)'),lwd = 2, col = c('red', 'blue'), bty = '', cex = 1)
text(10, 0.9,  'p-value = 0.009', cex=1.3) 



covariates <- c("age", "sex",  "ph.karno", "ph.ecog", "wt.loss")
univ_formulas <- sapply(colnames(vst_exp_surv1[, sig_univariate_res$Ensembl_id[1:10]]),
                        function(x) paste('+ vst_exp_surv1$', x , sep = ""))
temp = ""
for (i in univ_formulas){
    temp = paste(temp, i)
}
print(temp)

multi_res.cox <- coxph(Surv(ucec_prim_clinical_data$OS.time, 
                      ucec_prim_clinical_data$OS) ~ vst_exp_surv1$ENSG00000002587 + vst_exp_surv1$ENSG00000002726 + vst_exp_surv1$ENSG00000004809 + vst_exp_surv1$ENSG00000008300 + vst_exp_surv1$ENSG00000010379 + vst_exp_surv1$ENSG00000021826 + vst_exp_surv1$ENSG00000027644 + vst_exp_surv1$ENSG00000036473 + vst_exp_surv1$ENSG00000058866 + vst_exp_surv1$ENSG00000064270 + vst_exp_surv1$ENSG00000064655 + vst_exp_surv1$ENSG00000066230 + vst_exp_surv1$ENSG00000067840 + vst_exp_surv1$ENSG00000068078 + vst_exp_surv1$ENSG00000069020 + vst_exp_surv1$ENSG00000069764 + vst_exp_surv1$ENSG00000070748 + vst_exp_surv1$ENSG00000070915 + vst_exp_surv1$ENSG00000071073 + vst_exp_surv1$ENSG00000072041 + vst_exp_surv1$ENSG00000072133 + vst_exp_surv1$ENSG00000072195 + vst_exp_surv1$ENSG00000072657 + vst_exp_surv1$ENSG00000074211 + vst_exp_surv1$ENSG00000077264 + vst_exp_surv1$ENSG00000079215 + vst_exp_surv1$ENSG00000080224 + vst_exp_surv1$ENSG00000081479 + vst_exp_surv1$ENSG00000086062 + vst_exp_surv1$ENSG00000087253 + vst_exp_surv1$ENSG00000087586 + vst_exp_surv1$ENSG00000091664 + vst_exp_surv1$ENSG00000091704 + vst_exp_surv1$ENSG00000095303 + vst_exp_surv1$ENSG00000095777 + vst_exp_surv1$ENSG00000100170 + vst_exp_surv1$ENSG00000100417 + vst_exp_surv1$ENSG00000100626 + vst_exp_surv1$ENSG00000101210 + vst_exp_surv1$ENSG00000101224 + vst_exp_surv1$ENSG00000101255 + vst_exp_surv1$ENSG00000101638 + vst_exp_surv1$ENSG00000103056 + vst_exp_surv1$ENSG00000103546 + vst_exp_surv1$ENSG00000105205 + vst_exp_surv1$ENSG00000105409 + vst_exp_surv1$ENSG00000105613 + vst_exp_surv1$ENSG00000106123 + vst_exp_surv1$ENSG00000106633 + vst_exp_surv1$ENSG00000106688 + vst_exp_surv1$ENSG00000107611 + vst_exp_surv1$ENSG00000107954 + vst_exp_surv1$ENSG00000108602 + vst_exp_surv1$ENSG00000108813 + vst_exp_surv1$ENSG00000109193 + vst_exp_surv1$ENSG00000111181 + vst_exp_surv1$ENSG00000111371 + vst_exp_surv1$ENSG00000111962 + vst_exp_surv1$ENSG00000112309 + vst_exp_surv1$ENSG00000112319 + vst_exp_surv1$ENSG00000112742 + vst_exp_surv1$ENSG00000112981 + vst_exp_surv1$ENSG00000113924 + vst_exp_surv1$ENSG00000114113 + vst_exp_surv1$ENSG00000114200 + vst_exp_surv1$ENSG00000115419 + vst_exp_surv1$ENSG00000116039 + vst_exp_surv1$ENSG00000116254 + vst_exp_surv1$ENSG00000116675 + vst_exp_surv1$ENSG00000117215 + vst_exp_surv1$ENSG00000117507 + vst_exp_surv1$ENSG00000118094 + vst_exp_surv1$ENSG00000118402 + vst_exp_surv1$ENSG00000119915 + vst_exp_surv1$ENSG00000124140 + vst_exp_surv1$ENSG00000124564 + vst_exp_surv1$ENSG00000124568 + vst_exp_surv1$ENSG00000125772 + vst_exp_surv1$ENSG00000128039 + vst_exp_surv1$ENSG00000128683 + vst_exp_surv1$ENSG00000129951 + vst_exp_surv1$ENSG00000130035 + vst_exp_surv1$ENSG00000130055 + vst_exp_surv1$ENSG00000130066 + vst_exp_surv1$ENSG00000130508 + vst_exp_surv1$ENSG00000130589 + vst_exp_surv1$ENSG00000130707 + vst_exp_surv1$ENSG00000130829 + vst_exp_surv1$ENSG00000130876 + vst_exp_surv1$ENSG00000131055 + vst_exp_surv1$ENSG00000132164 + vst_exp_surv1$ENSG00000132437 + vst_exp_surv1$ENSG00000132518 + vst_exp_surv1$ENSG00000132677 + vst_exp_surv1$ENSG00000132874 + vst_exp_surv1$ENSG00000132932 + vst_exp_surv1$ENSG00000133216 + vst_exp_surv1$ENSG00000133742 + vst_exp_surv1$ENSG00000133878 + vst_exp_surv1$ENSG00000134014 + vst_exp_surv1$ENSG00000135069 + vst_exp_surv1$ENSG00000135220 + vst_exp_surv1$ENSG00000135318 + vst_exp_surv1$ENSG00000135454 + vst_exp_surv1$ENSG00000135917 + vst_exp_surv1$ENSG00000136943 + vst_exp_surv1$ENSG00000138030 + vst_exp_surv1$ENSG00000138696 + vst_exp_surv1$ENSG00000138744 + vst_exp_surv1$ENSG00000138769 + vst_exp_surv1$ENSG00000139044 + vst_exp_surv1$ENSG00000139540 + vst_exp_surv1$ENSG00000140057 + vst_exp_surv1$ENSG00000140297 + vst_exp_surv1$ENSG00000140505 + vst_exp_surv1$ENSG00000141639 + vst_exp_surv1$ENSG00000141736 + vst_exp_surv1$ENSG00000141744 + vst_exp_surv1$ENSG00000141934 + vst_exp_surv1$ENSG00000142319 + vst_exp_surv1$ENSG00000142494 + vst_exp_surv1$ENSG00000142619 + vst_exp_surv1$ENSG00000143753 + vst_exp_surv1$ENSG00000143797 + vst_exp_surv1$ENSG00000143882 + vst_exp_surv1$ENSG00000144278 + vst_exp_surv1$ENSG00000145626 + vst_exp_surv1$ENSG00000145675 + vst_exp_surv1$ENSG00000145692 + vst_exp_surv1$ENSG00000146039 + vst_exp_surv1$ENSG00000149150 + vst_exp_surv1$ENSG00000149452 + vst_exp_surv1$ENSG00000149527 + vst_exp_surv1$ENSG00000149599 + vst_exp_surv1$ENSG00000152270 + vst_exp_surv1$ENSG00000154027 + vst_exp_surv1$ENSG00000154269 + vst_exp_surv1$ENSG00000154928 + vst_exp_surv1$ENSG00000155897 + vst_exp_surv1$ENSG00000156219 + vst_exp_surv1$ENSG00000156475 + vst_exp_surv1$ENSG00000156689 + vst_exp_surv1$ENSG00000157087 + vst_exp_surv1$ENSG00000158089 + vst_exp_surv1$ENSG00000158104 + vst_exp_surv1$ENSG00000159714 + vst_exp_surv1$ENSG00000160200 + vst_exp_surv1$ENSG00000161798 + vst_exp_surv1$ENSG00000161896 + vst_exp_surv1$ENSG00000162174 + vst_exp_surv1$ENSG00000162571 + vst_exp_surv1$ENSG00000162669 + vst_exp_surv1$ENSG00000163558 + vst_exp_surv1$ENSG00000163629 + vst_exp_surv1$ENSG00000164100 + vst_exp_surv1$ENSG00000164181 + vst_exp_surv1$ENSG00000164574 + vst_exp_surv1$ENSG00000165695 + vst_exp_surv1$ENSG00000166840 + vst_exp_surv1$ENSG00000167103 + vst_exp_surv1$ENSG00000167311 + vst_exp_surv1$ENSG00000168032 + vst_exp_surv1$ENSG00000169418 + vst_exp_surv1$ENSG00000170324 + vst_exp_surv1$ENSG00000170439 + vst_exp_surv1$ENSG00000170482 + vst_exp_surv1$ENSG00000171094 + vst_exp_surv1$ENSG00000172264 + vst_exp_surv1$ENSG00000172828 + vst_exp_surv1$ENSG00000172955 + vst_exp_surv1$ENSG00000173083 + vst_exp_surv1$ENSG00000173627 + vst_exp_surv1$ENSG00000174156 + vst_exp_surv1$ENSG00000175063 + vst_exp_surv1$ENSG00000175229 + vst_exp_surv1$ENSG00000176463 + vst_exp_surv1$ENSG00000177108 + vst_exp_surv1$ENSG00000178445 + vst_exp_surv1$ENSG00000178538 + vst_exp_surv1$ENSG00000179455 + vst_exp_surv1$ENSG00000180176 + vst_exp_surv1$ENSG00000180537 + vst_exp_surv1$ENSG00000180767 + vst_exp_surv1$ENSG00000180815 + vst_exp_surv1$ENSG00000181409 + vst_exp_surv1$ENSG00000181856 + vst_exp_surv1$ENSG00000182050 + vst_exp_surv1$ENSG00000182333 + vst_exp_surv1$ENSG00000182870 + vst_exp_surv1$ENSG00000183023 + vst_exp_surv1$ENSG00000183196 + vst_exp_surv1$ENSG00000183317 + vst_exp_surv1$ENSG00000183654 + vst_exp_surv1$ENSG00000183921 + vst_exp_surv1$ENSG00000184154 + vst_exp_surv1$ENSG00000185527 + vst_exp_surv1$ENSG00000185974 + vst_exp_surv1$ENSG00000186204 + vst_exp_surv1$ENSG00000186526 + vst_exp_surv1$ENSG00000187210 + vst_exp_surv1$ENSG00000187714 + vst_exp_surv1$ENSG00000187980 + vst_exp_surv1$ENSG00000188338 + vst_exp_surv1$ENSG00000196090 + vst_exp_surv1$ENSG00000196368 + vst_exp_surv1$ENSG00000196632 + vst_exp_surv1$ENSG00000196839 + vst_exp_surv1$ENSG00000197142 + vst_exp_surv1$ENSG00000197208 + vst_exp_surv1$ENSG00000197241 + vst_exp_surv1$ENSG00000197977 + vst_exp_surv1$ENSG00000198099 + vst_exp_surv1$ENSG00000198400 + vst_exp_surv1$ENSG00000198842 + vst_exp_surv1$ENSG00000198910 + vst_exp_surv1$ENSG00000203805 + vst_exp_surv1$ENSG00000204099 + vst_exp_surv1$ENSG00000211448 + vst_exp_surv1$ENSG00000214617 + vst_exp_surv1$ENSG00000241635 + vst_exp_surv1$ENSG00000244067 + vst_exp_surv1$ENSG00000250305 + vst_exp_surv1$ENSG00000257335 + vst_exp_surv1$ENSG00000257594 + vst_exp_surv1$ENSG00000259075)

## 183th patient has 0 "OS.time", which is suitabe as input in glmnet package
fit = glmnet(as.matrix(vst_exp_surv[-c(183), sig_univariate_res$Ensembl_id]), # X matrix
                   Surv(ucec_prim_clinical_data[-c(183),]$OS.time, ucec_prim_clinical_data[-c(183), ]$OS), 
                   alpha = 1, # lasso: alpha = 1; ridge: alpha=0
                   family = "cox" ) # specify Cox PH model)
#plot(fit)
# print(fit)
# fit$lambda

# print(fit)

dim(summary(coef(fit, s = 0.018960)))

cv.fit = cv.glmnet(as.matrix(vst_exp_surv[-c(183), sig_univariate_res$Ensembl_id]), # X matrix
                   Surv(ucec_prim_clinical_data[-c(183),]$OS.time, ucec_prim_clinical_data[-c(183), ]$OS), 
                   alpha = 1, # lasso: alpha = 1; ridge: alpha=0
                   family = "cox" , type.measure = "C") # specify Cox PH model)

plot(cv.fit)

plot(cv.fit)

print(log(cv.fit$lambda.min)) # in log-scale
print(cv.fit$lambda.min)
print(cv.fit$lambda.1se)

print(log(cv.fit$lambda.min)) # in log-scale
print(cv.fit$lambda.min)
print(cv.fit$lambda.1se)

est.coef = coef(cv.fit, s = cv.fit$lambda.min) # returns the p length coefficient vector
                                               # of the solution corresponding to lambda 
active.k = which(est.coef != 0)
print(length(active.k))
active.k.vals = est.coef[active.k]

multi_res <- as.data.frame(summary(coef(fit, s = 0.01895937)))
print(dim(multi_res))

imp_coef <- c()
for (val in multi_res$i){
    imp_coef <- c(imp_coef, dimnames(est.coef)[[1]][val])
}
print(length(imp_coef))

multi_res$Genes <- gene_mapping[imp_coef, "external_gene_name"]
head(multi_res)

table(multi_res$x >0)

filter(multi_res, x<0)$Genes

filter(multi_res, x>0)$Genes

# print_genes(filter(multi_res, x>0)$Genes)
