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

ucec_clinical_data <- read_clinical_file(paste(dir_nam, files_nam[2], sep=''))
ucec_exp_data <- read_exp_file(paste(dir_nam, files_nam[1], sep=''))
ucec_prim_clinical_data = read.csv(ucec_prim_clinical_data_file)
head(ucec_prim_clinical_data, 5)
print(dim(ucec_prim_clinical_data))

#rownames(ucec_prim_clinical_data) <- ucec_prim_clinical_data$X
#ucec_prim_clinical_data$X <- NULL

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
