library(GEOquery)
library(limma)
library(umap)
library(NMF)
library(ggfortify)
library(dplyr)

hmr2_gene_file = '/home/akansha/downloads/metabolic_characterization/HMRdatabase2_00.xlsx - GENES.csv'
hmr2_genes = read.csv(hmr2_gene_file)
hmr2_genes <- unique(hmr2_genes$SHORT.NAME)
length(hmr2_genes)

# load series and platform data from GEO
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 4)

gset <- getGEO("GSE17025", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# gset 
print(dim(exprs(gset)))
print(dim(pData(gset)))
colnames(pData(gset))

colnames(pData(gset))[40] <- "grade"
colnames(pData(gset))[41] <- "histology"

table(gset$'grade')

table(gset$"histology")

ex <- exprs(gset)

# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 1)] <- NaN
  ex <- log2(ex) }

tumor_samples_id <- rownames(filter(pData(gset), grade != "n/a"))
tumor_exp <- as.data.frame(ex[, tumor_samples_id])
print(dim(tumor_exp))

# Remove nan values 
tumor_exp <- na.omit(tumor_exp)

print(dim(tumor_exp))
head(tumor_exp, 4)

apply(is.na(tumor_exp), 1, which)

# Check if dataframe has negative values or not
# tumor_exp[tumor_exp<0,]
has.neg <- apply(tumor_exp, 1, function(row) any(row < 0))
print(which(has.neg))

dir_n <- '/home/akansha/downloads/metabolic_characterization/met_subtypes_nmf/ucec_nmf_results_542sample_1000genes/'
probe_gene_mapping <- read.csv(paste(dir_n, "validation/GSE17025_probe_gene_mapping.csv", sep = ""))
rownames(probe_gene_mapping) <- probe_gene_mapping$ID
print(dim(probe_gene_mapping))
head(probe_gene_mapping, 2)

mad_genes <- read.csv(paste(dir_n, "mad_genes.csv", sep = ""))
print(dim(mad_genes))
head(mad_genes, 2)

mad_probes <-  filter(probe_gene_mapping, Gene.symbol %in% mad_genes$external_gene_name )
print(length(unique(mad_probes$Gene.symbol)))
print(dim(mad_probes))
head(mad_probes,3)

apply(is.na(mad_probes), 1, which)

# No of duplicated genes
print(length(unique(mad_probes[duplicated(mad_probes["Gene.symbol"]),]$Gene.symbol)))
mad_probes <- mad_probes[!duplicated(mad_probes["Gene.symbol"]),]
print(dim(mad_probes))
head(mad_probes,3)

mad_tumor_exp <- tumor_exp[intersect(mad_probes$ID, rownames(tumor_exp)),]
print(dim(mad_tumor_exp))
head(mad_tumor_exp,3)

merged_df <- merge(as.data.frame(t(mad_tumor_exp)), pData(gset)[c("histology","grade")], by = 'row.names', all = FALSE)
rownames(merged_df) <- merged_df$Row.names
merged_df$Row.names <- NULL
print(dim(merged_df))
head(merged_df, 3)

apply(is.na(mad_tumor_exp), 1, which)

nmf_run <- function(data, rank , method = 'brunet', seed_meth = "random", n_run = 10){
    res <- nmf(data, rank , method, seed= seed_meth, nrun = n_run, .options='t')
    return(res)
}

set.seed(1) 
nmf_result <- nmf_run( mad_tumor_exp , 2, method = 'offset', seed_meth ='random', n_run = 50)
nmf_result_summary <- t(as.data.frame(summary(nmf_result)))
nmf_result_summary

cp <- consensusmap(nmf_result, annCol=  merged_df[c("grade", "histology")],tracks =c(), Rowv = FALSE)

library("RColorBrewer")

consensusmap(ucec_best_res,
        annCol=ucec_clinical_data[colnames(ucec_prim_exp_data),][c('clinical_stage', 'histological_type', 'histological_grade')],
        annRow=ucec_clinical_data[colnames(ucec_prim_exp_data),][c('clinical_stage', 'histological_type', 'histological_grade')],
        annColors = list( c("darkgreen", "blue4", "orange", "red"), rainbow(3),   brewer.pal(n = 4, name = "Set1")),
        main='', 
        tracks=c(),
        info = FALSE,
        Rowv = FALSE,
        labRow = NA,
        labCol = NA,
        #color = "Blues",
        fontsize = 14
        ) 

options(repr.plot.width = 7, repr.plot.height = 5)
consensusmap(nmf_result, 
             annCol=  merged_df[c("grade", "histology")],
             annRow=  merged_df[c("grade", "histology")],
             annColors = list( brewer.pal(n = 3, name = "Set1"), brewer.pal(n = 4, name = "Dark2")),
             main='', 
        tracks=c(),
        info = FALSE,
        Rowv = FALSE,
        labRow = NA,
        labCol = NA,
        #color = "Blues",
        fontsize = 14)

consensus_plot_file <- file.path(paste('validation_consensus_plot.tiff', sep=''))

# Open a PNG file - width and height arguments control the size of the output
tiff(consensus_plot_file, height = 5, width = 7, units = 'in', res=300)

# Print your heatmap
consensusmap(nmf_result, 
             annCol=  merged_df[c("grade", "histology")],
             annRow=  merged_df[c("grade", "histology")],
             annColors = list( brewer.pal(n = 3, name = "Set1"), brewer.pal(n = 4, name = "Dark2")),
             main='', 
        tracks=c(),
        info = FALSE,
        Rowv = FALSE,
        labRow = NA,
        labCol = NA,
        #color = "Blues",
        fontsize = 14)

# Close the PNG file:
dev.off()

merged_df$nmf_cluster <- NA
merged_df[lapply(cut(cp$Rowv,0.5)$lower, function(l)rapply(l,function(i) i))[[1]], "nmf_cluster"] <- "cluster_1"
merged_df[lapply(cut(cp$Rowv,0.5)$lower, function(l)rapply(l,function(i) i))[[2]], "nmf_cluster"] <- "cluster_2" 
                                                           
merged_df$subtypes <- NA
merged_df[lapply(cut(cp$Rowv,0.5)$lower, function(l)rapply(l,function(i) i))[[1]], "subtypes"] <- "Metabolic_subtype-1"
merged_df[lapply(cut(cp$Rowv,0.5)$lower, function(l)rapply(l,function(i) i))[[2]], "subtypes"] <- "Metabolic_subtype-2"
                                                           
print(dim(merged_df))
head(merged_df, 3)

table(merged_df["nmf_cluster"])

table(merged_df["subtypes"])

table(merged_df[c("nmf_cluster", "histology")])

table(merged_df[c("nmf_cluster", "grade")])

get_probe_id <- function(gene_name){
    ids <- filter(probe_gene_mapping, Gene.symbol == gene_name)$ID
    return (ids)
}
get_probe_id("GLS")

copy_merged_df <- data.frame(merged_df)
colnames(copy_merged_df) <- c(filter(mad_probes, ID %in% colnames(merged_df))$Gene.symbol, "histology",
                              "grade", "nmf_cluster")
print(dim(copy_merged_df))
head(copy_merged_df,3)

gene_interest1 <- c("GLS", "ASS1", "PNMT", "ERBB2", "GLDC", "PSAT1", "SULT1E1", "UGT1A1",
                  "UGT3A1", "CYP1A2", "BHMT", "CBS", "TKTL1")
print(length(gene_interest1))

gene_intersect <- intersect(gene_interest1, colnames(copy_merged_df))
boxplot_df <- stack(copy_merged_df[gene_intersect])
boxplot_df$grade <- rep(copy_merged_df$grade, length(gene_intersect))
boxplot_df$nmf_cluster <- rep(copy_merged_df$nmf_cluster, length(gene_intersect))
boxplot_df$histology <- rep(copy_merged_df$histology, length(gene_intersect))
colnames(boxplot_df) = c("normalised_expression", "Genes", "grade", "nmf_cluster", "histology")
print(dim(boxplot_df))
head(boxplot_df, 4)

gene_intersect

ggplot(boxplot_df, aes(x= Genes, y = normalised_expression, fill = nmf_cluster)) + geom_boxplot() +
facet_wrap(~Genes, scale="free")

ggplot(boxplot_df, aes(x= Genes, y = normalised_expression, fill = grade)) + geom_boxplot() +
facet_wrap(~Genes, scale="free")

ggplot(boxplot_df, aes(x= Genes, y = normalised_expression, fill = histology)) + geom_boxplot() +
facet_wrap(~Genes, scale="free")

#heatmap(as.matrix(t(copy_merged_df[, intersect(gene_interest1, colnames(copy_merged_df)) ])), )

# set up design matrix
design <- model.matrix(~nmf_cluster+0, merged_df["nmf_cluster"])
colnames(design) <- c("cluster_1", "cluster_2")
print(dim(design))
head(design)

limma_exp_df <- tumor_exp[filter(probe_gene_mapping, Gene.symbol %in% hmr2_genes )$ID,]
print(dim(limma_exp_df))
head(limma_exp_df, 3)

# fit linear model
fit <- lmFit(limma_exp_df, design)  

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste("cluster_2","-", "cluster_1",sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

results <- decideTests(fit2, p.value=0.05)
vennDiagram(results) 

summary(results)

table(results)

ct <- 1
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
  highlight=length(which(results[,ct]!=0)), names=rep('+', nrow(fit2)))

# table of  significant genes (adj-p-value < 0.05)
sig_genes_res <- topTable(fit2,  number= 2769)
sig_genes_res <- merge(sig_genes_res, probe_gene_mapping,  by = "row.names")
rownames(sig_genes_res) <- sig_genes_res$Row.names
sig_genes_res$Row.names <- NULL
dim(sig_genes_res)
head(sig_genes_res, 4)
# write.table(tT, file=stdout(), row.names=F, sep="\t")

length(unique(sig_genes_res$Gene.symbol))

unq_sig_genes_res <- data.frame()
for (gene in unique(sig_genes_res$Gene.symbol)){
    #print(gene)
    temp_df <- filter(sig_genes_res, Gene.symbol == gene )
    #print(dim(temp_df)[1])
    if (dim(temp_df)[1] == 1){
        unq_sig_genes_res <- rbind(unq_sig_genes_res, temp_df)
    }
    else{
        unq_sig_genes_res <- rbind(unq_sig_genes_res, temp_df[which.min(temp_df$adj.P.Val),])
    }
}

print(dim(unq_sig_genes_res))
head(unq_sig_genes_res,4)

up_genes <- filter(unq_sig_genes_res, logFC> 0)
print(dim(up_genes))

down_genes <- filter(unq_sig_genes_res, logFC < 0)
print(dim(down_genes))

intersect(unq_sig_genes_res$Gene.symbol, gene_interest1)

intersect(up_genes$Gene.symbol, gene_interest1)

intersect(down_genes$Gene.symbol, gene_interest1)

print_genes <- function(lst_genes){
    for (i in lst_genes){
        cat(noquote(i), sep='\n')
    }
}

# print_genes(up_genes$Gene.symbol)
# print_genes(down_genes$Gene.symbol)
# print_genes(filter(up_genes, logFC>0.5)$Gene.symbol)
# print_genes(filter(down_genes, abs(logFC)>0.5)$Gene.symbol)
# print_genes(filter(up_genes, abs(logFC)>0.6)$Gene.symbol)
# print_genes(filter(down_genes, abs(logFC)>0.6)$Gene.symbol)

length(filter(up_genes, abs(logFC)>0.6)$Gene.symbol)

length(filter(down_genes, abs(logFC)>0.6)$Gene.symbol)

filter(unq_sig_genes_res, Gene.symbol %in% intersect(unq_sig_genes_res$Gene.symbol, gene_interest1))

gene_interest2 <- c("GCK", "TKTL1", "GLDC", "AGXT", "BHMT","ERBB2", "CBS", "GLS", "CYP1A2", "ELOVL2",
                   "UGT3A1", "SLC38A1")
gi2_res <- filter(unq_sig_genes_res, Gene.symbol %in% intersect(unq_sig_genes_res$Gene.symbol, gene_interest2))
gi2_res

gi2_exp_pheno_df <- merge(as.data.frame(t(limma_exp_df[rownames(gi2_res),])), 
                          merged_df[c("histology","grade", "nmf_cluster", "subtypes")], 
                         by = 'row.names', all = FALSE )
rownames(gi2_exp_pheno_df) <- gi2_exp_pheno_df$Row.names
gi2_exp_pheno_df$Row.names <- NULL
colnames(gi2_exp_pheno_df) <- c(gene_interest2, "histology","grade", "nmf_cluster", "subtypes")
print(dim(gi2_exp_pheno_df))
head(gi2_exp_pheno_df, 5)

boxplot_df2 <- stack(gi2_exp_pheno_df[gene_interest2])
boxplot_df2$grade <- rep(gi2_exp_pheno_df$grade, length(gene_interest2))
boxplot_df2$nmf_cluster <- rep(gi2_exp_pheno_df$nmf_cluster, length(gene_interest2))
boxplot_df2$subtypes <- rep(gi2_exp_pheno_df$subtypes, length(gene_interest2))
boxplot_df2$histology <- rep(gi2_exp_pheno_df$histology, length(gene_interest2))
colnames(boxplot_df2) = c("Normalized_expression", "Genes", "grade", "nmf_cluster", "subtypes", "histology")
print(dim(boxplot_df2))
head(boxplot_df2, 4)

ggplot(boxplot_df2, aes(x= Genes, y = Normalized_expression, fill = nmf_cluster)) + geom_boxplot()

library(tidyverse)
library(repr)

options(repr.plot.width = 10, repr.plot.height = 10 )
ggplot(boxplot_df2, aes(x= Genes, y = Normalized_expression, fill = subtypes)) + geom_boxplot() +
facet_wrap(~Genes, scale="free") +
theme(
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_rect(colour = "black", size = 0.6),
    axis.text.x.bottom = element_blank(),
    legend.title=element_blank(),
    axis.text = element_text(face="bold", size = 12),
    legend.text = element_text( size = 16, margin = margin(r = 30, unit = "pt")),
    legend.position="top",
    axis.title = element_text(size = 18)
  ) 

boxplot_plot_file <- file.path(paste( 'Boxplot-GSE17025-upDEGS_revision.tiff', sep =''))

# Open a PNG file - width and height arguments control the size of the output
tiff(boxplot_plot_file, height = 8, width = 8, units = 'in', res= 300)


# Print your heatmap
ggplot(boxplot_df2, aes(x= Genes, y = Normalized_expression, fill = subtypes)) + geom_boxplot() +
facet_wrap(~Genes, scale="free") +
theme(
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(colour = "black", size = 0.6),
    axis.text.x.bottom = element_blank(),
    legend.title=element_blank(),
    axis.text = element_text(face="bold", size = 12),
    legend.text = element_text( size = 16, margin = margin(r = 30, unit = "pt")),
    legend.position="top",
    axis.title = element_text(size = 18)
  ) 

# Close the PNG file:
dev.off()

library(vcd)

ucec_hist_type_f_test <- fisher.test(table(gi2_exp_pheno_df$histology, 
                                           gi2_exp_pheno_df$nmf_cluster))
print(ucec_hist_type_f_test)

assocstats(table(gi2_exp_pheno_df$histology, gi2_exp_pheno_df$nmf_cluster))

ucec_hist_grade_f_test <- fisher.test(table(gi2_exp_pheno_df$grade, 
                                           gi2_exp_pheno_df$nmf_cluster))
print(ucec_hist_grade_f_test)

assocstats(table(gi2_exp_pheno_df$grade, gi2_exp_pheno_df$nmf_cluster))


