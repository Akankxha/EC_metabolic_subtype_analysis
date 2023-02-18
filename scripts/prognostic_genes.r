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
