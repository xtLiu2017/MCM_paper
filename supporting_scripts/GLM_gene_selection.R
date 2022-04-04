# Author: Xiaotong Liu (liux4215@umn.edu)
# This first part of the scripts is to fit a GLM model for each gene to obtain the mean estimates from read count data
# This first part is probably the only part in the folder '/MCM_paper' that cannot run smoothly. Because
# I deleted the raw read count data (GSE88798) in order to avoid using others' data directly. 
# If you need to run this script, please go to NCBI to find the data first and run this script. 

setwd("~/Documents/MCM_paper/data")
Mydata1 = read.table('GSE88798_ReadCountTable_M001_348.txt',header = TRUE, row.names = 1)
Mydata2 = read.table('GSE88798_ReadCountTable_M349_366.txt',header = TRUE, row.names = 1)
gene_names1 <- rownames(Mydata1)
gene_names2 <- rownames(Mydata2)
gene_names <- intersect(gene_names1,gene_names2)
libraries1 <- colnames(Mydata1)
libraries2 <- colnames(Mydata2)
treatment <- c()
geno <- c()
time <- c()
for (mylib in c(libraries1,libraries2)){
  treatment <- c(treatment, strsplit(mylib,split = "_")[[1]][4])
  geno <- c(geno,strsplit(mylib,split = "_")[[1]][2])
  time <- c(time, strsplit(mylib,split = "_")[[1]][5])
}
Mydata <- cbind(Mydata1[gene_names,], Mydata2[gene_names,])
Mydata <- Mydata + 1
raw_genes_in_Ken_transcriptome = rownames(Mydata)
save(raw_genes_in_Ken_transcriptome,file = "~/Documents/MCM_paper/data/raw_genes_in_Ken_transriptome.Rdata")
#no 0
cval.90.Ken = apply(Mydata, 2, function(x) {
  quantile(x, probs = 0.9)
})
add.v = round((cval.90.Ken-1)/100) -1 # -1 because we already added 1 to all (line 18).
Mydata.add <- Mydata + data.frame(matrix(rep(add.v, ea=nrow(Mydata)), nrow=nrow(Mydata)))
cval.90.Ken.add = apply(Mydata.add, 2, function(x) {
  quantile(x, probs = 0.9)
})
offset.cval.90.add = log(cval.90.Ken.add/max(cval.90.Ken.add))  # as 23 added to 2444.
top.15.val.per.gene = apply(as.matrix(Mydata.add), 1, function(x) sort(x, decreasing=T)[15])
kept.genes = gene_names[top.15.val.per.gene > 40]
# mean(add.v) + 25
# 18442
# model fit, for the col data set (AvrRpt2, EV, mock)
# takes some time ~ 30 min
library(MASS)
glm_fixef = list() 
AIC_list <- c()
for (mygene in kept.genes) {
  gene_exp = as.numeric(Mydata.add[mygene,geno == "Col" & treatment != 'AvrRpm1'])
  mytreat <- as.factor(treatment[geno == "Col" & treatment != 'AvrRpm1'])
  mygeno <- as.factor(geno[geno == "Col" & treatment != 'AvrRpm1'])
  mytime <- as.factor(time[geno ==  "Col" & treatment != 'AvrRpm1'])
  myoffset <- offset.cval.90.add[geno == "Col" & treatment != 'AvrRpm1']
  glm_model = tryCatch({glm.nb(gene_exp ~ -1 + mytreat:mytime + offset(myoffset),
                               link=log)},error = function(err){
                                 return(2)
                               })
  if (inherits(glm_model,"numeric") != T){
    AIC_list <- c(AIC_list,AIC(glm_model))
    glm_coef = summary(glm_model)$coef
    glm.se = glm_coef[,2]
    glm_fixef[[mygene]] = cbind(Estimate=glm_coef[,1], Std.Error=glm.se) /log(2) 
  }else{
    AIC_list <- c(AIC_list,NA)
    glm_fixef[[mygene]] = cbind(Estimate=0, Std.Error=1) /log(2) 
  }
}
save(glm_fixef,AIC_list,file = "~/Documents/MCM_paper/data/glm_fixef_Ken_WT_with_pseudocounts.Rdata")

# This second part is to select genes that are significantly induced at at least one time point
# after AvrRpt2 treatment

library(qvalue)    # Storey FDR

# the directory where I can obtain the data
setwd("~/Documents/MCM_paper/data")
# GLM fitting of all genes
load("glm_fixef_Ken_WT_with_pseudocounts.Rdata")
# mock treatment labels
Mock_label = c("mytreatmock:mytime01h","mytreatmock:mytime02h","mytreatmock:mytime03h",
               "mytreatmock:mytime04h","mytreatmock:mytime06h","mytreatmock:mytime09h",
               "mytreatmock:mytime12h","mytreatmock:mytime16h","mytreatmock:mytime20h",
               "mytreatmock:mytime24h")
# AvrRpt2 treatment labels (Pto AvrRpt2)
AvrRpt2_label = c("mytreatAvrRpt2:mytime01h","mytreatAvrRpt2:mytime02h","mytreatAvrRpt2:mytime03h",
                  "mytreatAvrRpt2:mytime04h","mytreatAvrRpt2:mytime06h","mytreatAvrRpt2:mytime09h",
                  "mytreatAvrRpt2:mytime12h","mytreatAvrRpt2:mytime16h","mytreatAvrRpt2:mytime20h",
                  "mytreatAvrRpt2:mytime24h")
# Empty vector labels (Pto)
EV_label = c("mytreatEV:mytime01h","mytreatEV:mytime02h","mytreatEV:mytime03h","mytreatEV:mytime04h",
             "mytreatEV:mytime06h","mytreatEV:mytime09h","mytreatEV:mytime12h","mytreatEV:mytime16h",
             "mytreatEV:mytime20h","mytreatEV:mytime24h")
# compare the glm fit difference between AvrRpt2 and mock, between EV and mock
# Note that the difference is on log2 scale
AvrRpt2andEV_mock_diff <- lapply(glm_fixef,function(x){
  AvrRpt2_means <- x[AvrRpt2_label,1];AvrRpt2_stds <- x[AvrRpt2_label,2]
  EV_means <- x[EV_label,1];EV_stds <- x[EV_label,2]
  mock_means <- x[Mock_label,1];mock_stds <- x[Mock_label,2]
  AvrRpt2_diff_mean <- AvrRpt2_means - mock_means
  EV_diff_mean <- EV_means - mock_means
  AvrRpt2_diff_std <- sqrt(AvrRpt2_stds ^ 2 + mock_stds ^ 2)
  EV_diff_std <- sqrt(EV_stds ^ 2 + mock_stds ^ 2)
  AvrRpt2_diff_p <- 2 * pnorm(abs(AvrRpt2_diff_mean)/AvrRpt2_diff_std, lower.tail = F)
  EV_diff_p <- 2 * pnorm(abs(EV_diff_mean)/EV_diff_std, lower.tail = F)
  
  diff <- c(AvrRpt2_diff_mean,EV_diff_mean)
  std <- c(AvrRpt2_diff_std,EV_diff_std)
  p_values <- c(AvrRpt2_diff_p,EV_diff_p)
  return(cbind(diff,std,p_values))
})
# get all p values
p_AvrRpt2andEV_mock_diff = unlist(lapply(AvrRpt2andEV_mock_diff, '[', ,3))
# get all q values
q_AvrRpt2andEV_mock_diff = qvalue(p_AvrRpt2andEV_mock_diff)$qvalues
# pvalue cutoff for FDR = 0.05
p_cutoff = sort(p_AvrRpt2andEV_mock_diff)[sum(q_AvrRpt2andEV_mock_diff < 0.05)]
# criteria: log2FC > 1 and q < 0.05 at at least one time point
respond_mock <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[,3] < p_cutoff & abs(x[,1]) > 1 )})
respond_mock_genes <- names(which(unlist(respond_mock) != 0))
not_respond_genes = setdiff(names(glm_fixef),respond_mock_genes)
flat_genes_profile = lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[,3] > p_cutoff)})
# flat genes, no response at almost no time point for EV and AvrRpt2
flat_genes = names(which(unlist(flat_genes_profile) > 19))
save(flat_genes,file = "flat_genes.Rdata")


AvrRpt2_mock_positive <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[1:10,3] < p_cutoff & x[1:10,1] > 1 )})
AvrRpt2_mock_negative <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[1:10,3] < p_cutoff & x[1:10,1] < -1 )})

AvrRpt2_mock_positive_genes1 <- names(which(unlist(AvrRpt2_mock_positive) != 0))
AvrRpt2_mock_negative_genes1 <- names(which(unlist(AvrRpt2_mock_negative) != 0))
# positive responded (induced) in AvrRpt2, no time point of suppression
AvrRpt2_mock_positive_genes <- setdiff(AvrRpt2_mock_positive_genes1,AvrRpt2_mock_negative_genes1)
# negative responded (suppressed) in AvrRpt2, no time point of induction 
AvrRpt2_mock_negative_genes <- setdiff(AvrRpt2_mock_negative_genes1,AvrRpt2_mock_positive_genes1)

# induced at the first three time points
AvrRpt2EV_mock_positive_first_three_hours <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[c(1,2,3,11,12,13),3] < p_cutoff & x[c(1,2,3,11,12,13),1] > 1 )})
AvrRpt2EV_mock_positive_first_three_hours_genes <- names(which(unlist(AvrRpt2EV_mock_positive_first_three_hours) != 0))
# suppressed at the first three time points
AvrRpt2EV_mock_negative_first_three_hours <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[c(1,2,3,11,12,13),3] < p_cutoff & x[c(1,2,3,11,12,13),1] < -1 )})
AvrRpt2EV_mock_negative_first_three_hours_genes <- names(which(unlist(AvrRpt2EV_mock_negative_first_three_hours) != 0))

# get those genes, induced or suppressed at the first three time points
AvrRpt2_mock_positive_genes_respond_before_3h <- intersect(AvrRpt2_mock_positive_genes,AvrRpt2EV_mock_positive_first_three_hours_genes)
AvrRpt2_mock_negative_genes_respond_before_3h <- intersect(AvrRpt2_mock_negative_genes,AvrRpt2EV_mock_negative_first_three_hours_genes)

# induced or suppressed at later time points,
# These are what we wanted
AvrRpt2_mock_positive_genes <- setdiff(AvrRpt2_mock_positive_genes,AvrRpt2EV_mock_positive_first_three_hours_genes)
AvrRpt2_mock_negative_genes <- setdiff(AvrRpt2_mock_negative_genes,AvrRpt2EV_mock_negative_first_three_hours_genes)

All.genes = names(glm_fixef)
save(AvrRpt2_mock_positive_genes,AvrRpt2_mock_positive_genes_respond_before_3h,All.genes,
     file = "AvrRpt2_genes.Rdata")
save(AvrRpt2_mock_negative_genes,AvrRpt2_mock_negative_genes_respond_before_3h,All.genes,
     file = "AvrRpt2_negative_genes.Rdata")



