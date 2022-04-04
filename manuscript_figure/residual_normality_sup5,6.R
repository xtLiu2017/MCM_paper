# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate figure panels for distribution of the residuals
# and explain how we applied the altered MCM to the genes (peak + shoulders)
# I name it supplemental figure 5

#plot residual distribution
#The codes below generate all residual distributions

# directory of MCM project, /MCM_paper/data
setwd("~/Documents/MCM_paper/data")

#libraries
library("dplyr")
library("ggpubr")

# MCM synthetic profiles input set 1
load("profile_mat.Rdata")
# 3039 genes
load("AvrRpt2_genes.Rdata")
rownames(profiles_mat) = AvrRpt2_mock_positive_genes
colnames(profiles_mat) = c("early_4h","early_6h","early_9h","early_12h",
                           "early_16h","early_20h","early_24h","late_4h",
                           "late_6h","late_9h","late_12h","late_16h",
                           "late_20h","late_24h","res_4h","res_6h","res_9h",
                           "res_12h","res_16h","res_20h","res_24h")
# 1939 high quality genes
load("select.highquality.genes.Rdata")
#normality based on shapiro test for all time points
normality_4h = shapiro.test(profiles_mat[select.highquality.genes,15])
normality_6h = shapiro.test(profiles_mat[select.highquality.genes,16])
normality_9h = shapiro.test(profiles_mat[select.highquality.genes,17])
normality_12h = shapiro.test(profiles_mat[select.highquality.genes,18])
normality_16h = shapiro.test(profiles_mat[select.highquality.genes,19])
normality_20h = shapiro.test(profiles_mat[select.highquality.genes,20])
normality_24h = shapiro.test(profiles_mat[select.highquality.genes,21])

#starting the plots here
jpeg("~/Documents/MCM_paper/manuscript_figure/sup.5.jpeg", width = 350, height = 350, units = "mm",res = 400)
par(mfrow = c(5,4), mar = c(3,5,3,1), oma = c(3,8,3,3))
hist(profiles_mat[select.highquality.genes,15],breaks = seq(-2.2,2.2,0.02), main = "4h",ylab = "Number of genes",
     xlab = "",xlim = c(-2.2,2.2),ylim  = c(0,200),cex.lab = 2,cex.axis = 2, cex.main = 3)
mtext(side = 2, line = 11, text = "threshold = Inf \n 0 genes fixed", at = 100, padj = 1,cex = 1.5)
text(x = -1, y = 150, label = paste("W = ", round(normality_4h$statistic,digits = 3)), cex = 1.8)

hist(profiles_mat[select.highquality.genes,16],breaks = seq(-2.2,2.2,0.02), main = "6h", cex.main = 3,
     xlab = "",xlim = c(-2.2,2.2),ylim  = c(0,200),cex.lab = 2, yaxt = 'n', ylab = NULL,cex.axis = 2)
text(x = -1, y = 150, label = paste("W = ", round(normality_6h$statistic,digits = 3)), cex = 1.8)

hist(profiles_mat[select.highquality.genes,17],breaks = seq(-2.2,2.2,0.02), main = "9h", cex.main = 3,
     xlab = "",xlim = c(-2.2,2.2),ylim  = c(0,200),cex.lab = 2,yaxt = 'n', ylab = NULL,cex.axis = 2)
text(x = -1, y = 150, label = paste("W = ", round(normality_9h$statistic,digits = 3)), cex = 1.8)
hist(profiles_mat[select.highquality.genes,18],breaks = seq(-2.2,2.2,0.02), main = "12h",cex.main = 3,
     xlab = "",xlim = c(-2.2,2.2),ylim  = c(0,200),cex.lab = 2,yaxt = 'n', ylab = NULL,cex.axis = 2)
text(x = -1, y = 150, label = paste("W = ", round(normality_12h$statistic,digits = 3)), cex = 1.8)


# Altered MCM profiles (232 genes; in the end, not all 232 genes were fit with altered MCMs)
load("actual_fixed_mat.Rdata")
actual_fixed_genes = rownames(actual_fixed_mat)

# Classical MCM profiles with 1939 genes
select.profile_mat = profiles_mat[select.highquality.genes,]
all.residuals = c()
fixed.residuals = c()
#thresholds for residuals at 6 hours
all.thresholds = c(-1, -0.55, -0.35, -0.25)
genes_fixed_vec = c()
for (mycutoff1 in all.thresholds){
  #genes selected based on the threshold for residuals at 6 hours
  genes_aa = names(which(profiles_mat[select.highquality.genes,16] < mycutoff1))
  #genes should have positive residual at 9 hours
  genes_bb = names(which(profiles_mat[select.highquality.genes,17] > 0))
  # these genes with altered fitting
  genes_to_be_fixed = intersect(intersect(genes_aa,genes_bb),actual_fixed_genes)
  # number of genes with altered fitting
  no.of.genes = length(genes_to_be_fixed)
  genes_fixed_vec = c(genes_fixed_vec,no.of.genes )
  print(no.of.genes)
  # change the profile for the genes with altered MCM fitting
  select.profile_mat.fixed_percutoff = select.profile_mat
  select.profile_mat.fixed_percutoff[genes_to_be_fixed,] = profiles_mat_fixed[genes_to_be_fixed,]
  # residual after change
  all.residuals = rbind(all.residuals, select.profile_mat.fixed_percutoff[,15])
  all.residuals = rbind(all.residuals, select.profile_mat.fixed_percutoff[,16])
  all.residuals = rbind(all.residuals, select.profile_mat.fixed_percutoff[,17])
  all.residuals = rbind(all.residuals, select.profile_mat.fixed_percutoff[,18])
}
# This is the optimal threshold
mycutoff1 = -0.55
genes_aa = names(which(profiles_mat[select.highquality.genes,16] < mycutoff1))
genes_bb = names(which(profiles_mat[select.highquality.genes,17] > 0))
genes_to_be_fixed = intersect(intersect(genes_aa,genes_bb),actual_fixed_genes)
no.of.genes = length(genes_to_be_fixed)
fixed.residuals = cbind(profiles_mat[genes_to_be_fixed,c(15,16,17,18)],
                        profiles_mat_fixed[genes_to_be_fixed,c(15,16,17,18)])
# Now make the histograms for the residuals distributions at all time points for all thresholds
j = 1
for (i in 1:16){
  normality_i = round(shapiro.test(all.residuals[i,])$statistic,digits = 3)
  if (i %% 4 == 1){
    hist(all.residuals[i,],breaks = 100,main = "",xlim  = c(-2.2,2.2),xlab = '', 
         ylim = c(0,200),ylab = "Number of genes",cex.lab = 2, cex.axis = 2)
    mtext(side = 2, line = 11, text = paste("threshold = ",all.thresholds[j], "\n", genes_fixed_vec[j], "genes fixed"), 
                                            at = 100, padj = 1,cex = 1.5)
    j = j + 1
  }else{
    hist(all.residuals[i,],breaks = 100,main = "",xlim  = c(-2.2,2.2),xlab = '',
         ylim = c(0,200),ylab = "",cex.lab = 2, yaxt = 'n', cex.axis = 2)
  }
  text(paste('W = ',normality_i),x = -1,y = 150,cex = 1.8)
}
dev.off()

#plot residual distribution done
#The codes above generate supplemental figure 5 for all residual distributions
# --------- end ---------



#plot some example genes before and after revised fitting
#this is for supplemental figure 6
#make plots to see the original fitting and altered fitting 

# directory of MCM project, /MCM_paper/data
setwd("~/Documents/MCM_paper/data")

# 3039 genes
load("AvrRpt2_genes.Rdata")
# genes with altered fitting (optimal threshold obtained above)
load("genes_to_fix.Rdata")
# a function script 
source("~/Documents/MCM_paper/MCM_function.R")
# read conut data
load("/Users/liux/Documents/Ken_analysis/data/col_count_data.Rdata")
all.genes <- rownames(col_count_data)
# prepare explanatory variables
treatment <- c()
time <- c()
for (mylib in colnames(col_count_data)){
  treatment <- c(treatment, strsplit(mylib,split = "_")[[1]][4])
  time <- c(time, strsplit(mylib,split = "_")[[1]][5])
}
time <- as.numeric(substr(time,start=1,stop=2))
after_3_indexes <- which(time >= 3)
time_after3 <- time[after_3_indexes]
treatment_after3 <- treatment[after_3_indexes]
col_offset_after3 <- col_offset[after_3_indexes]
offset_mock_after3 <- col_offset_after3[treatment_after3 == "mock"]
offset_Rpt2_after3 <- col_offset_after3[treatment_after3 == "AvrRpt2"]

#k parameters for input set 1
k_vec1 = c(1.8, 0.5487, 1.1733, 1.1253, 0.6492, 0.6853, 0.6957, 0.6987, 0.4949,
           0.5065, 0.513, 0.5174, 0.2133)

compart_vec <- c(1,3,5,7,9,11)

# MCM parameters
load("merged_result.Rdata")
parameter_best_before = parameter_best
rownames(parameter_best_before) = AvrRpt2_mock_positive_genes
# altered MCM parameters
load("merged_result_16h_part.Rdata")
rm(parameter_best)

jpeg("~/Documents/MCM_paper/manuscript_figure/sup.6.jpeg",width = 5.8, height = 4.8, units = "cm", res = 700)
i = 1; j = 1
par(mfrow = c(3,2),mar = c(1.8,2.2,0.4,0.8))
for (mygene in genes_to_be_fixed[c(1,2,4)]){
  if (i %% 2 == 1){
    my_ylab = substitute(paste("log(", italic('Pto'), " AvrRpt2/mock)", sep =""))
  }else{
    my_ylab = " "
  }
  if (j > 4){
    my_xlab = "time (hour)"
  }else{
    my_xlab = " "
  }
  if (i == 1){
    my_main = 'before revised fitting'
  }else if (i == 2){
    my_main = 'after revised fitting'
  }else{
    my_main = ''
  }
 
  mygene_index = which(AvrRpt2_mock_positive_genes == mygene)
  # read count data
  y_after3 <- as.numeric(col_count_data[mygene,after_3_indexes])
  y_Rpt2 <- y_after3[treatment_after3 == "AvrRpt2"]
  y_mock <- y_after3[treatment_after3 == "mock"]
  # x -- time
  x <- as.numeric(time_after3[treatment_after3 == "mock"])
  # y -- log read count - offset
  log_Rpt2 <- log(y_Rpt2) - offset_Rpt2_after3
  log_mock <- log(y_mock) - offset_mock_after3
  
  compart = parameter_best_before[mygene, 12]
  # first input compartment, second input compartment
  compart1_indx <- compartment_function(compart)[1]
  compart2_indx <- compartment_function(compart)[2]
  compart1 <- compart_vec[compart1_indx]
  compart2 <- compart_vec[compart2_indx]
  # two functions based on two inputs
  func1_str <- func_string(compart1)
  func1 <- eval(parse(text =func1_str))
  func2_str <- func_string(compart2)
  func2 <- eval(parse(text =func2_str))
  # classical MCM parameter
  parameter_Rpt2 <- parameter_best_before[mygene,1:3]
  a1 = parameter_Rpt2[1]
  a2 = parameter_Rpt2[2]
  k1 = parameter_Rpt2[3]
  # mock parameters
  parameter_mock <- parameter_best_before[mygene,4:9]
  # mock fitting, polynomial
  pre_mock <- model_output(t=seq(3,24,0.01), p = c(1,1,1, parameter_mock), k_vec=k_vec1,mock = 1)
  # MCM fitting
  pre_Rpt2 <- model_output(t=seq(3,24,0.01), p = c(parameter_Rpt2,parameter_mock), k_vec=k_vec1,mock = 0)
  # fitting samples
  mock_to_subtract <- model_output(t=x, p = c(1,1,1,parameter_mock), k_vec=k_vec1,mock = 1)
  mock_to_subtract_glm = model_output(t=c(3,4,6,9,12,16,20,24), p = c(1,1,1,parameter_mock), k_vec=k_vec1,mock = 1)
  
  # range of y axis
  ymax_other <- max(log_Rpt2 - mock_to_subtract) + 0.5
  ymin_other <- min(log_Rpt2 - mock_to_subtract) - 0.5
  plot(x,log_Rpt2 - mock_to_subtract,xlab = "",ylab = ' ',main = my_main,
       ylim=c(ymin_other,ymax_other),xlim = c(0, 25), frame.plot = F, cex = 0.5,cex.lab = 0.35, 
       axes = F, lwd = 0.4, cex.main = 0.5)
  text(mygene,x = 13,y = 0.92 * ymax_other + 0.08 * ymin_other,cex = 0.4)
  axis(side = 1, lwd = 0.5, cex.axis = 0.5,tck = -0.08,mgp=c(2,0.1,0))
  axis(side = 2, lwd = 0.5, cex.axis = 0.5, tck = -0.06,mgp=c(2,0.2,0))
  title(ylab = my_ylab, cex.lab = 0.38, line = 1)
  title(xlab = my_xlab, cex.lab = 0.43, line = 0.7 )
  lines(seq(3,24,0.01),pre_Rpt2 - pre_mock,lwd = 0.5)
  i = i + 1; j = j + 1
  #-----
  if (i %% 2 == 1){
    my_ylab = substitute(paste("log(", italic('Pto'), " AvrRpt2/mock)", sep =""))
  }else{
    my_ylab = " "
  }
  if (j > 4){
    my_xlab = "time (hour)"
  }else{
    my_xlab = " "
  }
  if (i == 1){
    my_main = 'before revised fitting'
  }else if (i == 2){
    my_main = 'after revised fitting'
  }else{
    my_main = ''
  }
  # plot the altered MCM fitting
  compart = parameter_best_16h_part[mygene, 12]
  # again, two compartments, first input, second input
  compart1_indx <- compartment_function(compart)[1]
  compart2_indx <- compartment_function(compart)[2]
  compart1 <- compart_vec[compart1_indx]
  compart2 <- compart_vec[compart2_indx]
  func1_str <- func_string(compart1)
  func1 <- eval(parse(text =func1_str))
  func2_str <- func_string(compart2)
  func2 <- eval(parse(text =func2_str))
  
  # altered MCM parameters
  parameter_Rpt2.after <- parameter_best_16h_part[mygene,1:3]
  a1.after = parameter_Rpt2.after[1]
  a2.after = parameter_Rpt2.after[2]
  k1.after = parameter_Rpt2.after[3]
  # altered MCM fitting
  pre_Rpt2.after <- model_output(t=seq(3,24,0.01), p = c(parameter_Rpt2.after,parameter_mock), k_vec=k_vec1,mock = 0)
  # some data points not used (> 16 hours)
  colors = c()
  for (ii in 1:length(x)){
    if (x[ii] > 16){
      colors = c(colors,"gray85")
    }else{
      colors = c(colors,"black") 
    }
  }
  
  plot(x,log_Rpt2 - mock_to_subtract,xlab = "",ylab = "",ylim=c(ymin_other,ymax_other),yaxt="n",
       col = colors,xlim = c(0, 25), frame.plot = F, cex = 0.5,cex.lab = 0.35, axes = F, 
       lwd = 0.4, main = my_main, cex.main = 0.5)
  text(mygene,x = 13,y = 0.92 * ymax_other + 0.08 * ymin_other,cex = 0.4)
  axis(side = 1, lwd = 0.5, cex.axis = 0.5,tck = -0.08,mgp=c(2,0.1,0))
  axis(side = 2, lwd = 0.5, cex.axis = 0.5, tck = -0.06,mgp=c(2,0.2,0))
  title(ylab = my_ylab, cex.lab = 0.38, line = 1)
  title(xlab = my_xlab, cex.lab = 0.43, line = 0.7 )
  lines(seq(3,24,0.01),pre_Rpt2.after - pre_mock,lwd = 0.5)
  i = i + 1; j = j + 1
}
dev.off()

#the code above plot fixed genes before and after revised fitting
#for supplemental figure 6
#done
