# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate the panels of example gene fitting
# Previously, we put it in the main figures but now we make them in supplemental figure 2

# directory of MCM project, /MCM_paper/data
setwd("~/Documents/MCM_paper/data/")

# example genes, I only used the first three genes in the supplemental to show goodness-of-fit
# The other two genes are PROPER genes, used in anohter supplemental figure
example1 = "AT5G62570"
example2 = 'AT3G56400' 
example3 = 'AT3G48080'
example4 = 'AT5G64890'
example5 = 'AT5G64905'
example_labels = c("CBP60a", "WRKY70", "EDS1b",'PROPEP2','PROPEP3')
names(example_labels) = c(example1,example2,example3,example4, example5)

# This files contains all functions related to MCM fitting
# It is usedful here to show model fitting
source("/Users/liux/Documents/MCM_paper/MCM_function.R")
# read count data of wild type plant AvrRpt2 ETI
load("col_count_data.Rdata")
all.genes <- rownames(col_count_data)
#explanatory variables below
treatment <- c()
time <- c()
for (mylib in colnames(col_count_data)){
  treatment <- c(treatment, strsplit(mylib,split = "_")[[1]][4])
  time <- c(time, strsplit(mylib,split = "_")[[1]][5])
}
time <- as.numeric(substr(time,start=1,stop=2))
#remove before 3 hours data
after_3_indexes <- which(time >= 3)
time_after3 <- time[after_3_indexes]
treatment_after3 <- treatment[after_3_indexes]
col_offset_after3 <- col_offset[after_3_indexes]
offset_mock_after3 <- col_offset_after3[treatment_after3 == "mock"]
offset_Rpt2_after3 <- col_offset_after3[treatment_after3 == "AvrRpt2"]

# k values for input set 1
k_vec = c(1.8, 0.5487, 1.1733, 1.1253, 0.6492, 0.6853, 0.6957, 0.6987, 0.4949,
          0.5065, 0.513, 0.5174, 0.2133)
compart_vec <- c(1,3,5,7,9,11)
peak_time <- c(1,2,3,4.5,6,7.5,9,11,13,15,17,21) * 0.95
# parameters of MCM fit
load("merged_result.Rdata")
rownames(parameter_best) = rownames(col_count_data)

for (mygene in c(example1, example2,example3,example4, example5)){
  mylabel <- example_labels[mygene]
  #change the output directory accordingly
  jpeg(paste('/Users/liux/Documents/MCM_paper/manuscript_figure/sup.fig2.',mylabel,'.jpeg',sep = ""),
       width = 126, height = 100, unit = 'mm', res = 350)
  par(mar = c(5,5,2.5,3.5))
  y_after3 <- as.numeric(col_count_data[mygene,after_3_indexes])
  y_Rpt2 <- y_after3[treatment_after3 == "AvrRpt2"]
  y_mock <- y_after3[treatment_after3 == "mock"]
  #count data above
  
  #time
  x <- as.numeric(time_after3[treatment_after3 == "mock"])
  
  #log count - offset
  log_Rpt2 <- log(y_Rpt2) - offset_Rpt2_after3
  log_mock <- log(y_mock) - offset_mock_after3
  
  #compartment
  compart = parameter_best[mygene, 12]
  compart1_indx <- compartment_function(compart)[1]
  compart2_indx <- compartment_function(compart)[2]
  
  compart1 <- compart_vec[compart1_indx]
  compart2 <- compart_vec[compart2_indx]
  # first peak function, determined by compartment 1
  func1_str <- func_string(compart1)
  func1 <- eval(parse(text =func1_str))
  # second peak function, determined by compartment 2
  func2_str <- func_string(compart2)
  func2 <- eval(parse(text =func2_str))
  
  # parameters, a1,a2,k and mock parameters (polynomial function)
  parameter_Rpt2 <- parameter_best[mygene,1:3]
  a1 = parameter_Rpt2[1]
  a2 = parameter_Rpt2[2]
  k1 = parameter_Rpt2[3]
  parameter_mock <- parameter_best[mygene,4:9]
  
  # mock fitting with polynomial
  pre_mock <- model_output(t=seq(3,24,0.01), p = c(1,1,1, parameter_mock), k_vec=k_vec,mock = 1)
  # AvrRpt2 fitting with MCM
  pre_Rpt2 <- model_output(t=seq(3,24,0.01), p = c(parameter_Rpt2,parameter_mock), k_vec=k_vec,mock = 0)
  
  # mock fitting samples
  mock_to_subtract <- model_output(t=x, p = c(1,1,1,parameter_mock), k_vec=k_vec,mock = 1)
  
  #range of the y axis
  diff = max(log_Rpt2 - mock_to_subtract) - min(log_Rpt2 - mock_to_subtract)
  ymax_other <- max(log_Rpt2 - mock_to_subtract) + diff * 0.15
  ymin_other <- min(log_Rpt2 - mock_to_subtract) - diff * 0.15
  
  plot(x,log_Rpt2 - mock_to_subtract,ylab = "transcript level (log)", xlab = "time (hour)",
       ylim=c(ymin_other,ymax_other), frame.plot = F, xlim = c(0,25), cex = 1.5,
       col.lab = "black", axes = F,cex.lab = 1.4, lwd = 1.6)
  axis(side = 1, lwd = 2, col = "black",col.axis = "black",cex.axis = 1.2)
  axis(side = 2, lwd = 2, col = "black",col.axis = "black",cex.axis = 1.2)
  text(mylabel,x = 12,y = ymax_other,cex = 1.4,font = 4)
  lines(seq(3,24,0.01),pre_Rpt2 - pre_mock, lwd = 2.2)
  dev.off()
}

#-------------------------------------------------
# spline fitting of some immunity-related genes, 
# This is to show the double peak patterns are prevalent
# I wanted to make it supplemental figure 3
# This part is independent from the part above, so you can start running the code here

# directory of MCM project, /MCM_paper/data
setwd("~/Documents/MCM_paper/data")

# read count data
load("col_count_data.Rdata")
all.genes <- rownames(col_count_data)
# explanatory variables
treatment <- c()
time <- c()
for (mylib in colnames(col_count_data)){
  treatment <- c(treatment, strsplit(mylib,split = "_")[[1]][4])
  time <- c(time, strsplit(mylib,split = "_")[[1]][5])
}
time <- as.numeric(substr(time,start=1,stop=2))
offset_mock <- col_offset[treatment == "mock"]
offset_Rpt2 <- col_offset[treatment == "AvrRpt2"]
# mock parameters
# This mock parameters here are exactly the same as the merge_result.Rdata
load("mock_start.vals.Rdata")
rownames(mock_coef) = all.genes
gene_list = c('AT1G11310', 'AT1G19180', 
              'AT2G40750', 'AT1G19250', 'AT1G80590', 'AT1G09970')
gene_labels = c('MLO2','JAZ1','WRKY54','FMO1','WRKY66','RLK7')
names(gene_labels) = gene_list
glm_labels = c(rep(c('x_time1','x_time3','x_time4','x_time6','x_time9','x_time12','x_time16','x_time24'), time = 3),
rep(c('x_time2','x_time20'), each = 3))
library(MASS)
jpeg("~/Documents/MCM_paper/manuscript_figure/sup.fig3.jpeg",width = 5.3, height = 4.6, units = "cm", res = 600)
par(mfrow = c(3,2),mar = c(1.8,1.7,0.4,0.4))
i = 1
j = 1
for (mygene in gene_list){
  if (mygene %in% all.genes){
    # determine the x label and y label
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
    mygene_label = gene_labels[mygene]
    # count data
    y <- as.numeric(col_count_data[mygene,])
    y_Rpt2 <- y[treatment == "AvrRpt2"]
    y_mock <- y[treatment == "mock"]
    # time 
    myx_time <- time[treatment == "mock"]
    x_time <- as.character(time[treatment == "mock"])
    #glm fitting
    glm_mock = tryCatch({glm.nb(y_mock ~ -1 + x_time + offset(offset_mock),
                                link=log)},error = function(err){
                                  return(2)
                                })
    mock_offset = glm_mock$coefficients[glm_labels]
    y_Rpt2_subtract_mock_offset = log(y_Rpt2) - offset_Rpt2 - mock_offset
    # spline fitting
    mysp.fun = splinefun(sqrt(myx_time), y_Rpt2_subtract_mock_offset)
    mysp.val = mysp.fun(sqrt(seq(1,24,0.01)))
    # y axis range
    min_lim = min(y_Rpt2_subtract_mock_offset) - max(y_Rpt2_subtract_mock_offset) * 0.1
    max_lim = max(y_Rpt2_subtract_mock_offset) * 1.4
    
    plot(myx_time, y_Rpt2_subtract_mock_offset, xlab = "time (hour)",ylab = my_ylab,
         ylim = c(min_lim, max_lim), xlim = c(0, 25), frame.plot = F, cex = 0.5, 
         cex.lab = 0.35, axes = F, ann =F, lwd = 0.4)
    axis(side = 1, lwd = 0.5, cex.axis = 0.5,tck = -0.08,mgp=c(2,0.1,0))
    axis(side = 2, lwd = 0.5, cex.axis = 0.5, tck = -0.06,mgp=c(2,0.2,0))
    title(ylab = my_ylab, cex.lab = 0.38, line = 1)
    title(xlab = my_xlab, cex.lab = 0.6, line = 0.7 )
    lines(seq(1,24,0.01),mysp.val, lwd = 0.5)
    text(mygene_label,x = 14,y = max_lim * 0.93,cex = 0.5,font = 4)
    i = i + 1
    j = j + 1
  }
}
dev.off()



#** ------------------------------- **#
# The codes below are used for my defense slides
# Just ignore this for the manuscript
# plot glm fitting of the genes
setwd("~/Documents/MCM_paper/data")
load("col_count_data.Rdata")
load("glm_fixef_Ken_WT_with_pseudocounts.Rdata")
all.genes <- rownames(col_count_data)
treatment <- c()
time <- c()
for (mylib in colnames(col_count_data)){
  treatment <- c(treatment, strsplit(mylib,split = "_")[[1]][4])
  time <- c(time, strsplit(mylib,split = "_")[[1]][5])
}
time <- as.numeric(substr(time,start=1,stop=2))
offset_mock <- col_offset[treatment == "mock"]
offset_Rpt2 <- col_offset[treatment == "AvrRpt2"]

gene_list = c('AT1G11310', 'AT1G19180', 
              'AT2G40750', 'AT1G19250', 'AT1G80590', 'AT1G09970')
gene_labels = c('MLO2','JAZ1','WRKY54','FMO1','WRKY66','RLK7')
names(gene_labels) = gene_list
AvrRpt2_labels = c('mytreatAvrRpt2:mytime01h','mytreatAvrRpt2:mytime02h','mytreatAvrRpt2:mytime03h',
                   'mytreatAvrRpt2:mytime04h','mytreatAvrRpt2:mytime06h','mytreatAvrRpt2:mytime09h',
                   'mytreatAvrRpt2:mytime12h','mytreatAvrRpt2:mytime16h','mytreatAvrRpt2:mytime20h',
                   'mytreatAvrRpt2:mytime24h')
mock_labels = c('mytreatmock:mytime01h','mytreatmock:mytime02h','mytreatmock:mytime03h',
              'mytreatmock:mytime04h','mytreatmock:mytime06h','mytreatmock:mytime09h',
              'mytreatmock:mytime12h','mytreatmock:mytime16h','mytreatmock:mytime20h',
              'mytreatmock:mytime24h')
glm_time = c(1,2,3,4,6,9,12,16,20,24)
library(MASS)
jpeg("figures/ppt_glm.jpeg",width = 5.3, height = 4.6, units = "cm", res = 600)
par(mfrow = c(3,2),mar = c(1.8,1.7,0.4,0.4))
i = 1
j = 1
for (mygene in gene_list){
  if (mygene %in% all.genes){
    if (i %% 2 == 1){
      my_ylab = 'GLM estimate'
    }else{
      my_ylab = " "
    }
    if (j > 4){
      my_xlab = "time (hour)"
    }else{
      my_xlab = " "
    }
    mygene_label = gene_labels[mygene]
    y <- as.numeric(col_count_data[mygene,])
    y_Rpt2 <- y[treatment == "AvrRpt2"]
    y_mock <- y[treatment == "mock"]
    myx_time <- time[treatment == "mock"]
    x_time <- as.character(time[treatment == "mock"])
    
   
    y_Rpt2_subtract_offset = log(y_Rpt2) - offset_Rpt2 
    y_mock_subtract_offset = log(y_mock) - offset_mock
   
    min_lim = min(y_Rpt2_subtract_offset) - max(y_Rpt2_subtract_offset) * 0.1
    max_lim = max(y_Rpt2_subtract_offset) * 1.3
    
    plot(myx_time, y_Rpt2_subtract_offset, xlab = "time (hour)",ylab = my_ylab,
         ylim = c(min_lim, max_lim), xlim = c(0, 25), frame.plot = F, cex = 0.5, 
         cex.lab = 0.35, axes = F, ann =F, lwd = 0.4, col = "dark red")
    axis(side = 1, lwd = 0.5, cex.axis = 0.5,tck = -0.08, mgp=c(2,0.1,0))
    axis(side = 2, lwd = 0.5, cex.axis = 0.5, tck = -0.06, mgp=c(2,0.2,0))
    title(ylab = my_ylab, cex.lab = 0.5, line = 1)
    title(xlab = my_xlab, cex.lab = 0.6, line = 0.7 )
    
    text(mygene_label,x = 10,y = max_lim * 0.93,cex = 0.5,font = 4)
    
    AvrRpt2_glm = glm_fixef[[mygene]][AvrRpt2_labels,1] * log(2)
    lines(glm_time, AvrRpt2_glm, col = "dark red")
    points(myx_time, y_mock_subtract_offset,cex = 0.5, cex.lab = 0.35,lwd = 0.4,col = "orange" )
    mock_glm = glm_fixef[[mygene]][mock_labels,1] * log(2)
    lines(glm_time, mock_glm, col = "orange")
    if (i == 1 | i == 2){
      legend("topright", legend=c(expression(italic("Pto") * AvrRpt2), "mock"),
             col=c("dark red", "orange"), cex=0.38,lty=1, bty = 'n')
    }
   
    i = i + 1
    j = j + 1
  }
}
dev.off()



