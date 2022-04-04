#plot glm fitting of the genes
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

gene_list = c('AT1G74710', 'AT1G19180', 'AT2G40750', 'AT1G19250', 'AT1G80590', 'AT1G64280')
gene_labels = c('SID2','JAZ1','WRKY54','FMO1','WRKY66','NPR1')
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



