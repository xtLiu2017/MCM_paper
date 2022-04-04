# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate the histogram showing GLM vs MCM fitting
# Previously, we put it in the main figures but now we make it in supplemental figure 4

# directory of MCM project, /MCM_paper/data
setwd("~/Documents/MCM_paper/data")

# glm fitting of AvrRpt2 ETI data
load("glm_fixef_Ken_WT_with_pseudocounts.Rdata")
AR_labels = c('mytreatAvrRpt2:mytime03h','mytreatAvrRpt2:mytime04h','mytreatAvrRpt2:mytime06h',
              'mytreatAvrRpt2:mytime09h','mytreatAvrRpt2:mytime12h','mytreatAvrRpt2:mytime16h',
              'mytreatAvrRpt2:mytime20h','mytreatAvrRpt2:mytime24h')
MK_labels = c('mytreatmock:mytime03h','mytreatmock:mytime04h','mytreatmock:mytime06h',
              'mytreatmock:mytime09h','mytreatmock:mytime12h','mytreatmock:mytime16h',
              'mytreatmock:mytime20h','mytreatmock:mytime24h')
# 2435 genes
load("select.genes.Rdata")
# synthetic profiles, input 1
load('profile_mat.Rdata')
#3039  genes
load("AvrRpt2_genes.Rdata")
rownames(profiles_mat) = AvrRpt2_mock_positive_genes
# PCC between GLM and MCM
PCC_vec = c()
for (mygene in select.genes){
  glm_mk = glm_fixef[[mygene]][MK_labels, 1]
  glm_AR = glm_fixef[[mygene]][AR_labels, 1]
  glm_diff = glm_AR - glm_mk
  nlm = profiles_mat[mygene,1:7] + profiles_mat[mygene, 8:14]
  myPCC = cor(nlm, glm_diff[2:8])
  PCC_vec = c(PCC_vec, myPCC)
}
names(PCC_vec) = select.genes
#1939 genes
load("select.highquality.genes.Rdata")
jpeg('~/Documents/MCM_paper/manuscript_figure/sup.fig4.jpeg',width = 126, height = 100, unit = 'mm', res = 700)
par(mar = c(5,5,3.6,2))
hist(PCC_vec, breaks = 60, xlab = "PCC", ylab = "Number of genes", 
     main = "Pearson correlation coefficients (PCCs) \n between fit of GLMs and MCMs",
     col.lab = "black", axes = F, col = 'gray80',border = "gray20",cex.lab = 1.4, 
     cex.axis = 1.5, cex.main = 1.2)
axis(side = 1, col = "black",col.axis = "black",cex.axis = 1.2, lwd = 2)
axis(side = 2, col = "black", col.axis = 'black',cex.axis = 1.2, lwd = 2)
abline(v = median(PCC_vec),lty = 5,lwd = 2, col = "red")
text(x = 0.28, y = 460, labels = "Median PCC = 0.89",cex = 1.4, col = "red")
dev.off()



