# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate the panels in supplemental figure 1  for the manuscript

# loading the libraries
library(ComplexHeatmap) # heatmap package
library(circlize)  # color package

# directory of MCM project, /MCM_paper/data
setwd("/Users/liux/Documents/MCM_paper/data")

#glm of AvrRpt2 ETI data (Kenichi's data)
load("glm_fixef_Ken_WT_with_pseudocounts.Rdata")
#labels of glm
labels_Rpt2 = c("mytreatAvrRpt2:mytime01h","mytreatAvrRpt2:mytime02h","mytreatAvrRpt2:mytime03h",
           "mytreatAvrRpt2:mytime04h","mytreatAvrRpt2:mytime06h","mytreatAvrRpt2:mytime09h",
           "mytreatAvrRpt2:mytime12h","mytreatAvrRpt2:mytime16h","mytreatAvrRpt2:mytime20h",
           "mytreatAvrRpt2:mytime24h")
labels_mock = c("mytreatmock:mytime01h","mytreatmock:mytime02h","mytreatmock:mytime03h",
           "mytreatmock:mytime04h","mytreatmock:mytime06h","mytreatmock:mytime09h",
           "mytreatmock:mytime12h","mytreatmock:mytime16h","mytreatmock:mytime20h",
           "mytreatmock:mytime24h")

#make data mat, it takes ~ 3 min
dm = c()
for (mygene in names(glm_fixef)){
  dm = rbind(dm, (glm_fixef[[mygene]][labels_Rpt2,1] - glm_fixef[[mygene]][labels_mock,1]))
}

rownames(dm) = names(glm_fixef)
colnames(dm) = c('1', '2', "3", "4", "6", "9", "12", "16", "20", "24")
my_dist = dist(dm,method = "euclidean")
hc = hclust(my_dist, method = "average")

#not consider some extreme values for color scales in the heatmap
min_value = quantile(as.vector(dm),probs = 0.001)
max_value = quantile(as.vector(dm),probs = 0.999)
#color scale
col = colorRamp2(c(min_value, 0, 0.5 *max_value, max_value), c("dodgerblue","white","coral1","darkred"))

jpeg("/Users/liux/Documents/MCM_paper/manuscript_figure/sup.fig1a.jpeg", width = 220, height = 180, units = "mm", res = 350)
lgd = Legend(col_fun = col, title = expression(log[2]~FC),legend_height = unit(40,"mm"),legend_width = unit(18,"mm"),
             title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 17),title_gap = unit(4,"mm"))
# takes ~ 5 min
ht1 = Heatmap(dm, col = col, cluster_rows = hc,name = 'logFC',
              cluster_columns = FALSE, column_names_rot = 0, column_names_centered = T,
              show_row_names = F,show_column_dend = FALSE, show_row_dend = T, 
              na_col = "#E0E0E0", row_title = paste(length(glm_fixef), "\n genes"),
              row_title_gp = gpar(fontsize = 24),row_title_rot = 0,
              row_dend_width = unit(20, "mm"), use_raster = F, heatmap_width = unit(180,"mm"),
              column_names_gp = gpar(fontsize = 20),
              show_heatmap_legend = F)
draw(ht1, ht_gap = unit(c(1), "cm"), annotation_legend_list = lgd)
grid.text('time (hour)', x = unit(205,'mm'),y = unit(6,'mm'),gp = gpar(fontsize = 17))
dev.off()

#grab the panel and further merge into a single figure, takes ~ 4 min
sup1a_plot = grid.grabExpr(draw(ht1, ht_gap = unit(c(1), "cm"), annotation_legend_list = lgd))
#sup 1a done


#sup 1b, signal-to-noise histogram, starting from 3039 AvrRpt2-induced genes
setwd("/Users/liux/Documents/MCM_paper/data")
load("/Users/liux/Documents/MCM_paper/data/glm_fixef_Ken_WT_with_pseudocounts.Rdata")
load("/Users/liux/Documents/MCM_paper/data/AvrRpt2_genes.Rdata")
Rpt2_labels = c("mytreatAvrRpt2:mytime03h","mytreatAvrRpt2:mytime04h","mytreatAvrRpt2:mytime06h",
                "mytreatAvrRpt2:mytime09h","mytreatAvrRpt2:mytime12h","mytreatAvrRpt2:mytime16h",
                "mytreatAvrRpt2:mytime20h","mytreatAvrRpt2:mytime24h")
signal_to_noise.all = c()
for (mygene in AvrRpt2_mock_positive_genes){
  mean.est = glm_fixef[[mygene]][Rpt2_labels,"Estimate"]
  ste = glm_fixef[[mygene]][Rpt2_labels,"Std.Error"]
  sig.noi = (max(mean.est) - min(mean.est) )/ max(ste)
  signal_to_noise.all = c(signal_to_noise.all, sig.noi)
}
names(signal_to_noise.all) = AvrRpt2_mock_positive_genes

jpeg("/Users/liux/Documents/MCM_paper/manuscript_figure/sup.fig1.b.jpeg",width = 150, height = 140,unit = "mm",res = 400)
par(mar = c(5,5,2,2))
hist(signal_to_noise.all,xlab = "signal-to-noise ratio",main = "",ylab = "Number of genes",breaks = 50,
     col = 'gray80',border = "gray20",cex.lab = 1.8, cex.axis = 1.5)
abline(v = 6.5,lty = 5,lwd = 2, col = "red")
dev.off()
save(signal_to_noise.all,file = "signal_to_noise.all.Rdata")
#sup 1b done


#select genes with high signal-noise ratio for further analysis
setwd("/Users/liux/Documents/MCM_paper/data")
load("AvrRpt2_genes.Rdata")
load("signal_to_noise.all.Rdata")
select.genes = names(which(signal_to_noise.all > 6.5))
save(select.genes,file = "select.genes.Rdata")


#sup 1c start
#load 2435 upregulated genes after removing low signal-noise fit genes
load("select.genes.Rdata")
labels_Rpt2 = c("mytreatAvrRpt2:mytime03h","mytreatAvrRpt2:mytime04h","mytreatAvrRpt2:mytime06h",
                "mytreatAvrRpt2:mytime09h","mytreatAvrRpt2:mytime12h","mytreatAvrRpt2:mytime16h",
                "mytreatAvrRpt2:mytime20h","mytreatAvrRpt2:mytime24h")
labels_mock = c("mytreatmock:mytime03h","mytreatmock:mytime04h","mytreatmock:mytime06h",
                "mytreatmock:mytime09h","mytreatmock:mytime12h","mytreatmock:mytime16h",
                "mytreatmock:mytime20h","mytreatmock:mytime24h")
#glm of AvrRpt2 ETI data (Kenichi's data)
load("glm_fixef_Ken_WT_with_pseudocounts.Rdata")
dm = c()
for (mygene in select.genes){
  dm = rbind(dm, glm_fixef[[mygene]][labels_Rpt2,1] - glm_fixef[[mygene]][labels_mock,1])
}
rownames(dm) = select.genes
colnames(dm) = c("3", "4", "6", "9", "12", "16", "20", "24")
my_dist = dist(dm,method = "euclidean")
hc= hclust(my_dist, method = "average")
min_value = quantile(as.vector(dm),probs = 0.001)
max_value = quantile(as.vector(dm),probs = 0.999)
col = colorRamp2(c(min_value, 0, 0.5 *max_value, max_value), c("dodgerblue","white","coral1","darkred"))
jpeg("/Users/liux/Documents/MCM_paper/manuscript_figure/sup.fig1.c.jpeg", width = 220, height = 180, units = "mm", res = 350)
lgd = Legend(col_fun = col, title = expression(log[2]~FC),legend_height = unit(40,"mm"),legend_width = unit(8,"mm"),
             title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 15),title_gap = unit(4,"mm"))
ht1 = Heatmap(dm, col = col, cluster_rows = hc,name = "log2 scale",
              cluster_columns = FALSE, column_names_rot = 0,column_names_centered = T,
              show_row_names = F,show_column_dend = FALSE, show_row_dend = T, 
              na_col = "#E0E0E0", row_title = paste(length(select.genes), "\n genes"),
              row_title_gp = gpar(fontsize = 20),row_title_rot = 0,
              row_dend_width = unit(20, "mm"), use_raster = F, heatmap_width = unit(180,"mm"),
              column_names_gp = gpar(fontsize = 19),
              show_heatmap_legend = F)
draw(ht1, ht_gap = unit(c(1), "cm"), annotation_legend_list = lgd)
grid.text('time (hour)', x = unit(205,'mm'),y = unit(6,'mm'),gp = gpar(fontsize = 17))
dev.off()
sup_fig1c = grid.grabExpr(draw(ht1, ht_gap = unit(c(1), "cm"), annotation_legend_list = lgd))
#sup 1c done


#***************#
#sup 1d PCC between fit of input 1 and 2
setwd("/Users/liux/Documents/MCM_paper/data")

#synthetic profiles, input1
load("profile_mat.Rdata")
#2435 genes
load("select.genes.Rdata")
profile_set1 = profiles_mat
#synthetic profiles, input2
load("profile_mat_mid.Rdata")
profile_set2 = profiles_mat
#3039 genes
load("AvrRpt2_genes.Rdata")
rownames(profile_set1) = AvrRpt2_mock_positive_genes
rownames(profile_set2) = AvrRpt2_mock_positive_genes
select.profile.set1 = profile_set1[select.genes,]
select.profile.set2 = profile_set2[select.genes,]
cor_list = c()
for (i in 1:length(select.genes)){
  cor = cor(profile_set1[i,1:14],profile_set2[i,1:14])
  cor_list = c(cor_list, cor)
}
names(cor_list) = select.genes 
select.highquality.genes = names(which(cor_list > 0.9))
jpeg("/Users/liux/Documents/MCM_paper/manuscript_figure/sup.fig1.d.jpeg",
     width = 220, height = 150, unit = "mm",res = 400)
par(mar = c(5,5,2,4))
hist(cor_list,xlab = "PCC between synthetic profiles fit with input set1 and set2",
     main = "2435 genes",ylab = "Number of genes",breaks = 50,
     col = 'gray80',border = "gray20",cex.lab = 1.8, cex.axis = 1.8,cex.main = 2,
     ylim = c(0,760), xlim = c(-0.35, 1.2))
abline(v = 0.9,lty = 5,lwd = 2,col = "red")
text('1939 genes', x = 1.06, y = 750, cex = 1.5)
dev.off()
save(select.highquality.genes,file = "select.highquality.genes.Rdata")
#end




