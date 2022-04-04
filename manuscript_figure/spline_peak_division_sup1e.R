# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate the supplemental figure 1E 
# It shows an empirical way of decomposing two peak profiles based on spline method

# directory of MCM project, /MCM_paper/data
setwd("~/Documents/MCM_paper/data")

#packages
library(ComplexHeatmap)   #heatmap
library(circlize)         #colors

# 3039 significantly induced genes
load("AvrRpt2_genes.Rdata")
# MCM profile, input1
load("profile_mat.Rdata")
profile_mat_set1 = profiles_mat
rownames(profile_mat_set1) = AvrRpt2_mock_positive_genes
# MCM profiles, input2
load("profile_mat_mid.Rdata")
profile_mat_set2 = profiles_mat
rownames(profile_mat_set2) = AvrRpt2_mock_positive_genes
time_index = c()
for (i in c(4,6,9,12,16,20,24)){
  time_index = c(time_index, which(seq(3,24,0.01) ==i ))  }   

# 1939 high quality genes
load("select.highquality.genes.Rdata")
# spline profiles
load("spline_complete_synthetic_profiles.Rdata")
# we take the intersection because some genes may not be fit well with spline model
genes_final = intersect(rownames(spline_profiles),select.highquality.genes)   

# make the profiles together, input1, input2, spline
select.profile_mat = cbind(profile_mat_set1[genes_final,c(1:14)],
                           profile_mat_set2[genes_final,c(1:14)],
                           spline_profiles[genes_final,c(time_index,time_index + 2101)])
colnames(select.profile_mat) = c("peak_1_4h","peak_1_6h","peak_1_9h","peak_1_12h","peak_1_16h",
                                 "peak_1_20h","peak_1_24h","peak_2_4h","peak_2_6h","peak_2_9h","peak_2_12h",
                                 "peak_2_16h", "peak_2_20h","peak_2_24h","peak_1_4h","peak_1_6h","peak_1_9h","peak_1_12h","peak_1_16h",
                                 "peak_1_20h","peak_1_24h","peak_2_4h","peak_2_6h","peak_2_9h","peak_2_12h",
                                 "peak_2_16h", "peak_2_20h","peak_2_24h","peak_1_4h","peak_1_6h","peak_1_9h","peak_1_12h","peak_1_16h",
                                 "peak_1_20h","peak_1_24h","peak_2_4h","peak_2_6h","peak_2_9h","peak_2_12h",
                                 "peak_2_16h", "peak_2_20h","peak_2_24h")
# colors for legend
col = colorRamp2(c(-2, 0, 2, 6), c("cornflowerblue","white","coral1","darkred"))
dist_mat = dist(select.profile_mat)
# hierarchical clustering
mhclust = hclust(dist_mat,method = "average")
# legend
lgd = Legend(col_fun = col, title = "log(FC)",legend_height = unit(30,"mm"),legend_width = unit(8,"mm"),
             title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14),title_gap = unit(4,"mm"))                     
ht1 = Heatmap(select.profile_mat, col = col, cluster_rows = mhclust, name = "log(FC)",
              cluster_columns = FALSE, column_names_gp = gpar(fontsize = 7),
              row_dend_gp = gpar(col = "gray50"),border = "gray20",
              show_row_names = F,show_column_dend = FALSE, show_row_dend = T, 
              na_col = "#E0E0E0", row_title = paste(dim(select.profile_mat)[1],"genes"),
              column_gap = unit(2,"mm"), row_title_rot = 0,
              column_title_gp = gpar(fill = c("#FED9A6", '#B3CDE3', '#CCEBC5'),col = "gray20",fontsize=14),
              row_dend_width = unit(25, "mm"), use_raster = F,heatmap_height = unit(150,"mm"),
              heatmap_width = unit(190,"mm"), column_names_rot = 45,
              column_split = rep(c("MCM (set1)","MCM (set2)","Spline"),each=14),
              show_heatmap_legend = F)
jpeg("~/Documents/MCM_paper/manuscript_figure/sup.fig1e.jpeg",height = 160, width = 220, units = "mm",res = 400)
draw(ht1, annotation_legend_list = lgd)
dev.off()


# The bottom part is not used, it just compute the PCC between these profiles (input1, input2, spline)
cor_all.1 = c()
cor_all.2 = c()
for (mygene in rownames(select.profile_mat)){
  profile_set1_my = select.profile_mat[mygene,1:14]
  profile_set2_my = select.profile_mat[mygene,15:28]
  profile_spline_my = select.profile_mat[mygene,29:42]
  mycor1 = cor(profile_set1_my,profile_spline_my)
  mycor2 = cor(profile_set2_my,profile_spline_my)
  cor_all.1 = c(cor_all.1,mycor1)
  cor_all.2 = c(cor_all.2,mycor2)
}

hist(cor_all.1,xlab = "correlation between spline profiles and set1 profiles",main="",cex.lab = 1.2)
hist(cor_all.2,xlab = "correlation between spline profiles and set2 profiles",main="",cex.lab=1.2)


