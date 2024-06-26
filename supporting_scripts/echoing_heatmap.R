# Author: Xiaotong Liu (liux4215@umn.edu)
# plot the synthetic profiles, including genes with classical MCM fitting and altered MCM fitting

# packages
library(ComplexHeatmap)   # heatmap
library(circlize)         # color
library(RColorBrewer)     # more colors

# directory where the data are saved
setwd("~/Documents/MCM_paper/data")
# classical MCM profiles, input 1
load("profile_mat.Rdata")
# 3039 upregulated genes
load("AvrRpt2_genes.Rdata")
rownames(profiles_mat) = AvrRpt2_mock_positive_genes
colnames(profiles_mat) = c("early_4h","early_6h","early_9h","early_12h",
                           "early_16h","early_20h","early_24h","late_4h",
                           "late_6h","late_9h","late_12h","late_16h",
                           "late_20h","late_24h","res_4h","res_6h","res_9h",
                           "res_12h","res_16h","res_20h","res_24h")

# 1939 high quality genes
load("select.highquality.genes.Rdata")
select.profile_mat = profiles_mat[select.highquality.genes,c(1:14)]
colnames(select.profile_mat) = c("4","6","9","12","16","20","24",
                                 "4","6","9","12","16","20","24")

# genes fit with altered MCMs
load("genes_to_fix.Rdata")
# profiles of genes fit with altered MCMs 
load("actual_fixed_mat.Rdata")
select.profile_mat.fixed = select.profile_mat
select.profile_mat.fixed[genes_to_be_fixed,] = profiles_mat_fixed[genes_to_be_fixed,1:14]

# distance matrix
dist_mat = dist(select.profile_mat.fixed)
#hierarchical clustering
mhclust = hclust(dist_mat,method = "average")
# time point annotation
ha = HeatmapAnnotation(foo = anno_text(colnames(select.profile_mat.fixed), rot = 0, 
                                       just = 'center',gp = gpar(fontsize = 11), location = unit(1,'mm')),
                       show_annotation_name = F)

# colors
col = colorRamp2(c(-2, 0, 1, 2, 6), c("cornflowerblue","white","orange",'orangered',"red4"))
# legend
lgd = Legend( col_fun = col, title = "log(FC)",legend_height = unit(30,"mm"),grid_width = unit(4,"mm"),
              title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 11),title_gap = unit(4,"mm"))
# heatmap
ht1 = Heatmap(select.profile_mat.fixed, col = col, cluster_rows = mhclust,
              row_dend_gp = gpar(col = "gray50"),
              cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
              column_split = rep(c("first peak profiles","second peak profiles"),each=7),
              show_row_names = F, show_column_names = F, show_column_dend = FALSE, show_row_dend = T, 
              column_title_gp = gpar(fontsize=14), 
              border = "gray20", row_title_rot = 0, 
              na_col = "#E0E0E0", row_title = paste(dim(select.profile_mat.fixed)[1], "\n genes"),
              row_dend_width = unit(15, "mm"), use_raster = F,
              width = unit(90,'mm'), height = unit(130, 'mm'),
              show_heatmap_legend = F, bottom_annotation = ha)

#change the output directory as needed
jpeg("figures/fig3.a.jpeg",height = 150, width = 180, units = "mm",res = 300)
draw(ht1, annotation_legend_list = lgd)
grid.text('time (hour)', x = unit(156,'mm'),y = unit(5,'mm'))
dev.off()