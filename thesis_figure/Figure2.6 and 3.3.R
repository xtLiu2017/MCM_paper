# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate figure 2.6 and 3.3 for my PhD thesis
# The plots are the synthetic heatmaps. 
# Most genes were fit with classical MCMs, some were fit with altered MCMs

# loading the libraries
library(ComplexHeatmap)   # heatmap package
library(circlize)         # supporting package
library(RColorBrewer)     # supporting package, for colors
library(grid)             # supporting package, for merging multiple heatmaps in a page
library(gridExtra)        # same above

# directory of MCM project, /MCM_paper/data
setwd("~/Documents/MCM_paper/data/")

# load the data, 3039 genes MCM profiles
load("profile_mat.Rdata")
load("AvrRpt2_genes.Rdata")
rownames(profiles_mat) = AvrRpt2_mock_positive_genes

# plot the synthetic profiles without residuals for set1
# 1939 high quality genes (upregulated and modeled genes)
load("select.highquality.genes.Rdata")
select.profile_mat = profiles_mat[select.highquality.genes,c(1:14)]
colnames(select.profile_mat) = c("4","6","9","12","16","20","24",
                                 "4","6","9","12","16","20","24")

# genes with altered fitting
load("genes_to_fix.Rdata")
# the MCM profiles of altered fitting
load("actual_fixed_mat.Rdata")
select.profile_mat.fixed = select.profile_mat
select.profile_mat.fixed[genes_to_be_fixed,] = profiles_mat_fixed[genes_to_be_fixed,1:14]

#genes fit with classical MCMs
gene.not.fixed = setdiff(select.highquality.genes, genes_to_be_fixed)
#synthetic profile matrix (two peaks)
select.not.fixed.mat = select.profile_mat.fixed[gene.not.fixed,]
select.fixed.mat = select.profile_mat.fixed[genes_to_be_fixed,]
#residual matrix
residual.not.fixed.mat = profiles_mat[gene.not.fixed,15:21]
residual.fixed.mat = profiles_mat_fixed[genes_to_be_fixed,15:21]
#mean estimate matrix
mean.glm.not.fixed.mat = select.not.fixed.mat[,1:7] + select.not.fixed.mat[,8:14] + residual.not.fixed.mat
mean.glm.fixed.mat = select.fixed.mat[,1:7] + select.fixed.mat[,8:14] + residual.fixed.mat
#model fitting matrix
model.not.fixed.mat = select.not.fixed.mat[,1:7] + select.not.fixed.mat[,8:14]
model.fixed.mat = select.fixed.mat[,1:7] + select.fixed.mat[,8:14]

#clustering based on the synthetic profiles
dist_mat.not.fixed = dist(select.not.fixed.mat)
mhclust.not.fixed = hclust(dist_mat.not.fixed,method = "average")
dist_mat.fixed = dist(select.fixed.mat)
mhclust.fixed = hclust(dist_mat.fixed,method = "average")

#color scale
col = colorRamp2(c(-2, 0, 1, 2, 6), c("cornflowerblue","white","orange",'orangered',"red4"))
#annotations, for columns to show the time points
bot.anno = HeatmapAnnotation(foo = anno_text(colnames(mean.glm.not.fixed.mat), rot = 0, 
                                       just = 'center',gp = gpar(fontsize = 11), location = unit(1,'mm')),
                       show_annotation_name = F)
#mean estimate heatmaps
ht.mean.estimate.not.fixed = Heatmap(mean.glm.not.fixed.mat, col = col, cluster_rows = mhclust.not.fixed,
              row_dend_gp = gpar(col = "gray50"),
              cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 20),
              show_row_names = F, show_column_names = F, show_column_dend = FALSE, show_row_dend = F, 
              column_title_gp = gpar(fontsize=18), column_title = 'GLM mean estimate',
              border = "gray20", row_title_rot = 0, 
              row_title = paste(dim(mean.glm.not.fixed.mat)[1], "\n genes"),
              row_dend_width = unit(15, "mm"), use_raster = F,
              width = unit(60,'mm'), height = unit(130, 'mm'),
              show_heatmap_legend = F)
ht.mean.estimate.fixed = Heatmap(mean.glm.fixed.mat, col = col, cluster_rows = mhclust.fixed,
                                     row_dend_gp = gpar(col = "gray50"),
                                     cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 20),
                                     show_row_names = F, show_column_names = F, show_column_dend = FALSE, show_row_dend = F, 
                                     column_title_gp = gpar(fontsize=14), 
                                     border = "gray20", row_title_rot = 0, 
                                     row_title = paste(dim(mean.glm.fixed.mat)[1], "\n genes"),
                                     row_dend_width = unit(15, "mm"), use_raster = F,
                                     width = unit(60,'mm'), height = unit(4, 'mm'),
                                     show_heatmap_legend = F, bottom_annotation = bot.anno)
# grab this heatmap plotting, used for merging plots
mean_estimate_hm = grid.grabExpr(draw(ht.mean.estimate.not.fixed %v% ht.mean.estimate.fixed)) 


#model fitting heatmaps
ht.model.not.fixed = Heatmap(model.not.fixed.mat, col = col, cluster_rows = mhclust.not.fixed,
                                     row_dend_gp = gpar(col = "gray50"),
                                     cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 20),
                                     show_row_names = F, show_column_names = F, show_column_dend = FALSE, show_row_dend = F, 
                                     column_title_gp = gpar(fontsize=18), column_title = 'MCM fitting',
                                     border = "gray20", row_title_rot = 0, 
                                     row_dend_width = unit(15, "mm"), use_raster = F,
                                     width = unit(60,'mm'), height = unit(130, 'mm'),
                                     show_heatmap_legend = F)
ht.model.fixed = Heatmap(model.fixed.mat, col = col, cluster_rows = mhclust.fixed,
                             row_dend_gp = gpar(col = "gray50"),
                             cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
                             show_row_names = F, show_column_names = F, show_column_dend = FALSE, show_row_dend = F, 
                             column_title_gp = gpar(fontsize=14), 
                             border = "gray20", row_title_rot = 0, 
                             row_dend_width = unit(15, "mm"), use_raster = F,
                             width = unit(60,'mm'), height = unit(4, 'mm'),
                             show_heatmap_legend = F,bottom_annotation = bot.anno)
# grab this heatmap plotting, used for merging plots
model_hm = grid.grabExpr(draw(ht.model.not.fixed %v% ht.model.fixed)) 

#residual heatmaps
ht.resi.not.fixed = Heatmap(residual.not.fixed.mat, col = col, cluster_rows = mhclust.not.fixed,
                             row_dend_gp = gpar(col = "gray50"),
                             cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 20),
                             show_row_names = F, show_column_names = F, show_column_dend = FALSE, show_row_dend = F, 
                             column_title_gp = gpar(fontsize=18), 
                             border = "gray20", row_title_rot = 0, column_title = 'Discrepancy',
                             row_dend_width = unit(15, "mm"), use_raster = F,
                             width = unit(60,'mm'), height = unit(130, 'mm'),
                             show_heatmap_legend = F)
ht.resi.fixed = Heatmap(residual.fixed.mat, col = col, cluster_rows = mhclust.fixed,
                         row_dend_gp = gpar(col = "gray50"),
                         cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 20),
                         show_row_names = F, show_column_names = F, show_column_dend = FALSE, show_row_dend = F, 
                         column_title_gp = gpar(fontsize=18), 
                         border = "gray20", row_title_rot = 0, 
                         row_dend_width = unit(15, "mm"), use_raster = F,
                         width = unit(60,'mm'), height = unit(4, 'mm'),
                         show_heatmap_legend = F, bottom_annotation = bot.anno)
# grab this heatmap plotting, used for merging plots
residual_hm = grid.grabExpr(draw(ht.resi.not.fixed %v% ht.resi.fixed)) 

#synthetic heatmaps
bot_anno_synthetic = HeatmapAnnotation(foo = anno_text(colnames(select.not.fixed.mat), rot = 0, 
                                       just = 'center',gp = gpar(fontsize = 11), location = unit(1,'mm')),
                       show_annotation_name = F)

ht.synthetic.not.fixed = Heatmap(select.not.fixed.mat, col = col, cluster_rows = mhclust.not.fixed,
                            row_dend_gp = gpar(col = "gray50"),
                            cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
                            show_row_names = F, show_column_names = F, show_column_dend = FALSE, 
                            show_row_dend = T, row_dend_side = 'right',
                            column_title_gp = gpar(fontsize=14), 
                            border = "gray20", row_title_rot = 0, 
                            column_split = rep(c("first peak profiles","second peak profiles"),each=7),
                            row_title = paste(dim(mean.glm.not.fixed.mat)[1], "\n genes"),
                            row_dend_width = unit(15, "mm"), use_raster = F,
                            width = unit(95,'mm'), height = unit(130, 'mm'),
                            show_heatmap_legend = F)
ht.synthetic.fixed = Heatmap(select.fixed.mat, col = col, cluster_rows = mhclust.fixed,
                        row_dend_gp = gpar(col = "gray50"),
                        cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
                        show_row_names = F, show_column_names = F, show_column_dend = FALSE, 
                        show_row_dend = T, row_dend_side = 'right',
                        row_title = paste(dim(mean.glm.fixed.mat)[1], "\n genes"),
                        column_title_gp = gpar(fontsize=14), 
                        border = "gray20", row_title_rot = 0, 
                        column_split = rep(c("first peak profiles","second peak profiles"),each=7),
              
                        row_dend_width = unit(15, "mm"), use_raster = F,
                        width = unit(95,'mm'), height = unit(4, 'mm'),
                        show_heatmap_legend = F, bottom_annotation = bot_anno_synthetic)
# legend for the synthetic profile
lgd = Legend( col_fun = col, title = "log(FC)",legend_height = unit(30,"mm"),grid_width = unit(4,"mm"),
              title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 11),title_gap = unit(4,"mm"))
# grab this heatmap plotting, used for merging plots
synthetic_hm = grid.grabExpr(draw(ht.synthetic.not.fixed %v% ht.synthetic.fixed, 
                                  annotation_legend_list = lgd)) 


# layout matrix indicating how to organize the heatmaps in a page
lay = rbind(c(1,1,1,1,1,2,2,2,2,3,3,3,3),c(1,1,1,1,1,2,2,2,2,3,3,3,3))
# for below, change the output directory as you need
jpeg('~/Desktop/paper_figures/thesis/Figure2.6.jpeg',width = 340, height = 157, unit = 'mm',res = 500)
grid.arrange(mean_estimate_hm,model_hm,residual_hm,layout_matrix = lay)
dev.off()
jpeg('~/Desktop/paper_figures/thesis/Figure3.3.jpeg',width = 160, height = 180, unit = 'mm',res = 500)
draw(ht.synthetic.not.fixed %v% ht.synthetic.fixed, 
     annotation_legend_list = lgd)
dev.off()





