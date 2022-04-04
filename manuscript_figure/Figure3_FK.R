#plot the synthetic profile after fixing 
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(gridExtra)
setwd("/Users/liux/Documents/MCM_paper/data/Figure3_data/")
load("profile_mat.Rdata")
load("AvrRpt2_genes.Rdata")
rownames(profiles_mat) = AvrRpt2_mock_positive_genes
colnames(profiles_mat) = c("early_4h","early_6h","early_9h","early_12h",
                           "early_16h","early_20h","early_24h","late_4h",
                           "late_6h","late_9h","late_12h","late_16h",
                           "late_20h","late_24h","res_4h","res_6h","res_9h",
                           "res_12h","res_16h","res_20h","res_24h")

#plot the synthetic profiles without residuals for set1
load("select.highquality.genes.Rdata")
select.profile_mat = profiles_mat[select.highquality.genes,c(1:14)]
colnames(select.profile_mat) = c("4","6","9","12","16","20","24",
                                 "4","6","9","12","16","20","24")

load("genes_to_fix.Rdata")
load("actual_fixed_mat.Rdata")
select.profile_mat.fixed = select.profile_mat
select.profile_mat.fixed[genes_to_be_fixed,] = profiles_mat_fixed[genes_to_be_fixed,1:14]

load("peak_info_MCM_spline.Rdata")
load("peak_info_fixed.Rdata")
peak_info_after_fix = peak_info[select.highquality.genes,]
peak_info_after_fix[genes_to_be_fixed,] = peak_info_fixed[genes_to_be_fixed,]
peak_level_ratio =  peak_info_after_fix[,2] / peak_info_after_fix[,4]
gene_group1 = names(which(peak_level_ratio > 4))
gene_group2 = names(which(peak_level_ratio < 4 & peak_level_ratio > 1))
gene_group3 = names(which(peak_level_ratio < 1))

gene_name_pt_group1 = names(sort(peak_info_after_fix[gene_group1,1]))
gene_name_pt_group2 = names(sort(peak_info_after_fix[gene_group2,1]))
gene_name_pt_group3 = names(sort(peak_info_after_fix[gene_group3,1]))
gene_name_pt_group4 = setdiff(select.highquality.genes, c(gene_name_pt_group1,gene_name_pt_group2,gene_name_pt_group3))

mat1 = select.profile_mat.fixed[gene_name_pt_group1,]
mat2 = select.profile_mat.fixed[gene_name_pt_group2,]
mat3 = select.profile_mat.fixed[gene_name_pt_group3,]
mat4 = select.profile_mat.fixed[gene_name_pt_group4,]

bot_anno = HeatmapAnnotation(foo = anno_text(colnames(select.profile_mat.fixed), rot = 0, 
                                       just = 'center',gp = gpar(fontsize = 11), location = unit(1,'mm')),
                       show_annotation_name = F)
left_anno_1 = rowAnnotation(group = rep('1',length(gene_name_pt_group1)), 
                            col = list(group = c('1' = 'red','2'='blue','3' = 'darkgoldenrod1','4' = 'purple')),
                            show_annotation_name = F)
left_anno_2 = rowAnnotation(group = rep('2',length(gene_name_pt_group2)), 
                            col = list(group = c('1' = 'red','2'='blue','3' = 'darkgoldenrod1','4' = 'purple')),
                            show_annotation_name = F)
left_anno_3 = rowAnnotation(group = rep('3',length(gene_name_pt_group3)), 
                            col = list(group = c('1' = 'red','2'='blue','3' = 'darkgoldenrod1','4' = 'purple')),
                            show_annotation_name = F)

left_anno_4 = rowAnnotation(group = rep('4',length(gene_name_pt_group4)), 
                            col = list(group = c('1' = 'red','2'='blue','3' = 'darkgoldenrod1','4' = 'purple')),
                            show_annotation_name = F)


col = colorRamp2(c(-2, 0, 1, 2, 6), c("cornflowerblue","white","orange",'orangered',"red4"))
lgd = Legend( col_fun = col, title = "log(FC)",legend_height = unit(30,"mm"),grid_width = unit(4,"mm"),
              title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 11),title_gap = unit(4,"mm"))
ht_1 = Heatmap(mat1 , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'ACP-specific',
               column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
               column_split = rep(c("first peak profiles","second peak profiles"),each=7),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=14), 
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               row_title = paste(length(gene_name_pt_group1), "genes"),
               use_raster = F,width = unit(90,'mm'), show_heatmap_legend = F,
               height = unit(length(gene_name_pt_group1)/16, 'mm'),
               left_annotation = left_anno_1)
ht_2 = Heatmap(mat2 , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'echoing genes',
               column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
               column_split = rep(c("first peak profiles","second peak profiles"),each=7),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=14), 
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               row_title = paste(length(gene_name_pt_group2), "genes"),
               use_raster = F,width = unit(90,'mm'), show_heatmap_legend = F,
               height = unit(length(gene_name_pt_group2)/16, 'mm'),
               left_annotation = left_anno_2 )
ht_3 = Heatmap(mat3 , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'NACP-preferential',
               column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
               column_split = rep(c("first peak profiles","second peak profiles"),each=7),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=14), 
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               row_title = paste(length(gene_name_pt_group3), "genes"),
               use_raster = F,width = unit(90,'mm'), show_heatmap_legend = F,
               height = unit(length(gene_name_pt_group3)/16, 'mm'),
               left_annotation = left_anno_3)
ht_4 = Heatmap(mat4 , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'NACP-specific',
               column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
               column_split = rep(c("first peak profiles","second peak profiles"),each=7),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=14), 
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               row_title = paste(length(gene_name_pt_group4), "genes"),
               use_raster = F,width = unit(90,'mm'), show_heatmap_legend = F,
               height = unit(length(gene_name_pt_group4)/16, 'mm'),
               bottom_annotation = bot_anno,
               left_annotation = left_anno_4)
hm_AvrRpt2 = grid.grabExpr(draw(ht_1 %v% ht_2 %v% ht_3 %v% ht_4, 
                                annotation_legend_list = lgd))

load('glm.fixef.RData')
WT_labels = c('flg22_genotypeJEPS:flg22_time0','flg22_genotypeJEPS:flg22_time1',
              'flg22_genotypeJEPS:flg22_time2','flg22_genotypeJEPS:flg22_time3',
              'flg22_genotypeJEPS:flg22_time5','flg22_genotypeJEPS:flg22_time9',
              'flg22_genotypeJEPS:flg22_time18')
fls2_labels = c('flg22_genotypefls2:flg22_time0','flg22_genotypefls2:flg22_time1',
                'flg22_genotypefls2:flg22_time2','flg22_genotypefls2:flg22_time3',
                'flg22_genotypefls2:flg22_time5','flg22_genotypefls2:flg22_time9',
                'flg22_genotypefls2:flg22_time18')
flg22_mat = c()
for (mygene in select.highquality.genes){
  if (mygene %in% names(glm_fixef)){
    WT_exp = glm_fixef[[mygene]][WT_labels,1]
    fls2_exp = glm_fixef[[mygene]][fls2_labels,1]
    flg22_mat = rbind(flg22_mat, log(2) * (WT_exp - fls2_exp))
  }else{
    flg22_mat = rbind(flg22_mat, rep(NA, 7))
  }
}
rownames(flg22_mat) = select.highquality.genes
colnames(flg22_mat) = c('0','1','2','3','5','9','18')

### normalize flg22_mat each row by making the vector length to 4.7
flg22_ori = flg22_mat

flg22_mat = t( apply(flg22_ori[,-c(2,7)], 1, function(x) {
  x.size = sqrt(sum(x^2))
  x/x.size * 4.7
}) )

flg22_mat1 = flg22_mat[gene_name_pt_group1,]
flg22_mat2 = flg22_mat[gene_name_pt_group2,]
flg22_mat3 = flg22_mat[gene_name_pt_group3,]
flg22_mat4 = flg22_mat[gene_name_pt_group4,]

col_flg22 = colorRamp2(c(-3.5,-2.7,-1.5, 0, 2,2.7,3.5,5), c("cornflowerblue",'lightskyblue1','white','white','white','orange','orangered',"red4"))

bot_anno = HeatmapAnnotation(foo = anno_text(colnames(flg22_mat1[,1:5]), rot = 0, 
                                             just = 'center',gp = gpar(fontsize = 11), location = unit(1,'mm')),
                             show_annotation_name = F)
flg22_ht1 = Heatmap(flg22_mat1[,1:5], cluster_rows = F, col = col_flg22,
              name = 'flg22_1', column_title = 'flg22 response',
              cluster_columns = FALSE,show_row_names = F, 
              column_title_gp = gpar(fontsize=14), column_names_gp = gpar(fontface = 'bold'),
              border = "gray20", column_names_rot = 0,
              na_col = "gray90", use_raster = F,
              width = unit(32,'mm'), height = unit(length(gene_name_pt_group1)/16, 'mm'),
              show_heatmap_legend = F)
flg22_ht2 = Heatmap(flg22_mat2[,1:5], cluster_rows = F, col = col_flg22,
                    name = 'flg22_2',
                    cluster_columns = FALSE,show_row_names = F, 
                    column_title_gp = gpar(fontsize=14), column_names_gp = gpar(fontface = 'bold'),
                    border = "gray20", column_names_rot = 0,
                    na_col = "gray90", use_raster = F,
                    width = unit(32,'mm'), height = unit(length(gene_name_pt_group2)/16, 'mm'),
                    show_heatmap_legend = F)
flg22_ht3 = Heatmap(flg22_mat3[,1:5], cluster_rows = F, col = col_flg22,
                    name = 'flg22_3',
                    cluster_columns = FALSE,show_row_names = F, 
                    column_title_gp = gpar(fontsize=14), column_names_gp = gpar(fontface = 'bold'),
                    border = "gray20", column_names_rot = 0,
                    na_col = "gray90", use_raster = F,
                    width = unit(32,'mm'), height = unit(length(gene_name_pt_group3)/16, 'mm'),
                    show_heatmap_legend = F)

flg22_ht4 = Heatmap(flg22_mat4[,1:5], cluster_rows = F, col = col_flg22,
                    name = 'flg22_4', show_column_names = F,
                    cluster_columns = FALSE, show_row_names = F, 
                    column_title_gp = gpar(fontsize=14),
                    border = "gray20", column_names_rot = 0,
                    na_col = "gray90", use_raster = F,
                    width = unit(32,'mm'), height = unit(length(gene_name_pt_group4)/16, 'mm'),
                    show_heatmap_legend = F, bottom_annotation = bot_anno)
lgd_flg22 = Legend( col_fun = col_flg22, title = "log(FC)_flg22",legend_height = unit(30,"mm"),grid_width = unit(4,"mm"),
              title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 11),title_gap = unit(4,"mm"))
hm_flg22 = grid.grabExpr(draw(flg22_ht1 %v% flg22_ht2 %v% flg22_ht3 %v% flg22_ht4, 
                                annotation_legend_list = lgd_flg22))

jpeg("/Users/liux/Documents/MCM_paper/manuscript_figure/Figure3_FK.jpeg",height = 150, width = 250, units = "mm",res = 300)
grid.arrange(hm_AvrRpt2,hm_flg22, layout_matrix = rbind(c(1,1,1,2,2),c(1,1,1,2,2)))
grid.text('time (hour)', x = unit(140,'mm'),y = unit(6,'mm'))
grid.text('A', x = unit(5,'mm'),y = unit(145,'mm'), gp= gpar(fontsize = 9,fontface = 'bold'))
grid.text('B', x = unit(160,'mm'),y = unit(145,'mm'),gp= gpar(fontsize = 9,fontface = 'bold'))
dev.off()
save(gene_name_pt_group1,gene_name_pt_group2,gene_name_pt_group3,gene_name_pt_group4,file = "gene_grouped_by_peak_ratio.Rdata")

