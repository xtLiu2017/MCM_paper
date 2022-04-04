#plot the synthetic profile after fixing 
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
setwd("/Users/liux/Documents/MCM_paper/data")
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
dist_mat = dist(select.profile_mat)
mhclust = hclust(dist_mat,method = "average")
colnames(select.profile_mat) = c("4","6","9","12","16","20","24",
                                 "4","6","9","12","16","20","24")

load("genes_to_fix.Rdata")
load("actual_fixed_mat.Rdata")
select.profile_mat.fixed = select.profile_mat
select.profile_mat.fixed[genes_to_be_fixed,] = profiles_mat_fixed[genes_to_be_fixed,1:14]
dist_mat = dist(select.profile_mat.fixed)
mhclust = hclust(dist_mat,method = "average")
mytree = cutree(mhclust,k = 5)
ha = HeatmapAnnotation(foo = anno_text(colnames(select.profile_mat.fixed), rot = 0, 
                                       just = 'center',gp = gpar(fontsize = 11), location = unit(1,'mm')),
                       show_annotation_name = F)


col = colorRamp2(c(-2, 0, 1, 2, 6), c("cornflowerblue","white","orange",'orangered',"red4"))
col_flg22 = colorRamp2(c(-3,-2,0,1.5,4,7), c("dodgerblue2","cornflowerblue",'white','orange','orangered',"red4"))
lgd = Legend( col_fun = col, title = "log(FC)",legend_height = unit(30,"mm"),grid_width = unit(4,"mm"),
              title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 11),title_gap = unit(4,"mm"))
ht1 = Heatmap(select.profile_mat.fixed, col = col, cluster_rows = mhclust,
              row_dend_gp = gpar(col = "gray50"), name = 'echoing',
              cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
              column_split = rep(c("first peak profiles","second peak profiles"),each=7),
              show_row_names = F, show_column_names = F, show_column_dend = FALSE, show_row_dend = T, 
              column_title_gp = gpar(fontsize=14), 
              border = "gray20", row_title_rot = 0, 
              na_col = "#E0E0E0", row_title = paste(dim(select.profile_mat.fixed)[1], "\n genes"),
              row_dend_width = unit(15, "mm"), use_raster = F,
              width = unit(90,'mm'), height = unit(130, 'mm'),
              show_heatmap_legend = T, bottom_annotation = ha)
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
ht2 = Heatmap(flg22_mat[,c(1,2,4,5,6)], cluster_rows = mhclust, col = col_flg22,
              row_dend_gp = gpar(col = "gray50"), name = 'flg22',
              cluster_columns = FALSE,  row_title_gp = gpar(fontsize = 14),
              show_row_names = F, show_column_dend = FALSE, show_row_dend = T, 
              column_title_gp = gpar(fontsize=14), column_names_gp = gpar(fontface = 'bold'),
              border = "gray20", row_title_rot = 0, column_names_rot = 0,
              na_col = "gray80", row_title = paste(dim(select.profile_mat.fixed)[1], "\n genes"),
              row_dend_width = unit(15, "mm"), use_raster = F,column_title = 'flg22 response',
              width = unit(60,'mm'), height = unit(130, 'mm'),
              show_heatmap_legend = T)

load("glm_fixef_Ken_WT_with_pseudocounts.Rdata")
EV_labels = c('mytreatEV:mytime04h','mytreatEV:mytime06h','mytreatEV:mytime09h','mytreatEV:mytime12h',
              'mytreatEV:mytime16h','mytreatEV:mytime20h','mytreatEV:mytime24h')
mock_labels = c('mytreatmock:mytime04h','mytreatmock:mytime06h','mytreatmock:mytime09h','mytreatmock:mytime12h',
                'mytreatmock:mytime16h','mytreatmock:mytime20h','mytreatmock:mytime24h')
EV_mat = c()
for (mygene in select.highquality.genes){
  if (mygene %in% names(glm_fixef)){
    EV_exp = glm_fixef[[mygene]][EV_labels,1]
    mock_exp = glm_fixef[[mygene]][mock_labels,1]
    EV_mat = rbind(EV_mat, log(2) * (EV_exp - mock_exp))
  }else{
    EV_mat = rbind(EV_mat, rep(NA, 7))
  }
}
rownames(EV_mat) = select.highquality.genes
colnames(EV_mat) = c('4','6','9','12','16','20','24')
ht3 = Heatmap(EV_mat, cluster_rows = mhclust, col = col,
              row_dend_gp = gpar(col = "gray50"), name = 'Pto',
              cluster_columns = FALSE,  row_title_gp = gpar(fontsize = 14),
              show_row_names = F, show_column_dend = FALSE, show_row_dend = T, 
              column_title_gp = gpar(fontsize=14), 
              border = "gray20", row_title_rot = 0, column_names_rot = 0,
              na_col = "gray80", 
              row_dend_width = unit(15, "mm"), use_raster = F,column_title = 'Pto response',
              width = unit(60,'mm'), height = unit(130, 'mm'),
              show_heatmap_legend = T)


jpeg("figures/echoing+flg22.jpeg",height = 150, width = 280, units = "mm",res = 300)
draw(ht1 + ht2 + ht3)
grid.text('time (hour)', x = unit(265,'mm'),y = unit(6,'mm'))
dev.off()

library(org.At.tair.db)
x <- org.At.tairGO
mapped_genes <- mappedkeys(x)
xx <- as.list(org.At.tairGO2TAIR)
xx['GO:0001666']
xx['GO:0071456']
length(xx['GO:0001666'][[1]])
length(xx['GO:0071456'][[1]])
intersect(xx['GO:0001666'][[1]], sp_genes)
intersect(xx['GO:0071456'][[1]], sp_genes)

setwd("/Users/liux/Documents/MCM_paper/data")
load("peak_info_MCM_spline.Rdata")
load("select.highquality.genes.Rdata")
load("genes_to_fix.Rdata")
load("peak_info_fixed.Rdata")
colnames(peak_info) = c('pt1_i1','pl1_i1','pt2_i1','pl2_i1')
#look at overlap of genes with GO term: hypoxia response
sliding_window_func = function(sorted_peak_level_ratio, genes_GO, window){
  No_overlap_vec = c()
  peak_ratio = c()
  for (mystart in 1:(length(sort_peak_level_ratio) - window)){
    genes_window = names(sort_peak_level_ratio[mystart: (mystart + window - 1)])
    overlap_genes = intersect(genes_window, genes_GO)
    No_overlap_vec = c(No_overlap_vec, length(overlap_genes))
    peak_ratio = c(peak_ratio, sort_peak_level_ratio[mystart + window/2])
  }
  return(rbind(peak_ratio,No_overlap_vec))
}
peak_info_1939 = peak_info[select.highquality.genes,]
peak_info_1939[genes_to_be_fixed,] = peak_info_fixed[genes_to_be_fixed,]
peak_level_ratio = peak_info_1939[,'pl1_i1'] / peak_info_1939[,'pl2_i1']
sort_peak_level_ratio = sort(peak_level_ratio, na.last = NA)
window1 = 30
No.overlap_30 = sliding_window_func(sorted_peak_level_ratio = sort_peak_level_ratio, 
                                    genes_GO  = xx['GO:0071456'][[1]], window = window1)
df_hypoxia = data.frame(index = 1:dim(No.overlap_30)[2], ratio = No.overlap_30[1,], No = No.overlap_30[2,])
myat = c()
for (ii in 1/c(0.25,0.5,1,1.25,2.5)){
  myat = c(myat, which.min(abs(ii - No.overlap_30[1,])))
}
GO_gene_no. = length(intersect(xx['GO:0071456'][[1]],names(sort_peak_level_ratio)))
MCM_gene_no. = length(sort_peak_level_ratio)
lower_0.95 = qhyper(p=0.025, m = GO_gene_no., n=MCM_gene_no. - GO_gene_no., k = window1, lower.tail = TRUE, log.p = FALSE)
upper_0.95 = qhyper(p=0.975, m = GO_gene_no., n=MCM_gene_no. - GO_gene_no., k = window1, lower.tail = TRUE, log.p = FALSE)
mycol <- rgb(142,229,238, max = 255, alpha = 125)

fig_sliding_hypoxia = ggplot(df_hypoxia, aes(x = index, y = No)) + geom_line() + theme_bw() + xlab('peak 1 level / peak 2 level') +
  ylab('No of genes in GO: \n cellular response tp hypoxia') + 
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        plot.margin = unit(c(15, 15, 6, 10), 'mm'), axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), panel.grid = element_blank()) + 
  scale_x_continuous(breaks = myat, labels = 1/c(0.25,0.5,1,1.25,2.5)) + 
  geom_hline(yintercept=window1 * GO_gene_no. / MCM_gene_no. , linetype = "dashed", col = "blue") + 
  annotate("rect", xmin = -50, ymin = lower_0.95, xmax = 1600, ymax = upper_0.95, alpha=0.5, fill=mycol) +
  annotate("text", x = 500, y = 4, label = 'expected number with \n 95% confidence interval',
           size = 3.5, col = "blue") + 
  annotate('text',x = 700, y = 8,  label = 'window size = 30 genes', size = 5)


genes_a = names(which(sort_peak_level_ratio > 4))
genes_b = names(which(sort_peak_level_ratio < 4 & sort_peak_level_ratio > 1))
genes_c = names(which(sort_peak_level_ratio < 1))
ht_a = Heatmap(select.profile_mat.fixed[genes_a,], col = col, cluster_rows = T,
              row_dend_gp = gpar(col = "gray50"), name = 'first peak preferential',
              cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
              column_split = rep(c("first peak profiles","second peak profiles"),each=7),
              show_row_names = F, show_column_names = F, show_column_dend = FALSE, show_row_dend = T, 
              column_title_gp = gpar(fontsize=14), 
              border = "gray20", row_title_rot = 0, 
              na_col = "#E0E0E0", row_title = paste(length(genes_a), "\n genes"),
              row_dend_width = unit(15, "mm"), use_raster = F,
              width = unit(90,'mm'), height = unit(length(genes_a)/12, 'mm'),
              show_heatmap_legend = F)
ht_b = Heatmap(select.profile_mat.fixed[genes_b,], col = col, cluster_rows = T,
               row_dend_gp = gpar(col = "gray50"), name = 'echoing genes',
               cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
               column_split = rep(c("first peak profiles","second peak profiles"),each=7),
               show_row_names = F, show_column_names = F, show_column_dend = FALSE, show_row_dend = T, 
               column_title_gp = gpar(fontsize=14), 
               border = "gray20", row_title_rot = 0, 
               na_col = "#E0E0E0", row_title = paste(length(genes_b), "\n genes"),
               row_dend_width = unit(15, "mm"), use_raster = F,
               width = unit(90,'mm'), height = unit(length(genes_b)/12, 'mm'),
               show_heatmap_legend = F)
ht_c = Heatmap(select.profile_mat.fixed[genes_c,], col = col, cluster_rows = T,
               row_dend_gp = gpar(col = "gray50"), name = 'second peak preferential',
               cluster_columns = FALSE, column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 14),
               column_split = rep(c("first peak profiles","second peak profiles"),each=7),
               show_row_names = F, show_column_names = F, show_column_dend = FALSE, show_row_dend = T, 
               column_title_gp = gpar(fontsize=14), 
               border = "gray20", row_title_rot = 0, 
               na_col = "#E0E0E0", row_title = paste(length(genes_c), "\n genes"),
               row_dend_width = unit(15, "mm"), use_raster = F,
               width = unit(90,'mm'), height = unit(length(genes_c)/12, 'mm'),
               show_heatmap_legend = F, bottom_annotation = ha)
col = colorRamp2(c(-2, 0, 1, 2, 6), c("cornflowerblue","white","orange",'orangered',"red4"))
lgd = Legend( col_fun = col, title = "log(FC)",legend_height = unit(30,"mm"),grid_width = unit(4,"mm"),
              title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 11),title_gap = unit(4,"mm"))
seperate_hm = grid.grabExpr(draw(ht_a %v% ht_b %v% ht_c, annotation_legend_list = lgd))
jpeg("figures/seperate_heatmaps.jpeg",height = 155, width = 300, units = "mm",res = 300)
grid.arrange(seperate_hm, fig_sliding_hypoxia, 
             layout_matrix = rbind(c(1,2),c(1,2)))
grid.text('time (hour)', x = unit(145,'mm'),y = unit(5,'mm'))
dev.off()



