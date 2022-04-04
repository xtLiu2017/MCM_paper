#plot the synthetic profile after fixing 
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
setwd("~/Documents/MCM_paper/data")
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
lgd = Legend( col_fun = col, title = "log(FC)",legend_height = unit(30,"mm"),grid_width = unit(4,"mm"),
             title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 11),title_gap = unit(4,"mm"))
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
jpeg("figures/fig3.a.jpeg",height = 150, width = 180, units = "mm",res = 300)
draw(ht1, annotation_legend_list = lgd)
grid.text('time (hour)', x = unit(156,'mm'),y = unit(5,'mm'))
dev.off()
fig3a_plot = grid.grabExpr(draw(ht1,annotation_legend_list = lgd)) 


setwd("~/Documents/MCM_paper/data")
#data preprocessing
Daisuke_info = read.csv(file = "preliminary.gamma.distr.fit.a.pt.csv",header = T)
genes_Dais = Daisuke_info[,1]
load("peak_info_MCM_spline.Rdata")
load("select.highquality.genes.Rdata")
genes_intersect = intersect(genes_Dais,select.highquality.genes)
rownames(Daisuke_info) = genes_Dais
colnames(peak_info) = c('pt1_i1','pl1_i1','pt2_i1','pl2_i1')
colnames(peak_info_set2) = c('pt1_i2','pl1_i2','pt2_i2','pl2_i2')

peak_time_Dais = Daisuke_info[genes_intersect,3]
genes_good = genes_intersect[which(peak_time_Dais < 7)]
genes_bad = genes_intersect[which(peak_time_Dais > 7)]

load("genes_to_fix.Rdata")
load("peak_info_fixed.Rdata")
genes_to_fix_here = intersect(genes_to_be_fixed,genes_good)
genes_to_fix_info = peak_info_fixed[genes_to_fix_here,] 
peak_info_after_fix = peak_info[genes_good,]
peak_info_after_fix[genes_to_fix_here,] = genes_to_fix_info
colnames(peak_info_after_fix) = c("pt1_fix",'pl1_fix','pt2_fix','pl2_fix')

load("merged_result.Rdata")
parameter_best_1 = parameter_best
rownames(parameter_best_1) = rownames(peak_info)
load("merged_result_mid.Rdata")
parameter_best_2 = parameter_best
rownames(parameter_best_2) = rownames(peak_info_set2)
load("merged_result.Rdata")
rownames(parameter_best) = rownames(peak_info)
df = data.frame(peak_info[genes_good,], peak_info_set2[genes_good,],peak_info_after_fix,
                a1_k_i1 = parameter_best_1[genes_good,1]/parameter_best_1[genes_good,3],
                amp_Dais = Daisuke_info[genes_good,2], pt_Dais = Daisuke_info[genes_good,3])
length(genes_bad)
length(genes_good)
#plotting, first peak time, MCM vs Daisuke's
l2_l1 = data.matrix(df['pl2_fix']/df['pl1_fix'][,1])[,]
sp_genes = names(which(l2_l1 < 0.3 & l2_l1 > 0))
df['sp_genes'] = l2_l1 < 0.3 & l2_l1 > 0

library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(tidytext)
pcc_1 = cor(df["pt1_fix"],df["pt_Dais"],use = "pair")
scc_1 = cor(df["pt1_fix"],df["pt_Dais"],use = "pair",method = "spearman")
p1 = ggplot(df,aes(y = pt1_fix + 3, x= pt_Dais)) + geom_bin2d(bins = 100) + theme_bw() + 
  labs(y = "first peak time (MCM, set1)",x = substitute(paste('peak time', italic(' (In planta)'))),size = 20)+ 
  annotate("text", x = 4.2, y=12, label = paste("PCC =",round(pcc_1,2)),size = 4)+
  annotate("text", x = 4.2, y=11.1, label = paste("SCC =",round(scc_1,2)),size = 4)+
  scale_fill_distiller(palette = "Spectral",direction = -1, limits = c(1,7))+
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=13),
        axis.text = element_text(size = 11),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        plot.margin = unit(c(4, 4, 5, 3.5), 'mm'), legend.position = 'none') + 
  geom_abline(slope=1,intercept = 0,size=0.5,colour = "grey") + ggtitle("first peak")

pcc_2 = cor(df["pl1_fix"],df["amp_Dais"],use = "pair")
scc_2 = cor(df["pl1_fix"],df["amp_Dais"],use = "pair",method = "spearman")
p2 = ggplot(df,aes(y = pl1_fix, x= amp_Dais)) + geom_bin2d(bins = 100) + theme_bw() + 
  labs(y = "first peak level (MCM, set1)",x = substitute(paste('peak level', italic(' (In planta)'))),size = 20) + 
  annotate("text", x = 4, y=5.55, label = paste("PCC =",round(pcc_2,2)),size = 4)+
  annotate("text", x = 4, y=5, label = paste("SCC =",round(scc_2,2)),size = 4)+
  scale_fill_distiller(palette = "Spectral",direction = -1, limits = c(1,7))+
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=13),
        axis.text = element_text(size = 11),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        plot.margin = unit(c(1, 4, 8, 3.5), 'mm'),legend.position = "none") + ggtitle("first peak")


pcc_3 = cor(df["pt2_fix"],df["pt_Dais"],use = "pair")
scc_3 = cor(df["pt2_fix"],df["pt_Dais"],use = "pair",method = "spearman")
p3 = ggplot(df,aes(y = pt2_fix + 3, x= pt_Dais)) + geom_bin2d(bins = 100) + theme_bw() + 
  labs(y = "second peak time (MCM, set1)",x = substitute(paste('peak time', italic(' (In planta)'))),size = 20)+ 
  annotate("text", x = 4.2, y=23, label = paste("PCC =",round(pcc_3,2)),size = 4)+
  annotate("text", x = 4.2, y=21.5, label = paste("SCC =",round(scc_3,2)),size = 4)+
  scale_fill_distiller(palette = "Spectral",direction = -1, limits = c(1,7))+
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=13),
        axis.text = element_text(size = 11),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        plot.margin = unit(c(4, 2, 5, 5), 'mm')) +
  geom_abline(slope=1,intercept = 0,size=0.5,colour = "grey") + ggtitle("second peak")
pcc_4 = cor(df["pl2_fix"],df["amp_Dais"],use = "pair")
scc_4 = cor(df["pl2_fix"],df["amp_Dais"],use = "pair",method = "spearman")
p4 = ggplot(df,aes(y = pl2_fix, x= amp_Dais)) + geom_bin2d(bins = 100) + theme_bw() + 
  labs(y = "second peak level (MCM, set1)",x = substitute(paste('peak level', italic(' (In planta)'))),size = 20)+ 
  annotate("text", x = 4.2, y=3.2, label = paste("PCC =",round(pcc_4,2)),size = 4)+
  annotate("text", x = 4.2, y=2.87, label = paste("SCC =",round(scc_4,2)),size = 4)+
  scale_fill_distiller(palette = "Spectral",direction = -1, limits = c(1,7))+
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=13),
        axis.text = element_text(size = 11),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        plot.margin = unit(c(1, 2, 8, 5), 'mm')) + ggtitle("second peak")

# get the peak time distribution across 1939 genes,
# there are definitely genes without peaks
peak_info_1939 = peak_info[select.highquality.genes,]
fixed_genes_1939 = intersect(select.highquality.genes, genes_to_be_fixed)
peak_info_1939_fixed = peak_info_1939
peak_info_1939_fixed[fixed_genes_1939,] = peak_info_fixed[fixed_genes_1939,]

mydf_pt1pt2 = data.frame(cbind(c(peak_info_1939_fixed[,'pt1_i1'], peak_info_1939_fixed[,'pt2_i1']),
                               rep(c(1,2), each = 1939)))
colnames(mydf_pt1pt2) = c('pt','X')
mydf_pt1pt2$X = as.factor(mydf_pt1pt2$X)
first_peak_number = 1939 - sum(is.na(peak_info_1939_fixed[,1]))
second_peak_number = 1939 - sum(is.na(peak_info_1939_fixed[,3]))
p <- ggplot(mydf_pt1pt2, aes(x = X, y=pt)) + geom_boxplot(width = 0.5) + theme_bw() + xlab('') +
  ylab('peak time (hour)') + scale_x_discrete(breaks = c(1,2),labels = c('first peak','second peak')) +
  theme(axis.text = element_text(size = 10.5), 
        axis.text.y = element_text(colour = c('black','black','red','black','black','black','black'))) + 
  geom_hline(yintercept=9, linetype="dashed", color = "red", size=1) +
  scale_y_continuous(breaks = c(0,5,9,10,15,20,25), labels = c(0,5,9,10,15,20,25), limits = c(0,25)) + 
  annotate("text", x = 1, y = 23, label = paste(first_peak_number,'\n peaks'), size = 3) + 
  annotate("text", x = 2, y = 23, label = paste(second_peak_number,'\n peaks'),  size = 3) + 
  ggtitle('1939 well-fit genes')
  
  
jpeg('figures/peak_times.jpeg', height = 100, width = 60, res = 300, unit = 'mm')
p
dev.off()

#look at overlap of genes with Daisuke's data (expressed in cell pop 1)
sliding_window_func = function(sorted_peak_level_ratio, genes_Dais, window){
  No_overlap_vec = c()
  peak_ratio = c()
  for (mystart in 1:(length(sort_peak_level_ratio) - window)){
    genes_window = names(sort_peak_level_ratio[mystart: (mystart + window - 1)])
    overlap_genes = intersect(genes_window, genes_Dais)
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
                                    genes_Dais = genes_Dais, window = window1)

df_fig3b = data.frame(index = 1:dim(No.overlap_30)[2], ratio = No.overlap_30[1,], No = No.overlap_30[2,])
myat = c()
for (ii in 1/c(0.25,0.5,1,1.25,2.5)){
  myat = c(myat, which.min(abs(ii - No.overlap_30[1,])))
}
Dais_gene_no. = length(intersect(genes_Dais,names(sort_peak_level_ratio)))
MCM_gene_no. = length(sort_peak_level_ratio)
lower_0.95 = qhyper(p=0.025, m = Dais_gene_no., n=MCM_gene_no. - Dais_gene_no., k = window1, lower.tail = TRUE, log.p = FALSE)
upper_0.95 = qhyper(p=0.975, m = Dais_gene_no., n=MCM_gene_no. - Dais_gene_no., k = window1, lower.tail = TRUE, log.p = FALSE)
mycol <- rgb(142,229,238, max = 255, alpha = 125)

fig3b = ggplot(df_fig3b, aes(x = index, y = No)) + geom_line() + theme_bw() + xlab('peak 1 level / peak 2 level') +
  ylab(bquote(atop('No. of genes induced', 'by ' ~ italic('In planta') ~ ' AvrRpt2'))) + 
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        plot.margin = unit(c(15, 15, 6, 10), 'mm'), axis.title = element_text(size = 15),
        axis.text = element_text(size = 13), panel.grid = element_blank()) + 
  scale_x_continuous(breaks = myat, labels = 1/c(0.25,0.5,1,1.25,2.5)) + 
  geom_hline(yintercept=window1 * Dais_gene_no. / MCM_gene_no. , linetype = "dashed", col = "blue") + 
  annotate("rect", xmin = -50, ymin = lower_0.95, xmax = 1600, ymax = upper_0.95, alpha=0.5, fill=mycol) +
  annotate("text", x = 1270, y = 16.5, label = 'expected number with \n 95% confidence interval',
           size = 3.6, col = "blue") + 
  annotate('text',x = 800, y = 3.5,  label = 'window size = 30 genes', size = 5)

library(clusterProfiler)
library(org.At.tair.db)
sp_genes = names(which(peak_level_ratio > 4))
#for single peak genes
mean(peak_info_1939[sp_genes, 'pt1_i1']) + 3
sd(peak_info_1939[sp_genes, 'pt1_i1'])
#for all genes in the background
mean(peak_info_1939[, 'pt1_i1'],na.rm = T) + 3
sd(peak_info_1939[, 'pt1_i1'],na.rm = T) 

sp.BP <- enrichGO(gene = sp_genes,universe = select.highquality.genes,
                  OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                  pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = F)


top12_terms_sp = data.frame(sp.BP)[1:12,1]
df_GO_HSF = data.frame(terms_ID = c(top12_terms_sp), 
                       terms = data.frame(sp.BP)[top12_terms_sp,'Description'],
                       p_val = data.frame(sp.BP)[top12_terms_sp,'pvalue'])

ggplot_sp = ggplot(df_GO_HSF,aes(x = terms, y = -log(p_val))) + 
  geom_bar(stat="identity",position=position_dodge(),width = 0.7) +
  coord_flip() + theme_bw() + xlab('GO terms') + ylab('-log(p-value)') + 
  ggtitle("single peak genes") + 
  theme(panel.border = element_blank(), axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),plot.title = element_text(size=16,hjust=0.5,face = 'bold'),
        plot.margin = unit(c(0,6,4,0),'mm')) + 
  scale_x_discrete(limits=rev(df_GO_HSF$terms[1:12])) 


lay <- rbind(c(1,1,1,1,1,1,1,1,NA,3,3,3,3,4,4,4,4,4),
             c(1,1,1,1,1,1,1,1,NA,3,3,3,3,4,4,4,4,4),
             c(1,1,1,1,1,1,1,1,NA,3,3,3,3,4,4,4,4,4),
             c(1,1,1,1,1,1,1,1,NA,3,3,3,3,4,4,4,4,4),
             c(1,1,1,1,1,1,1,1,NA,3,3,3,3,4,4,4,4,4),
             c(1,1,1,1,1,1,1,1,NA,5,5,5,5,6,6,6,6,6),
             c(1,1,1,1,1,1,1,1,NA,5,5,5,5,6,6,6,6,6),
             c(1,1,1,1,1,1,1,1,NA,5,5,5,5,6,6,6,6,6),
             c(1,1,1,1,1,1,1,1,NA,5,5,5,5,6,6,6,6,6),
             c(2,2,2,2,2,2,2,2,NA,5,5,5,5,6,6,6,6,6),
             c(2,2,2,2,2,2,2,2,NA,7,7,7,7,7,7,7,7,7),
             c(2,2,2,2,2,2,2,2,NA,7,7,7,7,7,7,7,7,7),
             c(2,2,2,2,2,2,2,2,NA,7,7,7,7,7,7,7,7,7),
             c(2,2,2,2,2,2,2,2,NA,7,7,7,7,7,7,7,7,7),
             c(2,2,2,2,2,2,2,2,NA,7,7,7,7,7,7,7,7,7))
             

jpeg('figures/figure3.jpeg', width = 320, height = 240, unit = 'mm', res = 500)
grid.arrange(fig3a_plot,fig3b,p1,p3,p2,p4,ggplot_sp,layout_matrix = lay)
grid.text('hour', x = unit(132,'mm'),y = unit(98,'mm'))
grid.text('A', x = unit(10,'mm'),y = unit(235,'mm'),gp = gpar(fontface = 'bold'))
grid.text('B', x = unit(10,'mm'),y = unit(84,'mm'),gp = gpar(fontface = 'bold'))
grid.text('C', x = unit(150,'mm'),y = unit(235,'mm'),gp = gpar(fontface = 'bold'))
grid.text('D', x = unit(150,'mm'),y = unit(84,'mm'),gp = gpar(fontface = 'bold'))
dev.off()




