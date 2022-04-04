setwd("/Users/liux/Documents/MCM_paper/data")
load("pvalues_sliding_windows.Rdata")
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(tidytext)
padj_first = p.adjust(c(p_first_all,p_first_max), method = "BH")
padj_second = p.adjust(c(p_second_all,p_second_max), method = "BH")
FDR_second = -log10(c(p_second_all,p_second_max)[which.min(abs(padj_second - 0.05))])
FDR_first = -log10(c(p_first_all,p_first_max)[which.min(abs(padj_first - 0.05))])
FDR_mean = round(mean(c(FDR_first,FDR_second)), digits = 1)

df_first = data.frame(TF_labels = toupper(names(sort(p_first_max[prioritized_first]))), 
                      pval = -log10(sort(p_first_max[prioritized_first])), peak = 1)
df_second = data.frame(TF_labels = toupper(names(sort(p_second_max[prioritized_second]))), 
                       pval = -log10(sort(p_second_max[prioritized_second])), peak = 2)
df_fig4a = rbind(df_first,df_second)
tibble_fig4a = as.tibble(df_fig4a)
myfacet_label_4a = c('first peak (23 TFs, FDR = 0.05)','second peak (13 TFs, FDR = 0.05)')
names(myfacet_label_4a) = c(1,2)

fig4a_plot = tibble_fig4a %>% mutate(peak=as.factor(peak),
                                     TF_labels=reorder_within(TF_labels,pval, peak))  %>% 
  ggplot(aes(x = TF_labels, y=pval)) + 
  geom_bar(stat="identity",width = 0.7) + 
  facet_wrap( ~ peak , scales = 'free',labeller = labeller(peak = myfacet_label_4a), nrow = 2) + 
  coord_flip() + theme_bw() + ylab(expression('-log'[10]*' (p-value)')) + xlab("") + 
  ylim(0,9) + 
  theme(panel.border = element_blank(),axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7), axis.title.x = element_text(size = 12),
        strip.text.x = element_text(size=15), strip.background = element_blank()) + 
  scale_x_reordered() 

#---------
load("Annotation_mat.Rdata")
tfs = rownames(Annotation_mat)
tf_labels = unlist(lapply(tfs,function(x){strsplit(x,split = "_")[[1]][1]}))
names(tf_labels) = tfs
genes_union = union(prioritized_first, prioritized_second)
tf_union = c()
for (mytf in genes_union){
  tf_union = c(tf_union, names(which(tf_labels == mytf)))
}

Dense_sliding_window_func = function(x, y, window, start = 1){
  #x -- peak time, y -- binding site number, window -- local number of gene
  #x should be sorted, y should match x
  #check if x and y match
  if (sum(names(x) == names(y)) != length(x)){
    print("error")
  }
  #names of x and y should be the same
  ordered_gene_names = names(x)
  ordered_tf_number = y
  times = c()
  output = c()
  for (i in 1:(length(x) - window)){
    times = c(times, x[i + as.integer(window/2)] )
    gene_names_window = ordered_gene_names[i: (i + window - 1)]
    output = c(output, sum(y[gene_names_window] > 0))
  }
  return(rbind(times, output))
}
library(qvalue)
load("peak_info_all.Rdata")
load("select.highquality.genes.Rdata")
load("genes_to_fix.Rdata")
peak_info_update = peak_info
peak_info_update[genes_to_be_fixed,] = peak_info_fixed[genes_to_be_fixed,]

first_peaks = peak_info_update[select.highquality.genes,1]
second_peaks = peak_info_update[select.highquality.genes,3]
ordered_first_peaks = sort(first_peaks)
ordered_second_peaks = sort(second_peaks)

window1 = 150
time_mat = c()
gene_number_mat = c()
p_mat = c()
for (mytf in tf_union){
  my_TF_binding_first = Annotation_mat[mytf,names(ordered_first_peaks)]
  my_TF_binding_second = Annotation_mat[mytf,names(ordered_second_peaks)]
 
  my_TF_binding_all = Annotation_mat[mytf, select.highquality.genes]
  m = sum(my_TF_binding_all> 0) 
  n = length(my_TF_binding_all) - m
  
  time_first = Dense_sliding_window_func(x = ordered_first_peaks, y = my_TF_binding_first,window = window1,start = 1)[1,] + 3
  mydensity_first = Dense_sliding_window_func(x = ordered_first_peaks, y = my_TF_binding_first,window = window1,start = 1)[2,] 
  time_second = Dense_sliding_window_func(x = ordered_second_peaks, y = my_TF_binding_second,window = window1,start = 1)[1,] + 3
  mydensity_second = Dense_sliding_window_func(x = ordered_second_peaks, y = my_TF_binding_second,window = window1,start = 1)[2,]
  
  mymax = max(max(mydensity_first), max(mydensity_second))
  mymin = min(min(mydensity_first), min(mydensity_second))
  mydensity_first_sd = (mydensity_first - mymin )/  (mymax - mymin)
  mydensity_second_sd = (mydensity_second - mymin)/  (mymax - mymin)
  
  p_first = 1 - phyper(mydensity_first - 1, k = window1, m = m, n = n, lower.tail = TRUE, log.p = FALSE)
  max_p_first = 1 - phyper(max(mydensity_first) - 1, k = window1, m = m, n = n, lower.tail = TRUE, log.p = FALSE)
  max_p_first_pos = which.max(mydensity_first)
  max_time_first = time_first[max_p_first_pos] + 3
  
  p_second = 1 - phyper(mydensity_second - 1, k = window1, m = m, n = n, lower.tail = TRUE, log.p = FALSE)
  max_p_second = 1 - phyper(max(mydensity_second) - 1, k = window1, m = m, n = n, lower.tail = TRUE, log.p = FALSE)
  max_p_second_pos = which.max(mydensity_second)
  max_time_second = time_second[max_p_second_pos] + 3
  
  time_mat = rbind(time_mat, c(max_time_first,max_time_second))
  p_mat = rbind(p_mat, c( p_first,  p_second) )
  gene_number_mat = rbind(gene_number_mat, c(mydensity_first, mydensity_second))
}
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
genes_union_upp = toupper(genes_union)
rownames(gene_number_mat) = genes_union_upp
rownames(p_mat) = genes_union_upp
rowname_col = c()
for (i in genes_union_upp){
  if (grepl('HSF', i)){
    myrowname = 'red'
  }else if(grepl('ANAC',i) | grepl('NAP', i) | grepl('SND3',i)){
    myrowname = 'blue'
  }else if(grepl('WRKY',i)){
    myrowname = 'orange'
  }else{
    myrowname = 'gray30'
  }
  rowname_col = c(rowname_col, myrowname)
}

col = colorRamp2(seq(0,9, length.out = 8), brewer.pal(9,'Oranges')[1:8])
myat_first = c()
for (ii in c(4,5,6,7,8)){
  myat_first = c(myat_first, which.min(abs(ii - time_first)))
}
myat_second = c()
for (ii in c(10,14,17,20,22)){
  myat_second = c(myat_second, which.min(abs(ii - time_second)))
}
# get the labels for columns, time scale
myalpha = rep(0, dim(p_mat)[2])
myalpha[myat_first] = 1
myalpha[myat_second + length(time_first)] = 1
mylabels = rep('', dim(p_mat)[2])
mylabels[myat_first] = c('4','5','6','7','8')
mylabels[myat_second + length(time_first)] = c('10','14','17','20','22')

ha = HeatmapAnnotation(time = anno_lines(c(0.001,rep(0,dim(p_mat)[2]-1)), ylim = c(-1,1),
                                        border = F, axis = F, gp = gpar(lwd = 2),add_points = T,
                                        pt_gp = gpar(alpha = myalpha)),
                       foo = anno_text(mylabels, rot = 0, just = 'center',gp = gpar(fontsize = 11)),
                       show_annotation_name=F)
#clustering the rows
dist_mat = dist(-log10(p_mat))
mhclust = hclust(dist_mat,method = "average")  

lgd = Legend(col_fun = col, title = expression('-log'[10]*' (p-value)'),
             at = c(0, 2, 5, 7, 9, FDR_mean), labels = c('0','2','5','7','9',paste(FDR_mean,'(FDR = 0.05)')),
             legend_height = unit(30,"mm"),legend_width = unit(8,"mm"),title_position = 'topleft',
             title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 8),title_gap = unit(3,"mm"))
ht1 = Heatmap(-log10(p_mat), col = col,cluster_columns = FALSE, 
              heatmap_height = unit(14, "cm"),
              heatmap_width = unit(14.5, "cm"), cluster_rows = mhclust,
        name = 'j', row_dend_gp = gpar(), show_row_names = T, show_column_names = F,
        show_heatmap_legend = F, row_title = paste(dim(p_mat)[1],"TFs"), row_title_rot = 0, 
        use_raster = F,
        column_split = c(rep("first peak", times = length(p_first)),rep("second peak", times = length(p_second))),
        column_gap = unit(5,"mm"),bottom_annotation = ha,
        row_names_gp = gpar(col = rowname_col,fontsize = 10))

grob.fig4b = grid.grabExpr(draw(ht1,heatmap_legend_list = lgd)) 
jpeg('figures/figure4b.jpeg', width = 200, height = 160, unit = 'mm', res = 300)
grid.arrange(grob.fig4b)
grid.text('time (hour)', x = unit(150,'mm'),y = unit(13.1,'mm'))
dev.off()

col_number = colorRamp2(seq(0,max(gene_number_mat), length.out = 6), 
                 rev(c('red','orange','yellow','green','deepskyblue','blue')))
rownames(gene_number_mat) = genes_union_upp
number_of_columns = dim(gene_number_mat)[2]
PCC_vec = c()
for (i in 1:dim(gene_number_mat)[1]){
  early_one = gene_number_mat[i,1:(number_of_columns/2)]
  late_one = gene_number_mat[i,(number_of_columns/2+1):number_of_columns]
  PCC_vec = c(PCC_vec, cor(early_one, late_one))
}
jpeg('figures/sup.fig9.jpeg', width = 160, height = 160, unit = 'mm', res = 400)
lgd_number = Legend(col_fun = col_number, title = 'gene number',
             legend_height = unit(30,"mm"),legend_width = unit(8,"mm"),title_position = 'topleft',
             title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8),title_gap = unit(3,"mm"))

ht_number = Heatmap(gene_number_mat, col = col_number, cluster_columns = FALSE, name = 'gene number',
        heatmap_height = unit(14, "cm"),heatmap_width = unit(13, "cm"),show_heatmap_legend = F,
         row_dend_gp = gpar(), show_row_names = T, show_column_names = F,cluster_rows = mhclust,
        row_title = paste(dim(gene_number_mat)[1],"TFs"), row_title_rot = 0, use_raster = F,
         column_split = c(rep("first peak", times = length(time_first)),rep("second peak", times = length(time_second))),
         column_gap = unit(5,"mm"), row_names_gp = gpar(col = rowname_col,fontsize = 10),
        bottom_annotation = ha)
draw(ht_number,heatmap_legend_list = lgd_number)
decorate_annotation("foo", { 
  grid.text("time (hour)", y = unit(1, "mm"), x= unit(95,'mm'), just = "bottom") 
})
dev.off()



#generate temporal enrichment for HSF, ANAC and WRKY
HSFs = c('HSFC1_col_a','HSF3_col_a','HSFA6B_col','HSF6_col_a','HSF7_col_a')
ANACs = c('ANAC047_col','ANAC083_col_v3a','ANAC058_col_a','ANAC062_col','ANAC055_col','ANAC038_col_a')
WRKYs = c('WRKY7_col','WRKY75_col_a','WRKY15_col_b','WRKY70_col','WRKY24_col_a')
HSF_mat_1 = c()
HSF_mat_2 = c()
for (mytf in HSFs){
  my_TF_binding_first = Annotation_mat[mytf,names(ordered_first_peaks)]
  my_TF_binding_second = Annotation_mat[mytf,names(ordered_second_peaks)]
  
  time_first = Dense_sliding_window_func(x = ordered_first_peaks, y = my_TF_binding_first,window = window1,start = 1)[1,] + 3
  mydensity_first = Dense_sliding_window_func(x = ordered_first_peaks, y = my_TF_binding_first,window = window1,start = 1)[2,] 
  time_second = Dense_sliding_window_func(x = ordered_second_peaks, y = my_TF_binding_second,window = window1,start = 1)[1,] + 3
  mydensity_second = Dense_sliding_window_func(x = ordered_second_peaks, y = my_TF_binding_second,window = window1,start = 1)[2,]
  HSF_mat_1 = rbind(HSF_mat_1, mydensity_first)
  HSF_mat_2 = rbind(HSF_mat_2, mydensity_second)
}
rownames(HSF_mat_1) = HSFs
rownames(HSF_mat_2) = HSFs

ANAC_mat_1 = c()
ANAC_mat_2 = c()
for (mytf in ANACs){
  my_TF_binding_first = Annotation_mat[mytf,names(ordered_first_peaks)]
  my_TF_binding_second = Annotation_mat[mytf,names(ordered_second_peaks)]
  
  time_first = Dense_sliding_window_func(x = ordered_first_peaks, y = my_TF_binding_first,window = window1,start = 1)[1,] + 3
  mydensity_first = Dense_sliding_window_func(x = ordered_first_peaks, y = my_TF_binding_first,window = window1,start = 1)[2,] 
  time_second = Dense_sliding_window_func(x = ordered_second_peaks, y = my_TF_binding_second,window = window1,start = 1)[1,] + 3
  mydensity_second = Dense_sliding_window_func(x = ordered_second_peaks, y = my_TF_binding_second,window = window1,start = 1)[2,]
  ANAC_mat_1 = rbind(ANAC_mat_1, mydensity_first)
  ANAC_mat_2 = rbind(ANAC_mat_2, mydensity_second)
}
rownames(ANAC_mat_1) = ANACs
rownames(ANAC_mat_2) = ANACs

WRKY_mat_1 = c()
WRKY_mat_2 = c()
for (mytf in WRKYs){
  my_TF_binding_first = Annotation_mat[mytf,names(ordered_first_peaks)]
  my_TF_binding_second = Annotation_mat[mytf,names(ordered_second_peaks)]
  
  time_first = Dense_sliding_window_func(x = ordered_first_peaks, y = my_TF_binding_first,window = window1,start = 1)[1,] + 3
  mydensity_first = Dense_sliding_window_func(x = ordered_first_peaks, y = my_TF_binding_first,window = window1,start = 1)[2,] 
  time_second = Dense_sliding_window_func(x = ordered_second_peaks, y = my_TF_binding_second,window = window1,start = 1)[1,] + 3
  mydensity_second = Dense_sliding_window_func(x = ordered_second_peaks, y = my_TF_binding_second,window = window1,start = 1)[2,]
  WRKY_mat_1 = rbind(WRKY_mat_1, mydensity_first)
  WRKY_mat_2 = rbind(WRKY_mat_2, mydensity_second)
}
rownames(WRKY_mat_1) = WRKYs
rownames(WRKY_mat_2) = WRKYs

myat_first = c()
for (ii in c(4,5,6,7,8)){
  myat_first = c(myat_first, which.min(abs(ii - time_first)))
}
names(myat_first) = NULL
myat_second = c()
for (ii in c(10,14,16,20,22)){
  myat_second = c(myat_second, which.min(abs(ii - time_second)))
}
names(myat_second) = NULL

breaks_fun <- function(x) {
  if (max(x) > 1760) {
    myat_first
  } else {
    myat_second
  }
}
label_func = function(x){
  x_cha = as.character(x)
  a = c('4','5','6','7','8')
  names(a) = as.character(myat_first)
  b = c('10','14','16','20','22')
  names(b) = as.character(myat_second)
  c = c(a,b)
  c[x_cha]
}
myfacet_label = c('first peak','second peak')

HSF_mat_1_t = t(HSF_mat_1)
HSF_mat_2_t = t(HSF_mat_2)
HSF_df = data.frame(no = c(as.vector(HSF_mat_1_t), as.vector(HSF_mat_2_t)))
rownames(HSF_df) = NULL
HSF_df$time = as.vector(c(rep(time_first, dim(HSF_mat_1)[1]), rep(time_second, dim(HSF_mat_2)[1])))
HSF_df$peak = c(rep(1, length(time_first) * dim(HSF_mat_1)[1]), rep(2, length(time_second) * dim(HSF_mat_2)[1]))
HSF_df$gene = c(rep(HSFs, each = length(time_first)),rep(HSFs, each = length(time_second)))
HSF_df$index = c(rep(seq(1, length(time_first)), dim(HSF_mat_1)[1]), rep(seq(1, length(time_second)), dim(HSF_mat_2)[1]))

col_HSF = brewer.pal(9,'Reds')[3:7]

HSF_label = unlist(lapply(HSFs,function(x){strsplit(x,split = "_")[[1]][1]}))
names(HSF_label) = HSFs
names(myfacet_label) = c(1,2)

ANAC_mat_1_t = t(ANAC_mat_1)
ANAC_mat_2_t = t(ANAC_mat_2)
ANAC_df = data.frame(no = c(as.vector(ANAC_mat_1_t), as.vector(ANAC_mat_2_t)))
rownames(ANAC_df) = NULL
ANAC_df$time = as.vector(c(rep(time_first, dim(ANAC_mat_1)[1]), rep(time_second, dim(ANAC_mat_2)[1])))
ANAC_df$peak = c(rep(1, length(time_first) * dim(ANAC_mat_1)[1]), rep(2, length(time_second) * dim(ANAC_mat_2)[1]))
ANAC_df$gene = c(rep(ANACs, each = length(time_first)),rep(ANACs, each = length(time_second)))
ANAC_df$index = c(rep(seq(1, length(time_first)), dim(ANAC_mat_1)[1]), rep(seq(1, length(time_second)), dim(ANAC_mat_2)[1]))

col_ANAC = brewer.pal(9,'Blues')[3:8]

ANAC_label = unlist(lapply(ANACs,function(x){strsplit(x,split = "_")[[1]][1]}))
names(ANAC_label) = ANACs
names(myfacet_label) = c(1,2)

WRKY_mat_1_t = t(WRKY_mat_1)
WRKY_mat_2_t = t(WRKY_mat_2)
WRKY_df = data.frame(no = c(as.vector(WRKY_mat_1_t), as.vector(WRKY_mat_2_t)))
rownames(WRKY_df) = NULL
WRKY_df$time = as.vector(c(rep(time_first, dim(WRKY_mat_1)[1]), rep(time_second, dim(WRKY_mat_2)[1])))
WRKY_df$peak = c(rep(1, length(time_first) * dim(WRKY_mat_1)[1]), rep(2, length(time_second) * dim(WRKY_mat_2)[1]))
WRKY_df$gene = c(rep(WRKYs, each = length(time_first)),rep(WRKYs, each = length(time_second)))
WRKY_df$index = c(rep(seq(1, length(time_first)), dim(WRKY_mat_1)[1]), rep(seq(1, length(time_second)), dim(WRKY_mat_2)[1]))

col_WRKY = brewer.pal(9,'Oranges')[3:7]

WRKY_label = unlist(lapply(WRKYs,function(x){strsplit(x,split = "_")[[1]][1]}))
names(WRKY_label) = WRKYs
names(myfacet_label) = c(1,2)

HSF_df$family = 'HSF'
ANAC_df$family = 'ANAC'
WRKY_df$family = 'WRKY'
family_df = rbind(HSF_df,ANAC_df,WRKY_df)
family_df$family = factor(family_df$family, levels = c("HSF",'ANAC','WRKY'))

fig4c = ggplot(family_df, aes(x = index, y = no, color = gene)) + geom_line() + theme_bw() +
  facet_wrap( peak ~ family, scales = 'free_x',labeller = labeller(peak = myfacet_label)) +
  xlab("time (hour)") + ylab("No. of genes with at least 1 binding site") + 
  scale_x_continuous(breaks=breaks_fun, labels = label_func) +
  scale_color_manual(breaks = c(HSFs,ANACs,WRKYs), values=c(col_HSF, col_ANAC,col_WRKY), label = c(HSF_label, ANAC_label, WRKY_label)) +
  theme(legend.text = element_text(size=10),legend.title = element_blank(),
        strip.background = element_blank(),strip.text.x = element_text(size = 14),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10))



lay <- rbind(c(1,1,1,2,2,2,2,2,NA,NA,NA,NA,NA,NA),
             c(1,1,1,2,2,2,2,2,3,3,3,3,3,3),
             c(1,1,1,2,2,2,2,2,3,3,3,3,3,3),
             c(1,1,1,2,2,2,2,2,3,3,3,3,3,3),
             c(1,1,1,2,2,2,2,2,3,3,3,3,3,3),
             c(1,1,1,2,2,2,2,2,3,3,3,3,3,3),
             c(1,1,1,2,2,2,2,2,3,3,3,3,3,3),
             c(1,1,1,2,2,2,2,2,3,3,3,3,3,3),
             c(1,1,1,2,2,2,2,2,3,3,3,3,3,3),
             c(1,1,1,2,2,2,2,2,3,3,3,3,3,3),
             c(1,1,1,2,2,2,2,2,3,3,3,3,3,3),
             c(1,1,1,2,2,2,2,2,3,3,3,3,3,3),
             c(1,1,1,2,2,2,2,2,NA,NA,NA,NA,NA,NA))
jpeg('figures/figure4.jpeg', width = 490, height = 160, unit = 'mm', res = 500)
grid.arrange(fig4a_plot,grob.fig4b,fig4c,
             layout_matrix = lay)
grid.text('time (hour)', x = unit(241,'mm'),y = unit(13.1,'mm'))
grid.text('A', x = unit(6,'mm'),y = unit(156,'mm'),gp = gpar(fontface = 'bold'))
grid.text('B', x = unit(109,'mm'),y = unit(156,'mm'),gp = gpar(fontface = 'bold'))
grid.text('C', x = unit(283,'mm'),y = unit(156,'mm'),gp = gpar(fontface = 'bold'))
dev.off()


#plot done here


#genes with HSF binding site and peak time of enrichment
library(clusterProfiler)
library(org.At.tair.db)
HSF_genes1 = names(which(ordered_first_peaks < (5.5 - 3) & ordered_first_peaks > (4.5 - 3)))
HSF_genes2 = apply(Annotation_mat[HSFs, HSF_genes1],2,sum)
HSF_genes = names(which(HSF_genes2 > 0))
HSF.BP <- enrichGO(gene = HSF_genes,universe = select.highquality.genes,
                  OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                  pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = F)
HSF_genes_bs = names(which(apply(Annotation_mat[HSFs,select.highquality.genes],2,sum) > 0))
HSF_bs.BP <- enrichGO(gene = HSF_genes_bs,universe = select.highquality.genes,
                      OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                      pAdjustMethod = "BH", pvalueCutoff  = 0.2, qvalueCutoff  = 0.2, readable = F)
jpeg('figures/HSF.jpeg',height = 60, width = 60, unit = 'mm', res = 700)
plotGOgraph(HSF.BP)
dev.off()

top12_terms_HSF = data.frame(HSF.BP)[1:12,1]
df_GO_HSF = data.frame(terms_ID = c(top12_terms_HSF,top12_terms_HSF), 
                       terms = c(data.frame(HSF.BP)[top12_terms_HSF,'Description'],data.frame(HSF.BP)[top12_terms_HSF,'Description']),
                       p = c(data.frame(HSF.BP)[top12_terms_HSF,'pvalue'],
                       data.frame(HSF_bs.BP)[top12_terms_HSF,'pvalue']),
                       type = c(rep(c('1','2'),each = 12)))

ggplot_HSF = ggplot(df_GO_HSF,aes(x = terms, y = -log(p), fill = type)) + geom_bar(stat="identity",position=position_dodge(),width = 0.7) +
  coord_flip() + theme_bw() + xlab('GO terms') + ylab('-log(p-value)') + ggtitle("heat shock factors (HSFs)") + 
  theme(panel.border = element_blank(), axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),plot.title = element_text(size=7,hjust=0.5),
        legend.text = element_text(size=6), legend.key.size = unit(0.3, 'cm'),legend.position='none') + 
  guides(fill=guide_legend(title="")) + 
  scale_x_discrete(limits=rev(df_GO_HSF$terms[1:12])) + 
  scale_fill_discrete(limits = c('1','2'),labels = c('w/ peak time filtering','w/o peak time filtering'))



library("writexl")
write_xlsx(data.frame(HSF.BP),"HSF.BP.xlsx")

WRKY_genes1 = names(which(ordered_first_peaks < (6.5 - 3) & ordered_first_peaks > (5.5 - 3)))
WRKY_genes2 = names(which(ordered_second_peaks < (20 - 3) & ordered_second_peaks > (17 - 3)))
WRKY_genes3 = intersect(WRKY_genes1, WRKY_genes2)
WRKY_genes4 = apply(Annotation_mat[WRKYs, WRKY_genes3],2,sum)
WRKY_genes = names(which(WRKY_genes4 > 0))


myWRKY_genes1 = union(names(which(ordered_first_peaks < (4.8 - 3) & ordered_first_peaks > (4 - 3))),WRKY_genes1)
myWRKY_genes2 = union(names(which(ordered_second_peaks < (13 - 3) & ordered_second_peaks > (8 - 3))),WRKY_genes2)
myWRKY_genes3 = intersect(myWRKY_genes1, myWRKY_genes2 )
myWRKY_genes4 = apply(Annotation_mat[WRKYs, myWRKY_genes3],2,sum)
myWRKY_genes = names(which(myWRKY_genes4 > 0))


WRKY.BP <- enrichGO(gene = WRKY_genes,universe = select.highquality.genes,
                   OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                   pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = F)
myWRKY.BP <- enrichGO(gene = myWRKY_genes,universe = select.highquality.genes,
                      OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                      pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = F)


ANAC_genes1 = names(which(ordered_first_peaks > (7 - 3)))
ANAC_genes2 = names(which(ordered_second_peaks > (20 - 3)))
ANAC_genes3 = intersect(ANAC_genes1, ANAC_genes2)
ANAC_genes4 = apply(Annotation_mat[ANACs,ANAC_genes3], 2, sum)
ANAC_genes = names(which(ANAC_genes4 > 0))
ANAC.BP <- enrichGO(gene = ANAC_genes,universe = select.highquality.genes,
                    OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                    pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = F)

ERF_B1s = c('AT1G28160_col_a',"ERF4_col_a")
ERF_B1s_genes1 = names(which(ordered_first_peaks < (5.5 - 3) & ordered_first_peaks > (5 - 3)))
ERF_B1s_genes2 = apply(Annotation_mat[ERF_B1s, ERF_B1s_genes1],2,sum)
ERF_B1s_genes = names(which(ERF_B1s_genes2 > 0))
ERF_B1.BP = enrichGO(gene = ERF_B1s_genes,universe = select.highquality.genes,
                     OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                     pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = F)

MYB_likes = c('AT5G56840_col_a','At5g47390_col_b', 'At1g74840_col100_a' )
MYB_likes_genes1 = names(which(ordered_second_peaks < (15 - 3) & ordered_second_peaks > (10 - 3)))
MYB_likes_genes2 = apply(Annotation_mat[MYB_likes,MYB_likes_genes1],2,sum)
MYB_likes_genes = names(which(MYB_likes_genes2 > 0))
MYB_likes.BP = enrichGO(gene = MYB_likes_genes,universe = select.highquality.genes,
                        OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                        pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = F)


CAMTA1_genes1 = names(which(ordered_first_peaks < (6 - 3) & ordered_first_peaks > (4.5 - 3)))
CAMTA1_genes2 = Annotation_mat['CAMTA1_col_a',CAMTA1_genes1]
CAMTA1_genes = names(which(CAMTA1_genes2 > 0))
CAMTA1_BP = enrichGO(gene = CAMTA1_genes,universe = select.highquality.genes,
                     OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                     pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = F)
CAMTA1_genes_bs =  names(which(Annotation_mat['CAMTA1_col_a',select.highquality.genes] > 0))
CAMTA1_bs_BP = enrichGO(gene = CAMTA1_genes_bs,universe = select.highquality.genes,
                        OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                        pAdjustMethod = "BH", pvalueCutoff  = 0.3, qvalueCutoff  = 0.3, readable = F)

write_xlsx(data.frame(CAMTA1_BP),"CAMTA1.BP.xlsx")
jpeg('figures/CAMTA1.jpeg',height = 60, width = 60, unit = 'mm', res = 700)
plotGOgraph(CAMTA1_BP)
dev.off()

top11_terms_CAMTA1 = data.frame(CAMTA1_BP)[1:11,1]
df_GO_CAMTA1 = data.frame(terms_ID = c(top11_terms_CAMTA1,top11_terms_CAMTA1), 
                       terms = c(data.frame(CAMTA1_BP)[top11_terms_CAMTA1,'Description'],data.frame(CAMTA1_BP)[top11_terms_CAMTA1,'Description']),
                       p = c(data.frame(CAMTA1_BP)[top11_terms_CAMTA1,'pvalue'],
                             data.frame(CAMTA1_bs_BP)[top11_terms_CAMTA1,'pvalue']),
                       type = c(rep(c('1','2'),each = 11)))

ggplot_CAMTA1 = ggplot(df_GO_CAMTA1,aes(x = terms, y = -log(p), fill = type)) + geom_bar(stat="identity",position=position_dodge(),width = 0.7) +
  coord_flip() + theme_bw() + xlab('GO terms') + ylab('-log(p-value)') + ggtitle("calmodulin binding transcription activator 1 (CAMTA1)") + 
  theme(panel.border = element_blank(), axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),plot.title = element_text(size=7,hjust=0.5),
        legend.text = element_text(size=6), legend.key.size = unit(0.3, 'cm'),legend.position='bottom') + 
  guides(fill=guide_legend(title="")) + 
  scale_x_discrete(limits=rev(df_GO_CAMTA1$terms[1:11])) + 
  scale_fill_discrete(limits = c('1','2'),labels = c('w/ peak time filtering','w/o peak time filtering'))
margin_HSF = theme(plot.margin = unit(c(0.3,0.5,0.7,1), "cm"))
margin1_CAMTA1 = theme(plot.margin = unit(c(0.2,1,0.2,1), "cm"))
lay <- rbind(c(1,1,1,1,1),c(1,1,1,1,1),
             c(1,1,1,1,1),c(1,1,1,1,1),
             c(2,2,2,2,2),c(2,2,2,2,2),
             c(2,2,2,2,2),c(2,2,2,2,2),
             c(2,2,2,2,2))
jpeg('figures/HSF_GO_pvalues.jpeg',height = 120, width = 120, unit = 'mm', res = 300)
grid.arrange(ggplot_HSF ,ggplot_CAMTA1 ,nrow = 2, layout_matrix = lay,
             top = textGrob("Top significantly over-represented GO terms (FDR = 0.05)", gp=gpar(fontsize=11,fontface = 'bold')))
dev.off()

DEL2_genes1 = names(which(ordered_first_peaks < (7 - 3) & ordered_first_peaks > (5 - 3)))
DEL2_genes2 = Annotation_mat['DEL2_col_a', DEL2_genes1]
DEL2_genes = names(which(DEL2_genes2 > 0))
DEL2.BP = enrichGO(gene = DEL2_genes,universe = select.highquality.genes,
                   OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                   pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = F)

#save the genes regulated by TFs (families)
save(HSF_genes, WRKY_genes, ANAC_genes, ERF_B1s_genes,MYB_likes_genes, CAMTA1_genes,
     DEL2_genes, file = "TF_regulated_genes.Rdata")

x_gene_descri <- org.At.tairGENENAME
xx_gene_descri <- as.list(x_gene_descri[mappedkeys(x_gene_descri)])

x_gene_symb <- org.At.tairSYMBOL
xx_gene_symb <- as.list(x_gene_symb[mappedkeys(x_gene_symb)])


HSF_GO_terms_IDs = rownames(data.frame(HSF.BP))
HSF_GO_terms_Descriptions = data.frame(HSF.BP)[,'Description']
all.print_ready = c(paste('tair name','gene symbol','gene description','GO terms (HSF)',sep = "\t"))
for (mytair in HSF_genes){
  myGO_terms = ''
  for (myHSF_GO in HSF_GO_terms_IDs){
    genes_in_term = data.frame(HSF.BP)[myHSF_GO,'geneID']
    if (grepl(mytair, genes_in_term)){
      myterm = paste(myHSF_GO,":",data.frame(HSF.BP)[myHSF_GO,'Description'],'|',sep = "")
      myGO_terms = paste(myGO_terms, myterm, sep = '')
    }
  }
  mygene_symbol = paste(unlist(xx_gene_symb[mytair]), collapse = ' ')
  mygene_annota = xx_gene_descri[mytair][[1]][1]
  myGO_terms_print = substr(myGO_terms, 1, nchar(myGO_terms) - 1 )
  print_ready = paste(mytair,mygene_symbol,mygene_annota,myGO_terms_print,sep =' \t')
  all.print_ready = c(all.print_ready, print_ready)
}


fileConn<-file("HSF_genes.txt")
writeLines(all.print_ready, fileConn)
close(fileConn)

CAMTA1_GO_terms_IDs = rownames(data.frame(CAMTA1_BP))
CAMTA1_GO_terms_Descriptions = data.frame(CAMTA1_BP)[,'Description']
all.print_ready = c(paste('tair name','gene symbol','gene description','GO terms (HSF)',sep = "\t"))
for (mytair in CAMTA1_genes){
  myGO_terms = ''
  for (myCAMTA1_GO in CAMTA1_GO_terms_IDs){
    genes_in_term = data.frame(CAMTA1_BP)[myCAMTA1_GO,'geneID']
    if (grepl(mytair, genes_in_term)){
      myterm = paste(myCAMTA1_GO,":",data.frame(CAMTA1_BP)[myCAMTA1_GO,'Description'],'|',sep = "")
      myGO_terms = paste(myGO_terms, myterm, sep = '')
    }
  }
  mygene_symbol = paste(unlist(xx_gene_symb[mytair]), collapse = ' ')
  mygene_annota = xx_gene_descri[mytair][[1]][1]
  myGO_terms_print = substr(myGO_terms, 1, nchar(myGO_terms) - 1 )
  print_ready = paste(mytair,mygene_symbol,mygene_annota,myGO_terms_print,sep =' \t')
  all.print_ready = c(all.print_ready, print_ready)
}


fileConn<-file("CAMTA1_genes.txt")
writeLines(all.print_ready, fileConn)
close(fileConn)

