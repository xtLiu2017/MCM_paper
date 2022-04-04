# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate figure 3.6 and 3.7 for my PhD thesis
# It is the temporal enrichment of the 33 prioritized TFs

# directory of MCM project, /MCM_paper/data
setwd("~/Documents/MCM_paper/data")

# load the packages
library(ggplot2)    # plotting
library(gridExtra)  # grid.arrange for figure layout
library(grid)       # grid.arrange for figure layout
library(tidyverse) 
library(tidytext)
library(qvalue)     # Storey FDR
library(ComplexHeatmap) # heatmaps, same below
library(circlize)
library(RColorBrewer)

# ---Figure 3.6A starts
# The p values for the sliding windows for all TFs, based on both first and second peaks
load("pvalues_sliding_windows.Rdata")
padj_first = p.adjust(c(p_first_all,p_first_max), method = "BH")
padj_second = p.adjust(c(p_second_all,p_second_max), method = "BH")

# get -log10 q values
FDR_second = -log10(c(p_second_all,p_second_max)[which.min(abs(padj_second - 0.05))])
FDR_first = -log10(c(p_first_all,p_first_max)[which.min(abs(padj_first - 0.05))])
FDR_mean = round(mean(c(FDR_first,FDR_second)), digits = 1)

# for figure 3.6A, the barplots
df_first = data.frame(TF_labels = toupper(names(sort(p_first_max[prioritized_first]))), 
                      pval = -log10(sort(p_first_max[prioritized_first])), peak = 1)
df_second = data.frame(TF_labels = toupper(names(sort(p_second_max[prioritized_second]))), 
                       pval = -log10(sort(p_second_max[prioritized_second])), peak = 2)

# merge the data frame for first and second peaks
# I name it fig5a, it means figure 3.6A here
df_fig5a = rbind(df_first,df_second)
tibble_fig5a = as.tibble(df_fig5a)
myfacet_label_5a = c('first peak','second peak')
names(myfacet_label_5a) = c(1,2)

fig5a_plot = tibble_fig5a %>% mutate(peak=as.factor(peak),
                                     TF_labels=reorder_within(TF_labels,pval, peak))  %>% 
  ggplot(aes(x = TF_labels, y=pval)) + 
  geom_bar(stat="identity",width = 0.7) + 
  facet_wrap( ~ peak , scales = 'free',labeller = labeller(peak = myfacet_label_5a), nrow = 1) + 
  coord_flip() + theme_bw() + ylab(expression('-log'[10]*' (p-value)')) + xlab("") + 
  ylim(0,9) + 
  theme(panel.border = element_blank(),axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18), axis.title.x = element_text(size = 26),
        strip.text.x = element_text(size=28), strip.background = element_blank(), 
        plot.margin = unit(c(10,5,8,5),units = "mm")) + 
  scale_x_reordered() 
#---Figure 3.6A done


# ---Figure 3.7 start, logically it makes sense to code 3.7 first and then followed by 3.6B
# TF binding across a lot of genes (~20000 genes), we will take a subset of induced gene binding later
load("Annotation_mat.Rdata")
tfs = rownames(Annotation_mat)
tf_labels = unlist(lapply(tfs,function(x){strsplit(x,split = "_")[[1]][1]}))
names(tf_labels) = tfs
genes_union = union(prioritized_first, prioritized_second)
tf_union = c()
for (mytf in genes_union){
  tf_union = c(tf_union, names(which(tf_labels == mytf)))
}

# This function calculated the number of genes with binding site within the sliding windows
Dense_sliding_window_func = function(x, y, window, start = 1){
  #x -- peak time (genes), y -- binding site number, window -- local number of gene
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

# load the peak time, peak level information 
load("peak_info_all.Rdata")
# 1939 high quality (upregulated modelled genes)
load("select.highquality.genes.Rdata")
# genes with altered fitting
load("genes_to_fix.Rdata")
peak_info_update = peak_info
peak_info_update[genes_to_be_fixed,] = peak_info_fixed[genes_to_be_fixed,]

# first peak times, second peak times
first_peaks = peak_info_update[select.highquality.genes,1]
second_peaks = peak_info_update[select.highquality.genes,3]
ordered_first_peaks = sort(first_peaks)
ordered_second_peaks = sort(second_peaks)

# window size
window1 = 150
# the peak times of all sliding windows, based on the median peak time of the sliding window
time_mat = c()
# number of genes with at least one TF binding site in the sliding windows
gene_number_mat = c()
# p value of all sliding windows (hypergeometric distribution)
p_mat = c()
for (mytf in tf_union){
  # TF binding site of the genes with non-NA peak times 
  my_TF_binding_first = Annotation_mat[mytf,names(ordered_first_peaks)]
  my_TF_binding_second = Annotation_mat[mytf,names(ordered_second_peaks)]
  
  # TF binding site of the 1939 genes 
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
  
  # p values for the first peaks
  p_first = 1 - phyper(mydensity_first - 1, k = window1, m = m, n = n, lower.tail = TRUE, log.p = FALSE)
  max_p_first = 1 - phyper(max(mydensity_first) - 1, k = window1, m = m, n = n, lower.tail = TRUE, log.p = FALSE)
  max_p_first_pos = which.max(mydensity_first)
  max_time_first = time_first[max_p_first_pos] + 3
  
  # p values for the second peaks
  p_second = 1 - phyper(mydensity_second - 1, k = window1, m = m, n = n, lower.tail = TRUE, log.p = FALSE)
  max_p_second = 1 - phyper(max(mydensity_second) - 1, k = window1, m = m, n = n, lower.tail = TRUE, log.p = FALSE)
  max_p_second_pos = which.max(mydensity_second)
  max_time_second = time_second[max_p_second_pos] + 3
  
  time_mat = rbind(time_mat, c(max_time_first,max_time_second))
  p_mat = rbind(p_mat, c( p_first,  p_second) )
  gene_number_mat = rbind(gene_number_mat, c(mydensity_first, mydensity_second))
}

# color the row names(TFs) based on TF family
genes_union_upp = toupper(genes_union)
rownames(gene_number_mat) = genes_union_upp
rownames(p_mat) = genes_union_upp
rowname_col = c()
for (i in genes_union_upp){
  if (grepl('HSF', i)){
    myrowname = 'red'
  }else if(grepl('ANAC',i) | grepl('NAP', i) | grepl("SND3", i)){
    myrowname = 'blue'
  }else if(grepl('WRKY',i)){
    myrowname = 'orange'
  }else{
    myrowname = 'gray30'
  }
  rowname_col = c(rowname_col, myrowname)
}

#window index - time transformation for sliding windows
myat_first = c()
for (ii in c(4,5,6,7,8)){
  myat_first = c(myat_first, which.min(abs(ii - time_first)))
}
myat_second = c()
for (ii in c(10,14,17,20,22)){
  myat_second = c(myat_second, which.min(abs(ii - time_second)))
}

#heatmap color scale
col = colorRamp2(seq(0,9, length.out = 8), brewer.pal(9,'Oranges')[1:8])
myat_first = c()
for (ii in c(4,5,6,7,8)){
  myat_first = c(myat_first, which.min(abs(ii - time_first)))
}
myat_second = c()
for (ii in c(10,14,17,20,22)){
  myat_second = c(myat_second, which.min(abs(ii - time_second)))
}

# get the labels for columns, get their times, for plotting
myalpha = rep(0, dim(p_mat)[2])
myalpha[myat_first] = 1
myalpha[myat_second + length(time_first)] = 1
mylabels = rep('', dim(p_mat)[2])
mylabels[myat_first] = c('4','5','6','7','8')
mylabels[myat_second + length(time_first)] = c('10','14','17','20','22')

# annotations showing the time scale for the heatmap columns
ha = HeatmapAnnotation(time = anno_lines(c(0.001,rep(0,dim(p_mat)[2]-1)), ylim = c(-1,1),
                                         border = F, axis = F, gp = gpar(lwd = 2),add_points = T,
                                         pt_gp = gpar(alpha = myalpha)),
                       foo = anno_text(mylabels, rot = 0, just = 'center',gp = gpar(fontsize = 12)),
                       show_annotation_name=F)
#clustering the rows
dist_mat = dist(-log10(p_mat))
mhclust = hclust(dist_mat,method = "average")  

#legend 
lgd = Legend(col_fun = col, title = expression('-log'[10]*' (p-value)'),
             at = c(0, 2, 5, 7, 9, FDR_mean), labels = c('0','2','5','7','9',paste(FDR_mean,'(FDR = 0.05)')),
             legend_height = unit(40,"mm"),legend_width = unit(30,"mm"),title_position = 'topleft',
             title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 14),title_gap = unit(4,"mm"))
ht1 = Heatmap(-log10(p_mat), col = col,cluster_columns = FALSE, 
              heatmap_height = unit(15, "cm"),
              width = unit(12,'cm'),
              cluster_rows = mhclust,
              name = 'j', row_dend_gp = gpar(), show_row_names = T, show_column_names = F,
              show_heatmap_legend = F, row_title = paste(dim(p_mat)[1],"TFs"), row_title_rot = 0, 
              use_raster = F, column_title_gp = gpar(fontsize = 18),
              row_title_gp = gpar(fontsize = 16),
              column_split = c(rep("first peak", times = length(p_first)),rep("second peak", times = length(p_second))),
              column_gap = unit(5,"mm"),bottom_annotation = ha,
              row_names_gp = gpar(col = rowname_col,fontsize = 12))
# grab the heatmap
grob.fig5c = grid.grabExpr(draw(ht1,heatmap_legend_list = lgd)) 
# change the output directory as needed
jpeg('~/Desktop/paper_figures/thesis/figure3.7.jpeg', width = 230, height = 175, unit = 'mm', res = 400)
draw(ht1,heatmap_legend_list = lgd)
grid.text('time (hour)', x = unit(172,'mm'),y = unit(16,'mm'), gp = gpar(fontsize = 14))
dev.off()


#---Figure 3.7 done

#---Figure 3.6B start

#obtain data matrix for three TF families
HSFs = c('HSFC1_col_a','HSF3_col_a','HSFA6B_col','HSF6_col_a','HSF7_col_a')
ANACs = c('ANAC047_col','ANAC083_col_v3a','ANAC058_col_a','ANAC062_col','ANAC055_col','ANAC038_col_a')
WRKYs = c('WRKY7_col','WRKY75_col_a','WRKY15_col_b','WRKY70_col','WRKY24_col_a')
HSF_mat_1 = c()
HSF_mat_2 = c()
ANAC_mat_1 = c()
ANAC_mat_2 = c()
WRKY_mat_1 = c()
WRKY_mat_2 = c()
for (mytf in c(HSFs, ANACs, WRKYs)){
  my_TF_binding_first = Annotation_mat[mytf,names(ordered_first_peaks)]
  my_TF_binding_second = Annotation_mat[mytf,names(ordered_second_peaks)]
  
  time_first = Dense_sliding_window_func(x = ordered_first_peaks, y = my_TF_binding_first,window = window1,start = 1)[1,] + 3
  mydensity_first = Dense_sliding_window_func(x = ordered_first_peaks, y = my_TF_binding_first,window = window1,start = 1)[2,] 
  time_second = Dense_sliding_window_func(x = ordered_second_peaks, y = my_TF_binding_second,window = window1,start = 1)[1,] + 3
  mydensity_second = Dense_sliding_window_func(x = ordered_second_peaks, y = my_TF_binding_second,window = window1,start = 1)[2,]
  if (mytf %in% HSFs){
    HSF_mat_1 = rbind(HSF_mat_1, mydensity_first)
    HSF_mat_2 = rbind(HSF_mat_2, mydensity_second)
  }else if(mytf %in% ANACs){
    ANAC_mat_1 = rbind(ANAC_mat_1, mydensity_first)
    ANAC_mat_2 = rbind(ANAC_mat_2, mydensity_second)
  }else{
    WRKY_mat_1 = rbind(WRKY_mat_1, mydensity_first)
    WRKY_mat_2 = rbind(WRKY_mat_2, mydensity_second)
  }
}
rownames(HSF_mat_1) = HSFs
rownames(HSF_mat_2) = HSFs
rownames(ANAC_mat_1) = ANACs
rownames(ANAC_mat_2) = ANACs
rownames(WRKY_mat_1) = WRKYs
rownames(WRKY_mat_2) = WRKYs

#get the time scale
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

# This is a trick
# this is for using different time scales for the first and second peaks
# under the ggplot geom_facet() function
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

#prepare HSF data
HSF_mat_1_t = t(HSF_mat_1)
HSF_mat_2_t = t(HSF_mat_2)
HSF_df = data.frame(no = c(as.vector(HSF_mat_1_t), as.vector(HSF_mat_2_t)))
rownames(HSF_df) = NULL
HSF_df$time = as.vector(c(rep(time_first, dim(HSF_mat_1)[1]), rep(time_second, dim(HSF_mat_2)[1])))
HSF_df$peak = c(rep(1, length(time_first) * dim(HSF_mat_1)[1]), rep(2, length(time_second) * dim(HSF_mat_2)[1]))
HSF_df$gene = c(rep(HSFs, each = length(time_first)),rep(HSFs, each = length(time_second)))
HSF_df$index = c(rep(seq(1, length(time_first)), dim(HSF_mat_1)[1]), rep(seq(1, length(time_second)), dim(HSF_mat_2)[1]))
HSF_label = unlist(lapply(HSFs,function(x){strsplit(x,split = "_")[[1]][1]}))
names(HSF_label) = HSFs
names(myfacet_label) = c(1,2)

#colors for HSF TFs
col_HSF = brewer.pal(9,'Reds')[3:7]

#prepare ANAC data
ANAC_mat_1_t = t(ANAC_mat_1)
ANAC_mat_2_t = t(ANAC_mat_2)
ANAC_df = data.frame(no = c(as.vector(ANAC_mat_1_t), as.vector(ANAC_mat_2_t)))
rownames(ANAC_df) = NULL
ANAC_df$time = as.vector(c(rep(time_first, dim(ANAC_mat_1)[1]), rep(time_second, dim(ANAC_mat_2)[1])))
ANAC_df$peak = c(rep(1, length(time_first) * dim(ANAC_mat_1)[1]), rep(2, length(time_second) * dim(ANAC_mat_2)[1]))
ANAC_df$gene = c(rep(ANACs, each = length(time_first)),rep(ANACs, each = length(time_second)))
ANAC_df$index = c(rep(seq(1, length(time_first)), dim(ANAC_mat_1)[1]), rep(seq(1, length(time_second)), dim(ANAC_mat_2)[1]))
#colors for ANAC TFs
col_ANAC = brewer.pal(9,'Blues')[3:8]

ANAC_label = unlist(lapply(ANACs,function(x){strsplit(x,split = "_")[[1]][1]}))
names(ANAC_label) = ANACs
names(myfacet_label) = c(1,2)

#prepare WRKY data
WRKY_mat_1_t = t(WRKY_mat_1)
WRKY_mat_2_t = t(WRKY_mat_2)
WRKY_df = data.frame(no = c(as.vector(WRKY_mat_1_t), as.vector(WRKY_mat_2_t)))
rownames(WRKY_df) = NULL
WRKY_df$time = as.vector(c(rep(time_first, dim(WRKY_mat_1)[1]), rep(time_second, dim(WRKY_mat_2)[1])))
WRKY_df$peak = c(rep(1, length(time_first) * dim(WRKY_mat_1)[1]), rep(2, length(time_second) * dim(WRKY_mat_2)[1]))
WRKY_df$gene = c(rep(WRKYs, each = length(time_first)),rep(WRKYs, each = length(time_second)))
WRKY_df$index = c(rep(seq(1, length(time_first)), dim(WRKY_mat_1)[1]), rep(seq(1, length(time_second)), dim(WRKY_mat_2)[1]))
#colors for WRKY TFs
col_WRKY = brewer.pal(9,'Oranges')[3:7]

WRKY_label = unlist(lapply(WRKYs,function(x){strsplit(x,split = "_")[[1]][1]}))
names(WRKY_label) = WRKYs
names(myfacet_label) = c(1,2)


HSF_df$family = 'HSF'
ANAC_df$family = 'ANAC'
WRKY_df$family = 'WRKY'
# This is a data frame for all three families
family_df = rbind(HSF_df,ANAC_df,WRKY_df)
family_df$family = factor(family_df$family, levels = c("HSF",'WRKY','ANAC'))
#fig5B actually means figure 3.6B
fig5B = ggplot(family_df, aes(x = index, y = no, color = gene)) + geom_line() + theme_bw() +
  facet_wrap( family ~ peak, scales = 'free_x',labeller = labeller(peak = myfacet_label), ncol = 2) +
  xlab("time (hour)") + ylab("No. of genes with at least 1 binding site") + 
  scale_x_continuous(breaks=breaks_fun, labels = label_func) +
  scale_color_manual(breaks = c(HSFs,WRKYs,ANACs), values=c(col_HSF,col_WRKY, col_ANAC), label = c(HSF_label, ANAC_label, WRKY_label)) +
  theme(legend.text = element_text(size=14),legend.title = element_blank(),
        strip.background = element_blank(),strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 15),
        plot.margin = unit(c(3,8,3,8), 'mm'))

lay <- rbind(c(1,2),c(1,2))

#change output directory as needed
jpeg('~/Desktop/paper_figures/thesis/figure3.6A-B.jpeg', width = 450, height = 275, unit = 'mm', res = 400)
grid.arrange(fig5a_plot,fig5B,layout_matrix = lay)
grid.text('A', x = unit(7,'mm'),y = unit(270,'mm'),gp = gpar(fontface = 'bold', fontsize = 18))
grid.text('B', x = unit(225,'mm'),y = unit(270,'mm'),gp = gpar(fontface = 'bold',fontsize = 18))
dev.off()

