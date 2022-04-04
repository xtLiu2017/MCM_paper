# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate figure 4 for the manuscript

# directory of MCM project, /MCM_paper/data
setwd("~/Documents/MCM_paper/data")

# packages
library(ggplot2)      # plotting
library(gridExtra)    # grid.arrange, for figure panel layout
library(grid)         # grid.arrange
library(tidyverse)
library(tidytext)
library(clusterProfiler)  # GO enrichment for single peak genes
library(org.At.tair.db)   # Arabidopsis database

# import Daisuke's data
# The Daisuke's data was processed by Dr. Fumi Katagiri and it will be/has been published in another paper
load("peak.time.level.EdAvrRpt2.RData")
genes_Dais = names(peak.time)
# import MCM peak info
load("peak_info_MCM_spline.Rdata")
load("select.highquality.genes.Rdata")

# intersect genes
genes_intersect = intersect(genes_Dais,select.highquality.genes)

# remove the problematic gene !!!!!!!
# This gene has the smallest peak time in Daisuke's data but the largest peak time in MCM data
# The fitting looks good. So this is so weired.
genes_intersect = setdiff(genes_intersect, 'AT5G37540')
colnames(peak_info) = c('pt1_i1','pl1_i1','pt2_i1','pl2_i1')

# Daisuke genes with reasonable and unreasonable fit, if greater than 7.5 hours, probably bad fitting
# because no time points after 6 hours
# < 7.5 hours called genes_good, others called genes_bad
peak_time_Dais = peak.time[genes_intersect]
genes_good = genes_intersect[which(peak_time_Dais < 7.5)]
genes_bad = genes_intersect[which(peak_time_Dais > 7.5)]

# the genes with peak + shoulder, fit with altered MCMs
load("genes_to_fix.Rdata")
# the synthetic profiles of the genes fit with altered MCMs
load("peak_info_fixed.Rdata")
genes_to_fix_here = intersect(genes_to_be_fixed,genes_intersect)
genes_to_fix_info = peak_info_fixed[genes_to_fix_here,] 
peak_info_after_fix = peak_info[genes_intersect,]
peak_info_after_fix[genes_to_fix_here,] = genes_to_fix_info
colnames(peak_info_after_fix) = c("pt1_fix",'pl1_fix','pt2_fix','pl2_fix')

# data frame containing all infor for intersect genes
df = data.frame(peak_info_after_fix,
                pl_Dais = peak.level[genes_intersect], pt_Dais = peak.time[genes_intersect])
df_good = df[genes_good,]
df_bad = df[genes_bad,]

#Below are figure 4A-D
# peak time 1 of MCM vs peak time of Daisuke
pcc_1 = cor(df_good["pt1_fix"],df_good["pt_Dais"],use = "pair")
scc_1 = cor(df_good["pt1_fix"],df_good["pt_Dais"],use = "pair",method = "spearman")
p1 = ggplot(df_good,aes(y = pt1_fix + 3, x= pt_Dais)) + geom_bin2d(bins = 100) + theme_bw() + 
  labs(y = "first peak time (MCM, set1)",x = substitute(paste('peak time', italic(' (In planta)'))),size = 20)+ 
  annotate("text", x = 4.2, y=9.2, label = paste("PCC =",round(pcc_1,2)),size = 4)+
  annotate("text", x = 4.2, y=8.5, label = paste("SCC =",round(scc_1,2)),size = 4)+
  scale_fill_distiller(palette = "Spectral",direction = -1)+
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=13),
        axis.text = element_text(size = 11),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        plot.margin = unit(c(5,5,5,5), 'mm')) + 
  geom_abline(slope=1,intercept = 0,size=0.5,colour = "grey") + ggtitle("(first) peak time") +
  geom_boxplot(data = df_bad, aes(x = 8.5), color = "#3288BD",outlier.size=0.3) + 
  scale_x_continuous(breaks = c(4,5,6,7,8.5), labels = c('4','5','6','7','>7.5'), limits = c(3,9))
  
# peak level 1 of MCM vs peak level of Daisuke
pcc_2 = cor(df_good["pl1_fix"],df_good["pl_Dais"],use = "pair")
scc_2 = cor(df_good["pl1_fix"],df_good["pl_Dais"],use = "pair",method = "spearman")
p2 = ggplot(df,aes(y = pl1_fix, x= pl_Dais)) + geom_bin2d(bins = 100) + theme_bw() + 
  labs(y = "first peak level (MCM, set1)",x = substitute(paste('peak level', italic(' (In planta)'))),size = 20) + 
  annotate("text", x = 2.5, y=5.55, label = paste("PCC =",round(pcc_2,2)),size = 4)+
  annotate("text", x = 2.5, y=5, label = paste("SCC =",round(scc_2,2)),size = 4)+
  scale_fill_distiller(palette = "Spectral",direction = -1)+
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=13),
        axis.text = element_text(size = 11),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        plot.margin = unit(c(5,5,5,5), 'mm')) + ggtitle("(first) peak level")
  
# peak time 2 of MCM vs peak time of Daisuke
pcc_3 = cor(df_good["pt2_fix"],df_good["pt_Dais"],use = "pair")
scc_3 = cor(df_good["pt2_fix"],df_good["pt_Dais"],use = "pair",method = "spearman")
p3 = ggplot(df_good,aes(y = pt2_fix + 3, x= pt_Dais)) + geom_bin2d(bins = 100) + theme_bw() + 
  labs(y = "second peak time (MCM, set1)",x = substitute(paste('peak time', italic(' (In planta)'))),size = 20)+ 
  annotate("text", x = 4.3, y=23, label = paste("PCC =",round(pcc_3,2)),size = 4)+
  annotate("text", x = 4.3, y=21.5, label = paste("SCC =",round(scc_3,2)),size = 4)+
  scale_fill_distiller(palette = "Spectral",direction = -1, limits = c(1,7))+
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=13),
        axis.text = element_text(size = 11),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        plot.margin = unit(c(5,5,5,5), 'mm')) +
  ggtitle("(second) peak time") +
  geom_boxplot(data = df_bad, aes(x = 8.5), color = "#3288BD",outlier.size=0.3) + 
  scale_x_continuous(breaks = c(4,5,6,7,8.5), labels = c('4','5','6','7','>7.5'), limits = c(3,9))

# peak level 2 of MCM vs peak level of Daisuke
pcc_4 = cor(df_good["pl2_fix"],df_good["pl_Dais"],use = "pair")
scc_4 = cor(df_good["pl2_fix"],df_good["pl_Dais"],use = "pair",method = "spearman")
p4 = ggplot(df_good,aes(y = pl2_fix, x= pt_Dais)) + geom_bin2d(bins = 100) + theme_bw() + 
  labs(y = "second peak level (MCM, set1)",x = substitute(paste('peak level', italic(' (In planta)'))),size = 20)+ 
  annotate("text", x = 4.4, y=3.2, label = paste("PCC =",round(pcc_4,2)),size = 4)+
  annotate("text", x = 4.4, y=2.87, label = paste("SCC =",round(scc_4,2)),size = 4)+
  scale_fill_distiller(palette = "Spectral",direction = -1, limits = c(1,7))+
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=13),
        axis.text = element_text(size = 11),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        plot.margin = unit(c(5,5,5,5), 'mm')) + ggtitle("(second) peak level")

# look at overlap of genes with Daisuke's data (expressed in cell pop 1)
# This function calculated the number of genes in Daisuke's upregulated gene set within the sliding windows
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
# 1939 high quality (upregulated modelled genes)
peak_info_1939 = peak_info[select.highquality.genes,]
# genes with altered fitting
peak_info_1939[genes_to_be_fixed,] = peak_info_fixed[genes_to_be_fixed,]
peak_level_ratio = peak_info_1939[,'pl1_i1'] / peak_info_1939[,'pl2_i1']
sort_peak_level_ratio = sort(peak_level_ratio, na.last = NA)
window1 = 30  # window size
No.overlap_30 = sliding_window_func(sorted_peak_level_ratio = sort_peak_level_ratio, 
                                    genes_Dais = genes_Dais, window = window1)

df_fig3b = data.frame(index = 1:dim(No.overlap_30)[2], ratio = No.overlap_30[1,], No = No.overlap_30[2,])
# window index - peak level ratio transformation for sliding windows
myat = c()
for (ii in 1/c(0.25,0.5,1,1.25,2.5)){
  myat = c(myat, which.min(abs(ii - No.overlap_30[1,])))
}

# some background information, for hypergeometric test
Dais_gene_no. = length(intersect(genes_Dais,names(sort_peak_level_ratio)))
MCM_gene_no. = length(sort_peak_level_ratio)
lower_0.95 = qhyper(p=0.025, m = Dais_gene_no., n=MCM_gene_no. - Dais_gene_no., k = window1, lower.tail = TRUE, log.p = FALSE)
upper_0.95 = qhyper(p=0.975, m = Dais_gene_no., n=MCM_gene_no. - Dais_gene_no., k = window1, lower.tail = TRUE, log.p = FALSE)
mycol <- rgb(142,229,238, max = 255, alpha = 125)

# the figure of temporal 'enrichment', this is figure 4E (fig3b is a misleading notation below)
fig3b = ggplot(df_fig3b, aes(x = index, y = No)) + geom_line() + theme_bw() + xlab('peak 1 level / peak 2 level') +
  ylab(bquote(atop('No. of genes induced', 'by ' ~ italic('In planta') ~ ' AvrRpt2'))) + 
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), panel.grid = element_blank(),
        plot.margin = unit(c(8,8,8,8), 'mm')) + 
  scale_x_continuous(breaks = myat, labels = 1/c(0.25,0.5,1,1.25,2.5)) + 
  geom_hline(yintercept=window1 * Dais_gene_no. / MCM_gene_no. , linetype = "dashed", col = "blue") + 
  annotate("rect", xmin = -50, ymin = lower_0.95, xmax = 1600, ymax = upper_0.95, alpha=0.5, fill=mycol) +
  annotate("text", x = 1270, y = 10, label = 'expected number with \n 95% confidence interval',
           size = 3.6, col = "blue") + 
  annotate('text',x = 900, y = 3.5,  label = 'window size = 30 genes', size = 5)


sp_genes = names(which(peak_level_ratio > 4))
#for single peak genes
mean(peak_info_1939[sp_genes, 'pt1_i1']) + 3
sd(peak_info_1939[sp_genes, 'pt1_i1'])
#for all genes in the background
mean(peak_info_1939[, 'pt1_i1'],na.rm = T) + 3
sd(peak_info_1939[, 'pt1_i1'],na.rm = T) 
#GO analysis
sp.BP <- enrichGO(gene = sp_genes,universe = select.highquality.genes,
                  OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                  pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = F)

#extract some top significant terms
top12_terms_sp = data.frame(sp.BP)[1:12,1]
df_GO_HSF = data.frame(terms_ID = c(top12_terms_sp), 
                       terms = data.frame(sp.BP)[top12_terms_sp,'Description'],
                       p_val = data.frame(sp.BP)[top12_terms_sp,'pvalue'])

#plot these terms along with their -log p values
ggplot_sp = ggplot(df_GO_HSF,aes(x = terms, y = -log(p_val))) + 
  geom_bar(stat="identity",position=position_dodge(),width = 0.7) +
  coord_flip() + theme_bw() + xlab('GO terms') + ylab('-log(p-value)') + 
  ggtitle("GO enrichment of single peak genes") + 
  theme(panel.border = element_blank(), axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),plot.title = element_text(size=14,hjust=1.6,face = 'bold'),
        plot.margin = unit(c(2,2,2,2),'mm')) + 
  scale_x_discrete(limits=rev(df_GO_HSF$terms[1:12])) 

#layout matrix indicating how to organize the heatmaps in a page
lay <- rbind(c(1,1,2,2,3,3,3),
             c(4,4,5,5,6,6,6))
             
#change your directory as needed
jpeg('~/Documents/MCM_paper/manuscript_figure/Figure4.jpeg', width = 320, height = 190, unit = 'mm', res = 400)
grid.arrange(p1,p2,fig3b,p3,p4,ggplot_sp,layout_matrix = lay)
grid.text('A', x = unit(5,'mm'),y = unit(185,'mm'),gp = gpar(fontface = 'bold'))
grid.text('B', x = unit(97,'mm'),y = unit(185,'mm'),gp = gpar(fontface = 'bold'))
grid.text('C', x = unit(5,'mm'),y = unit(93,'mm'),gp = gpar(fontface = 'bold'))
grid.text('D', x = unit(97,'mm'),y = unit(93,'mm'),gp = gpar(fontface = 'bold'))
grid.text('E', x = unit(183,'mm'),y = unit(185,'mm'),gp = gpar(fontface = 'bold'))
grid.text('F', x = unit(183,'mm'),y = unit(93,'mm'),gp = gpar(fontface = 'bold'))
dev.off()
