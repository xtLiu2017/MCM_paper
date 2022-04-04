setwd("/Users/liux/Documents/MCM_paper/data")
load("EV_genes.Rdata")
load("select.highquality.genes.Rdata")
load("peak_info_MCM_spline.Rdata")
load("genes_to_fix.Rdata")
load("peak_info_fixed.Rdata")
genes_to_fix_info = peak_info_fixed[genes_to_be_fixed,] 
peak_info_after_fix = peak_info
peak_info_after_fix[genes_to_be_fixed,] = genes_to_fix_info
colnames(peak_info_after_fix) = c("pt1_fix",'pl1_fix','pt2_fix','pl2_fix')
length(intersect(select.highquality.genes,EV_mock_positive_genes))
length(intersect(EV_mock_positive_genes, rownames(peak_info) ))
overlap_genes = intersect(select.highquality.genes,EV_mock_positive_genes)
peak_level_ratio = peak_info_after_fix[select.highquality.genes,'pl1_fix'] / peak_info_after_fix[select.highquality.genes,'pl2_fix']
sort_peak_level_ratio = sort(peak_level_ratio, na.last = NA)

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
window1 = 50
No.overlap_50 = sliding_window_func(sorted_peak_level_ratio = sort_peak_level_ratio, 
                                    genes_Dais = EV_mock_positive_genes, window = window1)
df_sup8a = data.frame(index = 1:dim(No.overlap_50)[2], ratio = No.overlap_50[1,],
                      No = No.overlap_50[2,])
myat = c()
for (ii in 1/c(0.25,0.5,1,1.25,2.5)){
  myat = c(myat, which.min(abs(ii - No.overlap_50[1,])))
}
EV_gene_no. = length(intersect(EV_mock_positive_genes,names(sort_peak_level_ratio)))
MCM_gene_no. = length(sort_peak_level_ratio)
lower_0.95 = qhyper(p=0.025, m = EV_gene_no., n=MCM_gene_no. - EV_gene_no., k = window1, lower.tail = TRUE, log.p = FALSE)
upper_0.95 = qhyper(p=0.975, m = EV_gene_no., n=MCM_gene_no. - EV_gene_no., k = window1, lower.tail = TRUE, log.p = FALSE)
mycol <- rgb(142,229,238, max = 255, alpha = 125)


plot_sup8b = ggplot(df_sup8a, aes(x = index, y = No)) + geom_line() + theme_bw() + 
  xlab('peak 1 level / peak 2 level') +
  ylab('No. of genes induced in'~italic(Pto) ~ ' DC3000') + 
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        plot.margin = unit(c(15, 15, 6, 10), 'mm'), axis.title = element_text(size = 12),
        axis.text = element_text(size = 10), panel.grid = element_blank()) + 
  scale_x_continuous(breaks = myat, labels = 1/c(0.25,0.5,1,1.25,2.5)) + 
  geom_hline(yintercept=window1 * EV_gene_no. / MCM_gene_no. , linetype = "dashed", col = "blue") + 
  annotate("rect", xmin = -50, ymin = lower_0.95, xmax = 1600, ymax = upper_0.95, alpha=0.5, fill=mycol) +
  annotate("text", x = 1270, y = 16.5, label = 'expected number with \n 95% confidence interval',
           size = 3.6, col = "blue") + 
  annotate('text',x = 800, y = 10,  label = 'window size = 50 genes', size = 5)
#sup figure 8b done, showing the overlap genes changes with peak level ratio
#no clear trend

index_first_time_induce = unlist(EV_mock_positive_first_time)[overlap_genes]
index_max_time_induce_all = unlist(EV_mock_positive_max_time)
names(index_max_time_induce_all) = EV_mock_positive_genes
index_max_time_induce = index_max_time_induce_all[overlap_genes]
length(intersect(EV_mock_positive_genes, rownames(peak_info)))
length(EV_mock_positive_genes)
length(rownames(peak_info))
#80% of genes responding in DC3000 also respond in AvrRpt2

library(ggplot2)
library(gridExtra)
mydf =data.frame(parameter = c(peak_info_after_fix[overlap_genes,'pt1_fix'], peak_info_after_fix[overlap_genes,'pl1_fix'],
                               peak_info_after_fix[overlap_genes,'pt2_fix'], peak_info_after_fix[overlap_genes,'pl2_fix']))
mydf$gene = rep(overlap_genes, 4)
mydf$first_time_induced = rep(index_first_time_induce, 4) 
mydf$max_time_induced = rep(index_max_time_induce, 4)
mydf$type = c(rep('pt1',length(overlap_genes)), rep('pl1', length(overlap_genes)),
              rep('pt2',length(overlap_genes)), rep('pl2', length(overlap_genes)))

mydf$first_time_induced = as.factor(mydf$first_time_induced)
mydf$max_time_induced = as.factor(mydf$max_time_induced)
indices = c(1,2,3,4,5,6,7,8,9,10)
mytime = c(1,2,3,4,6,9,12,16,20,24)
myfacet_label = c('peak 1 time', 'peak 1 level', 'peak 2 time', 'peak 2 level')
names(myfacet_label) = c('pt1', 'pl1','pt2', 'pl2')
sup.fig8a = ggplot(mydf, aes(x = max_time_induced, y = parameter)) + geom_boxplot() + theme_bw() + 
  xlab('Time of maximum induction in' ~ italic('Pto') ~ ' DC3000 (hour)')  + ylab('peak parameters') + scale_x_discrete(breaks = indices, labels = mytime) + 
  facet_wrap( ~ type, scales = 'free_y',labeller = labeller(type = myfacet_label)) + 
  theme(legend.text = element_text(size=10),legend.title = element_blank(),
        strip.background = element_blank(),strip.text.x = element_text(size = 13),
        axis.title = element_text(size = 10), axis.text = element_text(size = 8))
lay <- rbind(c(1,1,1,2,2,2,2))
             
jpeg("figures/sup8.jpeg", height = 100, width = 250,res = 350, unit = "mm")
grid.arrange(sup.fig8a,plot_sup8b, nrow = 1, layout_matrix = lay)
grid.text('A', x = unit(6,'mm'),y = unit(96,'mm'),gp = gpar(fontface = 'bold', fontsize = 6))
grid.text('B', x = unit(115,'mm'),y = unit(96,'mm'),gp = gpar(fontface = 'bold',fontsize = 6))
dev.off()


