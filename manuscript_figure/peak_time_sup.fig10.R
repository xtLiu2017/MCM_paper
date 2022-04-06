# Author: Xiaotong Liu
# plot peak time comparison
# supplemental figure 10

# directory where the data are saved
setwd("~/Documents/MCM_paper/data")
# 1939 high quality genes
load("select.highquality.genes.Rdata")
# peak time, peak level of 3039 genes
load("peak_info_all.Rdata")
# genes fit with altered MCMs
load("genes_to_fix.Rdata")
# change the peak times and peak levels for those genes fit with altered MCMs
peak_info_update = peak_info
peak_info_update[genes_to_be_fixed,] = peak_info_fixed[genes_to_be_fixed,]

# first peak time of spline, MCM input 1, MCM input 2
peak_1_time_spline = spline_peak_info[select.highquality.genes,1] + 3
peak_1_time_MCM = peak_info_update[select.highquality.genes,1] + 3
peak_1_time_MCM_set2 = peak_info_set2[select.highquality.genes,1] + 3
# second peak time of spline, MCM input 1, MCM input 2
peak_2_time_spline = spline_peak_info[select.highquality.genes,3] + 3
peak_2_time_MCM = peak_info_update[select.highquality.genes,3] + 3
peak_2_time_MCM_set2 = peak_info_set2[select.highquality.genes,3] + 3
# pair-wise pcc
pcc_1_1 = cor(peak_1_time_spline,peak_1_time_MCM,use = "pair")
pcc_2_1 = cor(peak_2_time_spline,peak_2_time_MCM,use = "pair")
pcc_1_2 = cor(peak_1_time_spline,peak_1_time_MCM_set2,use = "pair")
pcc_2_2 = cor(peak_2_time_spline,peak_2_time_MCM_set2,use = "pair")
# pair-wise scc
scc_1_1 = cor(peak_1_time_spline,peak_1_time_MCM,use = "pair",method = "spearman")
scc_2_1 = cor(peak_2_time_spline,peak_2_time_MCM,use = "pair",method = "spearman")
scc_1_2 = cor(peak_1_time_spline,peak_1_time_MCM_set2,use = "pair",method = "spearman")
scc_2_2 = cor(peak_2_time_spline,peak_2_time_MCM_set2,use = "pair",method = "spearman")

# packages
library(tidyverse)
library(ggplot2)         # plotting 
library(RColorBrewer)    # colors
library(gridExtra)
color1 = brewer.pal(12, 'Paired')[1]
color2 = brewer.pal(12, "Paired")[2]
# plot MCM vs spline
# input set1
peak1_MCM_spline = data.frame(spline = peak_1_time_spline,MCM = peak_1_time_MCM)
peak2_MCM_spline = data.frame(spline = peak_2_time_spline,MCM = peak_2_time_MCM)
# input set2
peak1_MCM_spline_set2 = data.frame(spline = peak_1_time_spline,MCM = peak_1_time_MCM_set2)
peak2_MCM_spline_set2 = data.frame(spline = peak_2_time_spline,MCM = peak_2_time_MCM_set2)

p1 = ggplot(peak1_MCM_spline,aes(x = spline, y= MCM)) + geom_bin2d(bins = 150) + theme_bw() + 
  labs(x = "",y = "MCMs (input set1)")+ 
  annotate("text", x = 13, y=6, label = paste("PCC =",round(pcc_1_1,2)),size = 4.5)+
  annotate("text", x = 13, y=4.5, label = paste("SCC =",round(scc_1_1,2)),size = 4.5)+
  scale_fill_distiller(palette = "Spectral",direction = -1) + xlim(3,16) + ylim(3,16) +
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=14),
        axis.text = element_text(size = 12),legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.25, 'cm'), plot.title = element_text(hjust = 0.5, size = 16,face = 'bold')) +
  ggtitle('first peak time') + 
  geom_abline(slope=1,intercept = 0,size=0.5,colour = "grey")

p2 = ggplot(peak2_MCM_spline,aes(x = spline, y= MCM)) + geom_bin2d(bins = 150) + theme_bw() + 
  labs(x = "",y = "MCMs (input set1)") + xlim(5,25) + ylim(5,25) + 
  annotate("text", x = 21, y=9, label = paste("PCC =",round(pcc_2_1,2)),size = 4.5)+
  annotate("text", x = 21, y=7, label = paste("SCC =",round(scc_2_1,2)),size = 4.5)+
  scale_fill_distiller(palette = "Spectral",direction = -1)+ ggtitle('second peak time') + 
  theme(axis.line = element_line(colour = "black"),axis.title = element_text(size=14),
        axis.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        panel.border = element_blank(),panel.background = element_blank(),legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.25, 'cm')) +geom_abline(slope=1,intercept = 0,size=0.5,colour = "grey")

p3 = ggplot(peak1_MCM_spline_set2,aes(x = spline, y= MCM)) + geom_bin2d(bins = 150) + theme_bw() + 
  labs(x = "spline-based",y = "MCMs (input set2)")+ 
  annotate("text", x = 13, y=6, label = paste("PCC =",round(pcc_1_2,2)),size = 4.5)+
  annotate("text", x = 13, y=4.5, label = paste("SCC =",round(scc_1_2,2)),size = 4.5)+
  scale_fill_distiller(palette = "Spectral",direction = -1) + xlim(3,16) + ylim(3,16) +
  ggtitle('first peak time') + 
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=14),
        axis.text = element_text(size = 12),legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.25, 'cm'), plot.title = element_text(hjust = 0.5, size = 16, face = 'bold')) +
  geom_abline(slope=1,intercept = 0,size=0.5,colour = "grey")

p4 = ggplot(peak2_MCM_spline_set2,aes(x = spline, y= MCM)) + geom_bin2d(bins = 150) + theme_bw() + 
  labs(x = "spline-based",y = "MCMs (input set2)") + xlim(5,25) + ylim(5,25) +
  annotate("text", x = 21, y=9, label = paste("PCC =",round(pcc_2_2,2)),size = 4.5)+
  annotate("text", x = 21, y=7, label = paste("SCC =",round(scc_2_2,2)),size = 4.5)+
  scale_fill_distiller(palette = "Spectral",direction = -1) + ggtitle('second peak time') + 
  theme(axis.line = element_line(colour = "black"),axis.title = element_text(size=14),
        axis.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'), 
        panel.border = element_blank(),panel.background = element_blank(),
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.25, 'cm')) + geom_abline(slope=1,intercept = 0,size=0.5,colour = "grey")
jpeg("~/Documents/MCM_paper/manuscript_figure/sup.fig10a-d.jpeg",width = 180, height = 180, unit = "mm", res = 400)
grid.arrange(p1,p2,p3,p4, nrow = 2)
dev.off()



#plot first peak vs second peak (MCM input set 1)
peak_MCM = data.frame(first = peak_1_time_MCM,second = peak_2_time_MCM, color1 = rep('1',length(peak_1_time_MCM)),
                      color2 = rep('2', length(peak_1_time_MCM)))
pcc_1_2_peak = cor(peak_1_time_MCM,peak_2_time_MCM,use = "pair")
scc_1_2_peak = cor(peak_1_time_MCM,peak_2_time_MCM,use = "pair",method = "spearman")
pcc_1_2_peak_log = cor(log(peak_1_time_MCM),log(peak_2_time_MCM),use = "pair")
scc_1_2_peak_log = cor(log(peak_1_time_MCM),log(peak_2_time_MCM),use = "pair",method = "spearman")

p5 = ggplot(peak_MCM ,aes(x = first, y= second)) + geom_bin2d(bins = 150) + theme_bw() + 
  labs(x = "first peak time (MCM, set1)",y = "second peak time (MCM, set1)") + 
  annotate("text", x = 13, y=12, label = paste("PCC =",round(pcc_1_2_peak,2)), size = 5)+
  annotate("text", x = 13, y=10, label = paste("SCC =",round(scc_1_2_peak,2)), size = 5)+
  scale_fill_distiller(palette = "Spectral",direction = -1) +
  theme(axis.line = element_line(colour = "black"),axis.title = element_text(size=14),axis.text = element_text(size = 13),
        panel.border = element_blank(),panel.background = element_blank()) +
  theme(plot.margin = unit(c(5,5,5,5), "mm"))


p6 = ggplot(peak_MCM ,aes(x = first, y= second)) + theme_bw() + 
  geom_point(aes(colour = color2), size = 0.2) + 
  geom_jitter(aes(colour = color1), position=position_jitter(1),size = 0.2) +
  labs(x = "first peak time (MCM, set1)",y = "second peak time (MCM, set1)", colour = "data type") + 
  theme(axis.line = element_line(colour = "black"),axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        panel.border = element_blank(),panel.background = element_blank(),
        legend.title = element_text(size = 15),legend.text = element_text(size = 12),
        legend.position = "top",axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text( vjust = 4)) +
  theme(plot.margin = unit(c(5,5,5,5), "mm")) + 
  scale_color_manual(breaks = c('1','2'), values = c(color1,color2), labels = c("jittered","raw")) +
  guides(colour = guide_legend(override.aes = list(size=3)))

p7 = ggplot(peak_MCM ,aes(x = log(first), y= log(second))) + geom_bin2d(bins = 150) + theme_bw() + 
  labs(x = "first peak time (log; MCM, set1)",y = "second peak time(log; MCM, set1)")+ 
  annotate("text", x = 2.4, y=2.45, label = paste("PCC =",round(pcc_1_2_peak_log,2)), size = 5)+
  annotate("text", x = 2.4, y=2.3, label = paste("SCC =",round(scc_1_2_peak_log,2)), size = 5)+
  scale_fill_distiller(palette = "Spectral",direction = -1)+
  theme(axis.line = element_line(colour = "black"),axis.title = element_text(size=14),
        axis.text = element_text(size = 13),
        panel.border = element_blank(),panel.background = element_blank()) +
  theme(plot.margin = unit(c(5,5,5,5), "mm"))

p8 = ggplot(peak_MCM ,aes(x = log(first), y= log(second))) + theme_bw() + 
  geom_point(aes(colour = color2), size = 0.2) + 
  geom_jitter(aes(colour = color1), position=position_jitter(0.3),size = 0.2) +
  labs(x = "first peak time (log; MCM, set1)",y = "second peak time (log; MCM, set1)", colour = "data type") + 
  theme(axis.line = element_line(colour = "black"),axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        panel.border = element_blank(),panel.background = element_blank(),
        legend.title = element_text(size = 15),legend.text = element_text(size = 12),
        legend.position = "top",axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text( vjust = 4)) +
  theme(plot.margin = unit(c(5,5,5,5), "mm")) + 
  scale_color_manual(breaks = c('1','2'), values = c(color1,color2), labels = c("jittered","raw")) +
  guides(colour = guide_legend(override.aes = list(size=3)))


jpeg("~/Documents/MCM_paper/manuscript_figure/sup.fig10e-h.jpeg",width = 180, height = 180, unit = "mm", res = 400)
grid.arrange(p5,p6,p7,p8, nrow = 2)
dev.off()
# end here


# I don't not make a figure for the code below
# plot first vs second using spline fitting
peak_spline = data.frame(first = peak_1_time_spline,second = peak_2_time_spline)
pcc_1_2_peak_spline = cor(peak_1_time_spline,peak_2_time_spline,use = "pair")
scc_1_2_peak_spline = cor(peak_1_time_spline,peak_2_time_spline,use = "pair",method = "spearman")
pcc_1_2_peak_log_spline = cor(log(peak_1_time_spline),log(peak_2_time_spline),use = "pair")
scc_1_2_peak_log_spline = cor(log(peak_1_time_spline),log(peak_2_time_spline),use = "pair",method = "spearman")

jpeg("first_second_pkt_spline.jpeg")
ggplot(peak_spline ,aes(x = first, y= second)) + geom_bin2d(bins = 150) + theme_bw() + 
  labs(x = "first peak time (spline)",y = "second peak time (spline)")+ 
  annotate("text", x = 12.5, y=10.5, label = paste("PCC =",round(pcc_1_2_peak_spline,3)), size = 7)+
  annotate("text", x = 12.5, y=8.7, label = paste("SCC =",round(scc_1_2_peak_spline,3)), size = 7)+
  scale_fill_distiller(palette = "Spectral",direction = -1)+
  theme(axis.line = element_line(colour = "black"),axis.title = element_text(size=18),axis.text = element_text(size = 15),
        panel.border = element_blank(),panel.background = element_blank()) 
dev.off()
jpeg("first_second_pkt_spline_log.jpeg")
ggplot(peak_spline ,aes(x = log(first), y= log(second))) + geom_bin2d(bins = 150) + theme_bw() + 
  labs(x = "first peak time (log; spline)",y = "second peak time (log; spline)")+ 
  annotate("text", x = 2.5, y=2.35, label = paste("PCC =",round(pcc_1_2_peak_log_spline,3)), size = 7)+
  annotate("text", x = 2.5, y=2.16, label = paste("SCC =",round(scc_1_2_peak_log_spline,3)), size = 7)+
  scale_fill_distiller(palette = "Spectral",direction = -1)+
  theme(axis.line = element_line(colour = "black"),axis.title = element_text(size=18),axis.text = element_text(size = 15),
        panel.border = element_blank(),panel.background = element_blank()) 
dev.off()

