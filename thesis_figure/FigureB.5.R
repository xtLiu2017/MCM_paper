# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate figure B.5 for my PhD thesis
# It is the transcript responses of the HSF genes

# directory of MCM project, /MCM_paper/data
setwd("/Users/liux/Documents/MCM_paper/data")

#load the glm data
load("glm_fixef_Ken_WT_with_pseudocounts.Rdata")

#prepare HSF data
HSF_genes = c('AT3G22830', 'AT5G16820', 'AT3G24520', 'AT4G11660', 'AT5G62020')
names(HSF_genes) = c('HSFA6B', 'HSF3', 'HSFC1', 'HSF7', 'HSF6')
# below are the labels in glm fitting at different time points
AvrRpt2_labels = c('mytreatAvrRpt2:mytime03h', 'mytreatAvrRpt2:mytime04h', 'mytreatAvrRpt2:mytime06h',
                  'mytreatAvrRpt2:mytime09h', 'mytreatAvrRpt2:mytime12h', 'mytreatAvrRpt2:mytime16h',
                  'mytreatAvrRpt2:mytime20h', 'mytreatAvrRpt2:mytime24h')
mock_labels = c('mytreatmock:mytime03h', 'mytreatmock:mytime04h', 'mytreatmock:mytime06h',
                   'mytreatmock:mytime09h', 'mytreatmock:mytime12h', 'mytreatmock:mytime16h',
                   'mytreatmock:mytime20h', 'mytreatmock:mytime24h')
TF_labels = names(HSF_genes)
names(TF_labels) = HSF_genes
all.mean_diff = c()
all.std_diff = c()
TFs.glm.tair = c()
TFs.glm = c()
TF.family = c()
for (mytf in HSF_genes){
  if (mytf %in% names(glm_fixef)){
    mean_Rpt2 = glm_fixef[[mytf]][AvrRpt2_labels, 1]
    mean_mock = glm_fixef[[mytf]][mock_labels, 1]
    std_Rpt2 = glm_fixef[[mytf]][AvrRpt2_labels, 2]
    std_mock = glm_fixef[[mytf]][mock_labels, 2]
    all.mean_diff = c(all.mean_diff, mean_Rpt2 - mean_mock)
    all.std_diff = c(all.std_diff,sqrt(std_Rpt2 ^ 2 + std_mock ^ 2) )
    TFs.glm.tair = c(TFs.glm.tair, rep(mytf, 8))
    TFs.glm = c(TFs.glm, rep(TF_labels[mytf], 8))
    if (grepl('HSF',TF_labels[mytf])){
      TF.family = c(TF.family, rep('HSF', 8))
    }else if (grepl('ANAC',TF_labels[mytf] )){
      TF.family = c(TF.family, rep('ANAC', 8))
    }else{
      TF.family = c(TF.family, rep('WRKY', 8))
    }
  }
}

#load packages
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
#the data frame that saves the HSF TFs data
df_TFs = data.frame(time = rep(c(3,4,6,9,12,16,20,24), times = length(TF.family)/8) ,
                    mean_diff = all.mean_diff, std_diff = all.std_diff,
                    TF_family = TF.family, TF_label = TFs.glm)

#obtain some colors for HSF TFs
col_HSF = brewer.pal(9,'Reds')[4:7]


plot_B5 = ggplot(df_TFs, aes(x = time, y = mean_diff, color = TF_label)) + geom_line() + 
  geom_errorbar(aes(ymin=mean_diff- 0.5 * std_diff, ymax=mean_diff + 0.5 * std_diff)) + theme_bw() +
  ylab(expression(log[2]~(AvrRpt2/mock))) + 
  xlab("time (hour)")  + 
  scale_color_manual(breaks = matrix(TFs.glm,ncol = 8,byrow = T)[,1], values=col_HSF) +
  theme(legend.text = element_text(size=10),legend.title = element_blank(),
        strip.background = element_blank(),strip.text.x = element_text(size = 14),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10))
# change the output directory as needed
jpeg('/Users/liux/Desktop/paper_figures/thesis/figureB.5.jpeg', width = 120, height = 100, unit = 'mm', res = 500)
grid.arrange(plot_B5)
dev.off()
  





