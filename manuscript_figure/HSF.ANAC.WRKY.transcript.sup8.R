# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate supplemental figure 8
# The figure shows the transcript response of three families, HSF, WRKY, ANAC

# directory of MCM project, /MCM_paper/data
setwd("/Users/liux/Documents/MCM_paper/data")
# glm fitting, AvrRpt2 ETIs
load("glm_fixef_Ken_WT_with_pseudocounts.Rdata")

HSF_genes = c('AT3G22830', 'AT5G16820', 'AT3G24520', 'AT4G11660', 'AT5G62020')
names(HSF_genes) = c('HSFA6B', 'HSF3', 'HSFC1', 'HSF7', 'HSF6')
ANAC_genes = c('AT3G18400', 'AT5G13180', 'AT3G04070', 'AT3G15500', 'AT3G49530', 'AT2G24430')
names(ANAC_genes) = c('ANAC058', 'ANAC083', 'ANAC047', 'ANAC055', 'ANAC062', 'ANAC038')
WRKY_genes = c('AT4G24240','AT5G13080','AT2G23320','AT3G56400','AT5G41570')
names(WRKY_genes) = c('WRKY7', 'WRKY75', 'WRKY15', 'WRKY70', 'WRKY24')
# glm labels
AvrRpt2_labels = c('mytreatAvrRpt2:mytime03h', 'mytreatAvrRpt2:mytime04h', 'mytreatAvrRpt2:mytime06h',
                  'mytreatAvrRpt2:mytime09h', 'mytreatAvrRpt2:mytime12h', 'mytreatAvrRpt2:mytime16h',
                  'mytreatAvrRpt2:mytime20h', 'mytreatAvrRpt2:mytime24h')
mock_labels = c('mytreatmock:mytime03h', 'mytreatmock:mytime04h', 'mytreatmock:mytime06h',
                   'mytreatmock:mytime09h', 'mytreatmock:mytime12h', 'mytreatmock:mytime16h',
                   'mytreatmock:mytime20h', 'mytreatmock:mytime24h')
TF_labels = c(names(HSF_genes),names(ANAC_genes),names(WRKY_genes))
names(TF_labels) = c(HSF_genes,ANAC_genes,WRKY_genes)

# prepare data for ggplot
# plot mean, standard deviation
all.mean_diff = c()
all.std_diff = c()
TFs.glm.tair = c()
TFs.glm = c()
TF.family = c()
for (mytf in c(HSF_genes, ANAC_genes, WRKY_genes)){
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

# packages 
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
df_TFs = data.frame(time = rep(c(3,4,6,9,12,16,20,24), times = length(TF.family)/8) ,
                    mean_diff = all.mean_diff, std_diff = all.std_diff,
                    TF_family = TF.family, TF_label = TFs.glm)
col_HSF = brewer.pal(9,'Reds')[4:7]
col_ANAC = brewer.pal(9,'Blues')[4:7]
col_WRKY = brewer.pal(9,'Oranges')[3:6]
df_TFs$TF_family = factor(df_TFs$TF_family, levels = c("HSF",'ANAC','WRKY'))

plot_sup8 = ggplot(df_TFs, aes(x = time, y = mean_diff, color = TF_label)) + geom_line() + 
  geom_errorbar(aes(ymin=mean_diff- 0.5 * std_diff, ymax=mean_diff + 0.5 * std_diff)) + theme_bw() +
  facet_wrap( ~ TF_family, scales = 'free_x') + ylab(expression(log[2]~(AvrRpt2/mock))) + 
  xlab("time (hour)")  + 
  scale_color_manual(breaks = matrix(TFs.glm,ncol = 8,byrow = T)[,1], values=c(col_HSF, col_ANAC,col_WRKY)) +
  theme(legend.text = element_text(size=10),legend.title = element_blank(),
        strip.background = element_blank(),strip.text.x = element_text(size = 14),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10))
jpeg('/Users/liux/Documents/MCM_paper/manuscript_figure/sup.fig8.jpeg', width = 200, height = 100, unit = 'mm', res = 400)
grid.arrange(plot_sup8)
dev.off()
  





