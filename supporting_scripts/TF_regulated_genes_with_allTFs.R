setwd("~/Documents/MCM_paper/data")
load("Annotation_mat.Rdata")
load("select.highquality.genes.Rdata")
WRKY_TFs = c("WRKY14_col_a", "WRKY15_col_b", "WRKY18_col_a", "WRKY20_col_a", "WRKY21_col",
             "WRKY22_col", "WRKY24_col_a", "WRKY25_col_a", "WRKY27_col_a", "WRKY28_col_a", 
             "WRKY29_col_a", "WRKY3_col_a", "WRKY33_col_a","WRKY40_col", "WRKY45_col_a",
             "WRKY50_col_a", "WRKY55_col_a", "WRKY65_col_a", "WRKY7_col", "WRKY70_col",
             "WRKY71_col_a", "WRKY75_col_a","WRKY8_col")
WRKY_mat = Annotation_mat[WRKY_TFs,select.highquality.genes]
WRKY_genes = names(which(apply(WRKY_mat, 2, sum) > 0))

ANAC_TFs = c("ANAC004_col_a", "ANAC005_col_a", "ANAC013_col_a", "ANAC016_col_a",
             "ANAC017_col", "ANAC020_col_v31", "ANAC028_col_a", "ANAC034_col",          
             "ANAC038_col_a", "ANAC042_col_a", "ANAC045_col_a", "ANAC046_col",
             "ANAC047_col", "ANAC050_col_v3a", "ANAC053_col", "ANAC055_col",
             "ANAC057_col_a", "ANAC058_col_a", "ANAC062_col", "ANAC070_col",
             "ANAC071_col_a", "ANAC075_col_a", "ANAC079_col", "ANAC083_col_v3a",
             "ANAC087_col_v3a", "ANAC092_col_v3a", "ANAC094_col_a", "ANAC096_col",
             "ANAC103_col_a")
ANAC_mat = Annotation_mat[ANAC_TFs,select.highquality.genes]
ANAC_genes = names(which(apply(ANAC_mat, 2, sum) > 0))

HSF_TFs = c("HSF21_col_a", "HSF3_col_a", "HSF6_col_a",           
            "HSF7_col_a", "HSFA6B_col", "HSFC1_col_a")
HSF_mat = Annotation_mat[HSF_TFs,select.highquality.genes]
HSF_genes = names(which(apply(HSF_mat, 2, sum) > 0))

CAMTA_TFs = c("CAMTA1_col_a", "CAMTA5_col_a")
CAMTA_mat = Annotation_mat[CAMTA_TFs,select.highquality.genes]
CAMTA_genes = names(which(apply(CAMTA_mat, 2, sum) > 0))


load("peak_info_all.Rdata")
load("genes_to_fix.Rdata")
peak_info_update = peak_info
peak_info_update[genes_to_be_fixed,] = peak_info_fixed[genes_to_be_fixed,]
general_df = data.frame(pkt1 = peak_info_update[select.highquality.genes,1],
                        pkt2 = peak_info_update[select.highquality.genes,3],
                        pl1_pl2 = peak_info_update[select.highquality.genes,2] / (abs(peak_info_update[select.highquality.genes,2]) + abs(peak_info_update[select.highquality.genes,4])))

library(ggplot2)
library(gridExtra)
HSF_df = data.frame(pkt1 = peak_info_update[HSF_genes,1], pkt2 = peak_info_update[HSF_genes,3],
                    pl1_pl2 = peak_info_update[HSF_genes,2] / (abs(peak_info_update[HSF_genes,2]) + abs(peak_info_update[HSF_genes,4])))

p1 = ggplot(data = HSF_df, aes(x = pkt1 + 3)) + 
  geom_histogram(data = general_df, aes(x = pkt1 + 3, fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = "1"),bins = 100, size = 0.1) + theme_bw() + 
  xlab('peak 1 time (hour)') + ylab("Number of genes") + scale_fill_manual(values=c("red", "grey"),labels = c("HSF-regulated genes", "Background genes")) + 
  labs(fill = "Gene set") + theme(legend.position = 'none',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm")) + 
  ggtitle("HSF regulated genes") + 
  geom_vline(xintercept = median(HSF_df$pkt1, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pkt1,na.rm = T) + 3, linetype="dashed", color = "gray20") 

p2 = ggplot(data = HSF_df, aes(x = pkt2 + 3)) + 
  geom_histogram(data = general_df, aes(x = pkt2 + 3, fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = "1"),bins = 100, size = 0.1) + theme_bw() + 
  xlab('peak 2 time (hour)') + ylab("Number of genes")  + scale_fill_manual(values=c("red", "grey"),labels = c("HSF-regulated genes", "Background genes")) + 
  labs(fill = "Gene set") + theme(legend.position = 'none',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm")) + 
  geom_vline(xintercept = median(HSF_df$pkt2, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pkt2, na.rm = T) + 3, linetype="dashed", color = "gray20")

p3 = ggplot(data = HSF_df, aes(x = pl1_pl2 + 3)) + 
  geom_histogram(data = general_df, aes(x = pl1_pl2 + 3,fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = '1'),bins = 100,size = 0.1 ) + theme_bw() + 
  xlab('peak 1 level / (peak 1 level + peak 2 level)') + ylab("Number of genes") + 
  scale_fill_manual(values=c("red", "grey"),labels = c("HSF-regulated genes", "Background genes")) + 
  labs(fill = " ") + theme(legend.position = 'bottom',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),legend.key.size = unit(0.3, 'cm')) + 
  geom_vline(xintercept = median(HSF_df$pl1_pl2, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pl1_pl2,na.rm = T) + 3, linetype="dashed", color = "gray20")


ANAC_df = data.frame(pkt1 = peak_info_update[ANAC_genes,1], pkt2 = peak_info_update[ANAC_genes,3],
                    pl1_pl2 = peak_info_update[ANAC_genes,2] / (abs(peak_info_update[ANAC_genes,2]) + abs(peak_info_update[ANAC_genes,4])))

p4 = ggplot(data = ANAC_df, aes(x = pkt1 + 3)) + 
  geom_histogram(data = general_df, aes(x = pkt1 + 3, fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = "1"),bins = 100, size = 0.1) + theme_bw() + 
  xlab('peak 1 time (hour)') + ylab("Number of genes") + scale_fill_manual(values=c("red", "grey"),labels = c("ANAC-regulated genes", "Background genes")) + 
  labs(fill = "Gene set") + theme(legend.position = 'none',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm")) + 
  ggtitle("ANAC regulated genes") + 
  geom_vline(xintercept = median(ANAC_df$pkt1, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pkt1,na.rm = T) + 3, linetype="dashed", color = "gray20") 

p5 = ggplot(data = ANAC_df, aes(x = pkt2 + 3)) + 
  geom_histogram(data = general_df, aes(x = pkt2 + 3, fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = "1"),bins = 100, size = 0.1) + theme_bw() + 
  xlab('peak 2 time (hour)') + ylab("Number of genes")  + scale_fill_manual(values=c("red", "grey"),labels = c("ANAC-regulated genes", "Background genes")) + 
  labs(fill = "Gene set") + theme(legend.position = 'none',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm")) + 
  geom_vline(xintercept = median(ANAC_df$pkt2, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pkt2, na.rm = T) + 3, linetype="dashed", color = "gray20")

p6 = ggplot(data = ANAC_df, aes(x = pl1_pl2 + 3)) + 
  geom_histogram(data = general_df, aes(x = pl1_pl2 + 3,fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = '1'),bins = 100,size = 0.1 ) + theme_bw() + 
  xlab('peak 1 level / (peak 1 level + peak 2 level)') + ylab("Number of genes") + 
  scale_fill_manual(values=c("red", "grey"),labels = c("ANAC-regulated genes", "Background genes")) + 
  labs(fill = " ") + theme(legend.position = 'bottom',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),legend.key.size = unit(0.3, 'cm')) + 
  geom_vline(xintercept = median(ANAC_df$pl1_pl2, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pl1_pl2,na.rm = T) + 3, linetype="dashed", color = "gray20")


CAMTA_df = data.frame(pkt1 = peak_info_update[CAMTA_genes,1], pkt2 = peak_info_update[CAMTA_genes,3],
                     pl1_pl2 = peak_info_update[CAMTA_genes,2] / (abs(peak_info_update[CAMTA_genes,2]) + abs(peak_info_update[CAMTA_genes,4])))

p7 = ggplot(data = CAMTA_df, aes(x = pkt1 + 3)) + 
  geom_histogram(data = general_df, aes(x = pkt1 + 3, fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = "1"),bins = 100, size = 0.1) + theme_bw() + 
  xlab('peak 1 time (hour)') + ylab("Number of genes") + scale_fill_manual(values=c("red", "grey"),labels = c("CAMTA-regulated genes", "Background genes")) + 
  labs(fill = "Gene set") + theme(legend.position = 'none',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm")) + 
  ggtitle("CAMTA regulated genes") + 
  geom_vline(xintercept = median(CAMTA_df$pkt1, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pkt1,na.rm = T) + 3, linetype="dashed", color = "gray20") 

p8 = ggplot(data = CAMTA_df, aes(x = pkt2 + 3)) + 
  geom_histogram(data = general_df, aes(x = pkt2 + 3, fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = "1"),bins = 100, size = 0.1) + theme_bw() + 
  xlab('peak 2 time (hour)') + ylab("Number of genes")  + scale_fill_manual(values=c("red", "grey"),labels = c("CAMTA-regulated genes", "Background genes")) + 
  labs(fill = "Gene set") + theme(legend.position = 'none',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm")) + 
  geom_vline(xintercept = median(CAMTA_df$pkt2, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pkt2, na.rm = T) + 3, linetype="dashed", color = "gray20")

p9 = ggplot(data = CAMTA_df, aes(x = pl1_pl2 + 3)) + 
  geom_histogram(data = general_df, aes(x = pl1_pl2 + 3,fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = '1'),bins = 100,size = 0.1 ) + theme_bw() + 
  xlab('peak 1 level / (peak 1 level + peak 2 level)') + ylab("Number of genes") + 
  scale_fill_manual(values=c("red", "grey"),labels = c("CAMTA-regulated genes", "Background genes")) + 
  labs(fill = " ") + theme(legend.position = 'bottom',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),legend.key.size = unit(0.3, 'cm')) + 
  geom_vline(xintercept = median(CAMTA_df$pl1_pl2, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pl1_pl2,na.rm = T) + 3, linetype="dashed", color = "gray20")


WRKY_df = data.frame(pkt1 = peak_info_update[WRKY_genes,1], pkt2 = peak_info_update[WRKY_genes,3],
                      pl1_pl2 = peak_info_update[WRKY_genes,2] / (abs(peak_info_update[WRKY_genes,2]) + abs(peak_info_update[WRKY_genes,4])))

p10 = ggplot(data = WRKY_df, aes(x = pkt1 + 3)) + 
  geom_histogram(data = general_df, aes(x = pkt1 + 3, fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = "1"),bins = 100, size = 0.1) + theme_bw() + 
  xlab('peak 1 time (hour)') + ylab("Number of genes") + scale_fill_manual(values=c("red", "grey"),labels = c("WRKY-regulated genes", "Background genes")) + 
  labs(fill = "Gene set") + theme(legend.position = 'none',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm")) + 
  ggtitle("WRKY regulated genes") + 
  geom_vline(xintercept = median(WRKY_df$pkt1, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pkt1,na.rm = T) + 3, linetype="dashed", color = "gray20") 

p11 = ggplot(data = WRKY_df, aes(x = pkt2 + 3)) + 
  geom_histogram(data = general_df, aes(x = pkt2 + 3, fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = "1"),bins = 100, size = 0.1) + theme_bw() + 
  xlab('peak 2 time (hour)') + ylab("Number of genes")  + scale_fill_manual(values=c("red", "grey"),labels = c("WRKY-regulated genes", "Background genes")) + 
  labs(fill = "Gene set") + theme(legend.position = 'none',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm")) + 
  geom_vline(xintercept = median(WRKY_df$pkt2, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pkt2, na.rm = T) + 3, linetype="dashed", color = "gray20")

p12 = ggplot(data = WRKY_df, aes(x = pl1_pl2 + 3)) + 
  geom_histogram(data = general_df, aes(x = pl1_pl2 + 3,fill = '2'),alpha = .8, bins = 150) + 
  geom_histogram(aes(fill = '1'),bins = 100,size = 0.1 ) + theme_bw() + 
  xlab('peak 1 level / (peak 1 level + peak 2 level)') + ylab("Number of genes") + 
  scale_fill_manual(values=c("red", "grey"),labels = c("WRKY-regulated genes", "Background genes")) + 
  labs(fill = " ") + theme(legend.position = 'bottom',plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),legend.key.size = unit(0.3, 'cm')) + 
  geom_vline(xintercept = median(WRKY_df$pl1_pl2, na.rm = T) + 3, linetype="dashed", color = "red") + 
  geom_vline(xintercept = median(general_df$pl1_pl2,na.rm = T) + 3, linetype="dashed", color = "gray20")

lay_mat = rbind(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),
                c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),
                c(1,2,3,4),
                c(5,6,7,8),c(5,6,7,8),c(5,6,7,8),
                c(5,6,7,8),c(5,6,7,8),c(5,6,7,8),
                c(9,10,11,12),c(9,10,11,12),c(9,10,11,12),
                c(9,10,11,12),c(9,10,11,12),c(9,10,11,12),
                c(9,10,11,12),c(9,10,11,12))

jpeg("figures/TF_regulated_genes_allTFs.jpeg",width = 350, height = 150, unit = 'mm',res = 400 )
grid.arrange(p1,p4,p7,p10,p2,p5,p8,p11,p3,p6,p9,p12, nrow = 3,layout_matrix = lay_mat)
dev.off()








