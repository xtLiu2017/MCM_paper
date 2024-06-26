# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate figure 6 as well as two supplemetal figure panels (Supplental figure 9A and B) for the manuscript

# set up the directory where I can get the data
setwd("~/Documents/MCM_paper/data")

#load packages
library(readxl)
library(ComplexHeatmap)    # make heatmaps, same below
library(circlize)
library(RColorBrewer)      # for colors
library(lsa)               # cosine function
library(vegan)             # flip the trees for Figure 6C

# 3039 genes upregulated at at least one time point 
load("AvrRpt2_genes.Rdata")
# preparing data
# WRKY subgroup - WRKY gene matching
WRKY_phyloge = list('1' = c('WRKY1','WRKY2','WRKY3','WRKY4','WRKY10','WRKY19','WRKY20','WRKY25','WRKY26',
                            'WRKY32','WRKY33','WRKY34','WRKY44','WRKY45','WRKY58','WRKY73'),
                    '2a' = c('WRKY18', 'WRKY40', 'WRKY60'), 
                    '2b' = c('WRKY6',	'WRKY9','WRKY31','WRKY36','WRKY42','WRKY47','WRKY61','WRKY72'),
                    '2c' = c('WRKY8', 'WRKY12', 'WRKY13', 'WRKY23', 'WRKY24', 'WRKY28', 'WRKY43', 
                             'WRKY48', 'WRKY49', 'WRKY50', 'WRKY51', 'WRKY56', 'WRKY57', 'WRKY59', 
                             'WRKY68', 'WRKY71', 'WRKY75'),
                    '2d' = c('WRKY7', 'WRKY11', 'WRKY15', 'WRKY17','WRKY21', 'WRKY39', 'WRKY74'),
                    '2e' = c('WRKY14', 'WRKY16', 'WRKY22', 'WRKY27','WRKY29', 'WRKY35', 'WRKY65', 'WRKY69'),
                    '3' = c('WRKY30', 'WRKY38', 'WRKY41', 'WRKY46', 'WRKY52', 'WRKY53', 'WRKY54', 'WRKY55', 
                            'WRKY62',	'WRKY63','WRKY64', 'WRKY66', 'WRKY67', 'WRKY70'))
# same thing above, but in a different R format
WRKY_phyloge_df = data.frame(wrky_name = c(WRKY_phyloge[['1']],WRKY_phyloge[['2a']],WRKY_phyloge[['2b']],
                                           WRKY_phyloge[['2c']],WRKY_phyloge[['2d']],WRKY_phyloge[['2e']],
                                           WRKY_phyloge[['3']]),
                             family = c(rep('1', length(WRKY_phyloge[['1']])),rep('2a', length(WRKY_phyloge[['2a']])),
                                        rep('2b', length(WRKY_phyloge[['2b']])), rep('2c', length(WRKY_phyloge[['2c']])),
                                        rep('2d', length(WRKY_phyloge[['2d']])), rep('2e', length(WRKY_phyloge[['2e']])),
                                        rep('3', length(WRKY_phyloge[['3']]))))
rownames(WRKY_phyloge_df) = WRKY_phyloge_df$wrky_name
# five prioritized WRKY genes 
WRKY_genes = c('AT2G23320','AT5G13080','AT5G41570','AT4G24240','AT3G56400')
WRKY_family = c("WRKY14", "WRKY15", "WRKY18", "WRKY20","WRKY21","WRKY22",
                "WRKY24", "WRKY25", "WRKY27","WRKY28", "WRKY29","WRKY3" , 
                "WRKY33","WRKY40", "WRKY45", "WRKY50", "WRKY55", "WRKY65",
                "WRKY7","WRKY70", "WRKY71", "WRKY75", "WRKY8" )

# This is downloaded from TAIR, gene family, TAIR name - family name matching
ALL.info = data.frame(read_xlsx(path = "gene_families_sep_29_09_update.xlsx"))

# i was planning to do WRKY and ANAC, but later I only focus on WRKY
# get the WRKY tair names, and gene symbols, make sure they match each other
ALL.WRKY_genes = c()
ALL.WRKY_names = c()
for (j in 1:dim(ALL.info)[1]){
  gene_tair_name = ALL.info[j, 6]
  gene_family_name = ALL.info[j, 1]
  gene_symbol = ALL.info[j, 3]
  if (grepl('WRKY', gene_family_name, ignore.case = T)){
    if(gene_tair_name %in% c(AvrRpt2_mock_positive_genes,AvrRpt2_mock_positive_genes_respond_before_3h)){
      ALL.WRKY_genes = c(ALL.WRKY_genes, gene_tair_name)
      ALL.WRKY_names = c(ALL.WRKY_names, gene_symbol)
    }
  }
}

#WRKY
#use upper letters
ALL.WRKY_genes = toupper(ALL.WRKY_genes)
#remove 'AT' substring and use upper letters
ALL.WRKY_names = unlist(lapply(toupper(ALL.WRKY_names),function(x){substr(x,start = 3, stop = nchar(x))}))
names(ALL.WRKY_names) = ALL.WRKY_genes
names(ALL.WRKY_genes) = ALL.WRKY_names

#prepare data above 
#done


#Figure 6A, start
#get the glm estimate of AvrRpt2 
load("glm_fixef_Ken_WT_with_pseudocounts.Rdata")
#change the name of glm
glm_fixef_ETI = glm_fixef
#labels for glm
AvrRpt2_labels = c('mytreatAvrRpt2:mytime03h','mytreatAvrRpt2:mytime04h','mytreatAvrRpt2:mytime06h',
                   'mytreatAvrRpt2:mytime09h','mytreatAvrRpt2:mytime12h','mytreatAvrRpt2:mytime16h',
                   'mytreatAvrRpt2:mytime20h','mytreatAvrRpt2:mytime24h')
mock_labels = c('mytreatmock:mytime03h','mytreatmock:mytime04h','mytreatmock:mytime06h',
                   'mytreatmock:mytime09h','mytreatmock:mytime12h','mytreatmock:mytime16h',
                   'mytreatmock:mytime20h','mytreatmock:mytime24h')
# WRKY glm (ETI) mat
WRKY_mat = c()
for (mygene in ALL.WRKY_genes){
  WRKY_mat = rbind(WRKY_mat, glm_fixef_ETI[[mygene]][AvrRpt2_labels,1] - glm_fixef_ETI[[mygene]][mock_labels,1])
}
# There are 35 wrky genes in the 3039 genes (significantly upregulated at at least one time point) 

WRKY_names_add_group = c()
for (i in ALL.WRKY_names){
  WRKY_names_add_group = c(WRKY_names_add_group, paste(i,WRKY_phyloge_df[i,"family"],sep = "_"))
}
rownames(WRKY_mat) = WRKY_names_add_group
colnames(WRKY_mat) = c('3h','4h','6h','9h','12h','16h','20h','24h')

six_colors = brewer.pal(7,'Set1')
names(six_colors) = c('1','2a','2b','2c','2d','2e','3')
row_colors = six_colors[WRKY_phyloge_df[ALL.WRKY_names, 'family']]

col_glm_avrrpt2 = colorRamp2(c(-4,0,8), c('deepskyblue','white','red'))

dist_mat_cosine = as.dist(1 - cosine(t(WRKY_mat)))
mhclust_cosine = hclust(dist_mat_cosine,method = "average")
ht_cosine = Heatmap(WRKY_mat,name = "log2(FC)",col = col_glm_avrrpt2, show_heatmap_legend = F,
                       cluster_columns = F, cluster_rows = mhclust_cosine,
                       column_title = expression(italic(Pto) ~ 'AvrRpt2-induced ETI'),
                       row_names_gp = gpar(fontsize = 8, col = row_colors), 
                       column_names_rot = 0,
                       row_title = paste(dim(WRKY_mat)[1],"genes"),row_title_rot = 0,
                       row_title_gp = gpar(fontsize = 16),
                       column_names_centered = T,column_title_gp = gpar(fontsize = 16),
                       heatmap_width = unit(105,'mm'))
# The above is figure 6A
# Figure 6A, end

#Figure 6B, start
load("glm.fixef.RData")
glm_fixef_flg22 = glm_fixef
WT_labels = c('flg22_genotypeJEPS:flg22_time0','flg22_genotypeJEPS:flg22_time1',
              'flg22_genotypeJEPS:flg22_time2','flg22_genotypeJEPS:flg22_time3',
              'flg22_genotypeJEPS:flg22_time5','flg22_genotypeJEPS:flg22_time9',
              'flg22_genotypeJEPS:flg22_time18')
fls2_labels = c('flg22_genotypefls2:flg22_time0','flg22_genotypefls2:flg22_time1',
                'flg22_genotypefls2:flg22_time2','flg22_genotypefls2:flg22_time3',
                'flg22_genotypefls2:flg22_time5','flg22_genotypefls2:flg22_time9',
                'flg22_genotypefls2:flg22_time18')
flg22_mat = c()
for (mygene in ALL.WRKY_genes){
  if (mygene %in% names(glm_fixef_flg22)){
    WT_exp = glm_fixef_flg22[[mygene]][WT_labels,1]
    fls2_exp = glm_fixef_flg22[[mygene]][fls2_labels,1]
    flg22_mat = rbind(flg22_mat, (WT_exp - fls2_exp))
  }else{
    flg22_mat = rbind(flg22_mat, rep(NA, 7))
  }
}
colnames(flg22_mat) = c('0h','1h','2h','3h','5h','9h','18h')
ht_flg22 = Heatmap(flg22_mat,name = "log2(FC)",col = col_glm_avrrpt2, 
                    cluster_columns = F, cluster_rows = mhclust_cosine,
                    column_title = 'flg22 response', na_col = 'gray90',show_heatmap_legend = F,
                    row_names_gp = gpar(fontsize = 8, col = row_colors), column_names_rot = 0,
                    row_title = paste(dim(WRKY_mat)[1],"genes"),row_title_rot = 0,row_title_gp = gpar(fontsize = 12),
                    column_names_centered = T,column_title_gp = gpar(fontsize = 16),
                    heatmap_width = unit(65,'mm'))
#Figure 6B, end

#Figure 6C and Supplemental Figure 
#look at binding profiles of WRKY genes
# Below is the binding profiles of 349 TFs across 18442 genes (all genes before induction gene selection)
load('Annotation_mat_18442genes.Rdata')
WRKY_binding_profile = Annotation_mat[,ALL.WRKY_genes]

# remove the WRKY TFs with no binding site in any of the WRKY genes
TF_kept = names(which(rowSums(WRKY_binding_profile) != 0))
# remove the ANAC TFs with no binding site in any of the WRKY genes
TF_kept_ANAC = names(which(rowSums(ANAC_binding_profile) != 0))
# get the binding profiles of WRKY TFs
WRKY_binding_profile_kept = WRKY_binding_profile[TF_kept,]
# get the binding profiles of ANAC TFs
ANAC_binding_profile_kept = WRKY_binding_profile[TF_kept,]
TF_names = unlist(lapply(rownames(WRKY_binding_profile_kept),function(x){strsplit(x,split = "_")[[1]][1]}))
rownames(WRKY_binding_profile_kept) = TF_names
colnames(WRKY_binding_profile_kept) = ALL.WRKY_names
rownames(ANAC_binding_profile_kept) = TF_names
#get WRKY-WRKY profiles
WRKY_WRKY_binding_profile = t(WRKY_binding_profile_kept[275:297,])
#get ANAC-WRKY profiles
WRKY_ANAC_binding_profile = t(ANAC_binding_profile_kept[275:297,])

#WRKY indicator, if the TF is WRKY, 1, otherwise 0.
WRKY_TF_indicator = c()
for (i in TF_names){
  if (grepl('WRKY', i)){
    WRKY_TF_indicator = c(WRKY_TF_indicator, 1)
  }else{
    WRKY_TF_indicator = c(WRKY_TF_indicator, 0)
  }
}
names(WRKY_TF_indicator) = rownames(WRKY_binding_profile_kept)
col_WRKY_indicator = colorRamp2(c(0,1), c('white','red'))
# This is for supplemental figures showing all TF binding profiles across WRKY genes
WRKY_TF_anno = columnAnnotation(WRKY = WRKY_TF_indicator, col = list(WRKY = col_WRKY_indicator),
                             annotation_name_rot = 0, annotation_name_side = 'left',
                             annotation_name_gp = gpar(fontsize = 7, col = "red"), 
                             annotation_label = "WRKY",show_legend = F )

#color scale of the heatmap                      
col_wrky_bf = colorRamp2(c(0,1,2,3,4), brewer.pal(9,'Oranges')[c(1,3,5,7,9)])
#color of columns, showing WRKY name subgroup
column_colors = six_colors[WRKY_phyloge_df[colnames(WRKY_WRKY_binding_profile), 'family']]

WRKY_binding_profile_kept_t = t(WRKY_binding_profile_kept)

#304 TFs - 35 WRKY genes, supplemental figure
ht_WRKY_bf_cosine = Heatmap(WRKY_binding_profile_kept_t,col = col_wrky_bf,name = "BS No.",
                               cluster_rows = mhclust_cosine, show_row_dend = F,
                        
                               row_title = paste(dim(WRKY_binding_profile_kept)[1],"TFs"),row_title_rot = 0,
                               row_names_gp = gpar(fontsize = 8, col = row_colors,fontface = 'italic'),row_title_gp = gpar(fontsize = 10),
                               bottom_annotation = WRKY_TF_anno, show_row_names = T,show_column_names = F,
                                column_names_rot = 45,column_title = 'ALL TF',
                               heatmap_width = unit(110,'mm'))
#23 ANACs - 35 WRKYs, supplemental figure
ht_wrky_anac_bf = Heatmap(WRKY_ANAC_binding_profile, name = "BS number", col = col_wrky_bf,
                          row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
                          column_names_gp = gpar(fontsize = 10),
                          column_names_rot = 45)

#cluster the columns and sort by number of BS
# sum of WRKY binding sites across WRKY genes for all WRKY TFs 
binding_site_no = colSums(WRKY_WRKY_binding_profile)
dist_mat_column_cosine = as.dist(1 - cosine(WRKY_WRKY_binding_profile))
mhclust_column = hclust(dist_mat_column_cosine, method = "average")
# reorder the dendrogram based on the total binding site number
mhclust_column_reorder = reorder(mhclust_column, -binding_site_no)
# Below is for figure 6C
ht_WRKY_WRKY_bf_cosine = Heatmap(WRKY_WRKY_binding_profile,col = col_wrky_bf,
                                    name = "BS No.", cluster_columns = mhclust_column_reorder,
                                 show_column_dend = F,
                                    cluster_rows = mhclust_cosine, show_row_dend = F,show_heatmap_legend = F,
                                    row_names_gp = gpar(fontsize = 12, col = row_colors,fontface = 'italic'),
                                    column_names_gp = gpar(fontsize = 10, col = column_colors), column_title = 'WRKY TF',
                                    show_row_names = T,show_column_names = T,
                                 column_title_gp = gpar(fontsize = 16),
                                    column_names_rot = 45,heatmap_width = unit(130,'mm'))
#Figure 6C and supplemental Figure end

#prepare individual legends for (A, B) and C
lgd_flg22 = Legend( col_fun = col_glm_avrrpt2, title = expression(log[2] ~ (FC)),legend_height = unit(30,"mm"),grid_width = unit(4,"mm"),
                    title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 11),title_gap = unit(2,"mm"))
lgd_BS = Legend( col_fun = col_wrky_bf, title = "BS number",legend_height = unit(30,"mm"),grid_width = unit(4,"mm"),
                 title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 11),title_gap = unit(2,"mm"))
pd = packLegend(lgd_flg22,lgd_BS, gap = unit(5, "mm"), direction = "vertical")
#Figure 6 combined
jpeg('~/Documents/MCM_paper/manuscript_figure/Figure6.jpeg',width = 360, height = 180, unit = 'mm', res = 400)
draw(ht_cosine + ht_flg22 + ht_WRKY_WRKY_bf_cosine,ht_gap = unit(0.8, "cm"), 
     heatmap_legend_list =pd)
grid.text('A', x = unit(7,'mm'),y = unit(175,'mm'),gp = gpar(fontface = 'bold', fontsize = 12))
grid.text('B', x = unit(119,'mm'),y = unit(175,'mm'),gp = gpar(fontface = 'bold',fontsize = 12))
grid.text('C', x = unit(189,'mm'),y = unit(175,'mm'),gp = gpar(fontface = 'bold',fontsize = 12))
dev.off()
#end

#Figure B.6 for thesis, also can be used for supplemental figure of the manuscript
# Please change the output name and directory as needed
# for thesis
jpeg('~/Desktop/paper_figures/thesis/FigureB.6.jpeg',width = 150, height = 120, unit = 'mm', res = 400)
draw(ht_WRKY_bf_cosine)
dev.off()
#Figure B.7 for thesis, also can be used for supplemental figure of the manuscript
jpeg('~/Desktop/paper_figures/thesis/FigureB.7.jpeg',width = 150, height = 120, unit = 'mm', res = 400)
draw(ht_wrky_anac_bf)
dev.off()

# for manuscript
jpeg('~/Documents/MCM_paper/manuscript_figure/sup.fig9A.jpeg',width = 150, height = 120, unit = 'mm', res = 400)
draw(ht_WRKY_bf_cosine)
dev.off()

jpeg('~/Documents/MCM_paper/manuscript_figure/sup.fig9B.jpeg',width = 150, height = 120, unit = 'mm', res = 400)
draw(ht_wrky_anac_bf)
dev.off()
