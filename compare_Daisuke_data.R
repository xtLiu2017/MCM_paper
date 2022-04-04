setwd("/Users/liux/Documents/MCM_paper/data")
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

peak_info_1939 = peak_info[select.highquality.genes,]
peak_info_1939[genes_to_be_fixed,] = peak_info_fixed[genes_to_be_fixed,]
peak_level_ratio = peak_info_1939[,'pl1_i1'] / peak_info_1939[,'pl2_i1']
sort_peak_level_ratio = sort(peak_level_ratio, na.last = NA)

library(clusterProfiler)
library(org.At.tair.db)
keytypes(org.At.tair.db)
data(geneList,package = "DOSE")
sp_genes = names(which(peak_level_ratio > 4))
mean(peak_info_1939[sp_genes, 'pt1_i1']) + 3
sd(peak_info_1939[sp_genes, 'pt1_i1'])

mean(peak_info_1939[, 'pt1_i1'],na.rm = T) + 3
sd(peak_info_1939[, 'pt1_i1'],na.rm = T) 

dp_genes = names(which(peak_level_ratio < 0.5))
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
  theme(panel.border = element_blank(), axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),plot.title = element_text(size=7,hjust=0.5),
        legend.text = element_text(size=6), legend.key.size = unit(0.3, 'cm'),legend.position='none') + 
  scale_x_discrete(limits=rev(df_GO_HSF$terms[1:12])) 

sp.CC <- enrichGO(gene = sp_genes,universe = select.highquality.genes,
                  OrgDb = "org.At.tair.db", ont  = "CC", keyType = "TAIR",
                  pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
sp.MF <- enrichGO(gene = sp_genes,universe = select.highquality.genes,
                  OrgDb = "org.At.tair.db", ont  = "MF", keyType = "TAIR",
                  pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
dp.BP <- enrichGO(gene = dp_genes,universe = select.highquality.genes,
                  OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                  pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
dp.CC <- enrichGO(gene = dp_genes,universe = select.highquality.genes,
                  OrgDb = "org.At.tair.db", ont  = "CC", keyType = "TAIR",
                  pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
dp.MF <- enrichGO(gene = dp_genes,universe = select.highquality.genes,
                  OrgDb = "org.At.tair.db", ont  = "MF", keyType = "TAIR",
                  pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
write(sp_genes,file = 'sp_genes.txt')
write(dp_genes,file = 'dp_genes.txt')
library("writexl")
write_xlsx(data.frame(sp.BP),"sp.BP.GP.xlsx")
load("Annotation_mat.Rdata")
all.tfs = rownames(Annotation_mat)
all.p = c()
all.direction = c()
for (mytf in all.tfs){
  set1 = Annotation_mat[mytf,sp_genes]
  set2 = Annotation_mat[mytf,setdiff(select.highquality.genes,sp_genes)]
  all.p = c(all.p, wilcox.test(set1,set2)$p.value)
  all.direction = c(all.direction, mean(set1) > mean(set2))
}
names(all.p) = all.tfs
library(qvalue)
myqvalues = p.adjust(all.p, method = 'BH')
indices = which(myqvalues < 0.05)
myqvalues[indices]
all.direction[indices]
#check out binding sites of HSF in hypoxia genes
hypoxia_genes = strsplit(sp.BP$geneID[1],'/')[[1]]
heat_response_genes = strsplit(sp.BP$geneID[7],'/')[[1]]
#"HSF21_col_a"           "HSF3_col_a"           
# "HSF6_col_a"            "HSF7_col_a"            "HSFA6B_col"            "HSFC1_col_a"    

Annotation_mat['HSF3_col_a',hypoxia_genes]
Annotation_mat['HSF6_col_a',hypoxia_genes]
Annotation_mat['HSFA6B_col',hypoxia_genes]
Annotation_mat['HSFC1_col_a',hypoxia_genes]

Annotation_mat['HSF3_col_a',heat_response_genes]
Annotation_mat['HSF6_col_a',heat_response_genes]
Annotation_mat['HSFA6B_col',heat_response_genes]
Annotation_mat['HSFC1_col_a',heat_response_genes]


Annotation_mat['HSF3_col_a',setdiff(hypoxia_genes,heat_response_genes)]
#There are 5 overlapping genes between heat response and hypoxia response, and they are bound 
#by HSF TFs. But other genes in hypoxia response do not have HSF binding, such as Camodulin 
# binding genes. So no evidence showing HSF binding derives the enrichment of hypoxia response genes.


#ERF-VII binding 
Annotation_mat['RAP212_col_a',hypoxia_genes]
Annotation_mat['HSF6_col_a',hypoxia_genes]
