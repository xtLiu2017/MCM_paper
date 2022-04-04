setwd("/Users/liux/Documents/MCM_paper/data")
library(Biostrings)
library(tidyverse)
library(biclust)
patterns = c('AAACCAGGGGC','AAACCAGGCGC','AAACCAGCGGC','AAACCACGGGC',
             'AAACCACCGGC','AAACCACGCGC','AAACCAGCCGC','AAACCACCCGC')
rev_patterns = c('GCCCCTGGTTT','GCGCCTGGTTT','GCCGCTGGTTT', 'GCCCGTGGTTT',
                 'GCCGGTGGTTT','GCGCGTGGTTT','GCGGCTGGTTT','GCGGGTGGTTT')
fasta = readDNAStringSet("sp_genes.fasta")
df.fasta = data.frame(fasta)
expla_vec = rep(0,66)
for (i in 1:length(fasta)){
  sequence = fasta[[i]]
  for (mypattern in patterns){
    if (grepl(mypattern, sequence)){
      expla_vec[i] = 1
      break
    }
  }
  for (mypattern in rev_patterns){
    if (grepl(mypattern, sequence)){
      expla_vec[i] = 1
      break
    }
  }
}
load("genes_to_fix.Rdata")
load("peak_info_MCM_spline.Rdata")
load("select.highquality.genes.Rdata")
load("peak_info_fixed.Rdata")
colnames(peak_info) = c('pt1_i1','pl1_i1','pt2_i1','pl2_i1')
peak_info_1939 = peak_info[select.highquality.genes,]
peak_info_1939[genes_to_be_fixed,] = peak_info_fixed[genes_to_be_fixed,]
peak_level_ratio = peak_info_1939[,'pl1_i1'] / peak_info_1939[,'pl2_i1']
library(clusterProfiler)
library(org.At.tair.db)
sp_genes = names(which(peak_level_ratio > 4))
#for single peak genes
mean(peak_info_1939[sp_genes, 'pt1_i1']) + 3
sd(peak_info_1939[sp_genes, 'pt1_i1'])
#for all genes in the background
mean(peak_info_1939[, 'pt1_i1'],na.rm = T) + 3
sd(peak_info_1939[, 'pt1_i1'],na.rm = T) 

sp.BP <- enrichGO(gene = sp_genes,universe = select.highquality.genes,
                  OrgDb = "org.At.tair.db", ont  = "BP", keyType = "TAIR",
                  pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = F)
x = org.At.tairSYMBOL
gene_symbol = c()
for (mygene in sp_genes){
  if (length(as.list(x[mygene])) == 0){
    gene_symbol = c(gene_symbol, NA)
  }else{
    gene_symbol = c(gene_symbol, as.list(x[mygene])[[1]][1])
  }
}
sp.BP$Description
sp.BP$geneID

heat_vec = c()
oxygen_vec = c()
for (mygene in sp_genes){
  heat = grepl(mygene, sp.BP$geneID)[c(7,12)]
  oxygen = grepl(mygene, sp.BP$geneID)[1:6]
  heat_vec = c(heat_vec, sum(heat))
  oxygen_vec = c(oxygen_vec, sum(oxygen))
}

heat_vec = binarize(heat_vec, threshold = 0.4)
oxygen_vec = binarize(oxygen_vec, threshold = 0.4)

load("Annotation_mat.Rdata")
table1 = data.frame(gene_names = sp_genes,
                    gene_symbol = gene_symbol,
                   heat_GO_terms = heat_vec, hypoxia_GO_terms = oxygen_vec, 
                   HSF21_binding = Annotation_mat['HSF21_col_a',sp_genes],
                   HSF3_binding = Annotation_mat['HSF3_col_a',sp_genes],
                   HSF6_binding = Annotation_mat['HSF6_col_a',sp_genes],
                   HSF7_binding = Annotation_mat['HSF7_col_a',sp_genes],
                   HSF6B_binding = Annotation_mat['HSFA6B_col',sp_genes],
                   HSFC1_binding = Annotation_mat['HSFC1_col_a',sp_genes],
                   ERF4_binding = Annotation_mat['ERF4_col_a',sp_genes],
                   AT1G28160_binding = Annotation_mat['AT1G28160_col_a', sp_genes])

write.table(table1, file = 'table1.txt',quote = F,row.names=FALSE)

