library(qvalue)
setwd("/Users/liux/Documents/MCM_paper/data")
load("glm_fixef_Ken_WT_with_pseudocounts.Rdata")
Mock_label = c("mytreatmock:mytime01h","mytreatmock:mytime02h","mytreatmock:mytime03h",
               "mytreatmock:mytime04h","mytreatmock:mytime06h","mytreatmock:mytime09h",
               "mytreatmock:mytime12h","mytreatmock:mytime16h","mytreatmock:mytime20h",
               "mytreatmock:mytime24h")
AvrRpt2_label = c("mytreatAvrRpt2:mytime01h","mytreatAvrRpt2:mytime02h","mytreatAvrRpt2:mytime03h",
                  "mytreatAvrRpt2:mytime04h","mytreatAvrRpt2:mytime06h","mytreatAvrRpt2:mytime09h",
                  "mytreatAvrRpt2:mytime12h","mytreatAvrRpt2:mytime16h","mytreatAvrRpt2:mytime20h",
                  "mytreatAvrRpt2:mytime24h")
EV_label = c("mytreatEV:mytime01h","mytreatEV:mytime02h","mytreatEV:mytime03h","mytreatEV:mytime04h",
             "mytreatEV:mytime06h","mytreatEV:mytime09h","mytreatEV:mytime12h","mytreatEV:mytime16h",
             "mytreatEV:mytime20h","mytreatEV:mytime24h")
AvrRpt2andEV_mock_diff <- lapply(glm_fixef,function(x){
  AvrRpt2_means <- x[AvrRpt2_label,1];AvrRpt2_stds <- x[AvrRpt2_label,2]
  EV_means <- x[EV_label,1];EV_stds <- x[EV_label,2]
  mock_means <- x[Mock_label,1];mock_stds <- x[Mock_label,2]
  AvrRpt2_diff_mean <- AvrRpt2_means - mock_means
  EV_diff_mean <- EV_means - mock_means
  AvrRpt2_diff_std <- sqrt(AvrRpt2_stds ^ 2 + mock_stds ^ 2)
  EV_diff_std <- sqrt(EV_stds ^ 2 + mock_stds ^ 2)
  AvrRpt2_diff_p <- 2 * pnorm(abs(AvrRpt2_diff_mean)/AvrRpt2_diff_std, lower.tail = F)
  EV_diff_p <- 2 * pnorm(abs(EV_diff_mean)/EV_diff_std, lower.tail = F)
  
  diff <- c(AvrRpt2_diff_mean,EV_diff_mean)
  std <- c(AvrRpt2_diff_std,EV_diff_std)
  p_values <- c(AvrRpt2_diff_p,EV_diff_p)
  return(cbind(diff,std,p_values))
})
p_AvrRpt2andEV_mock_diff = unlist(lapply(AvrRpt2andEV_mock_diff, '[', ,3))
q_AvrRpt2andEV_mock_diff = qvalue(p_AvrRpt2andEV_mock_diff)$qvalues
p_cutoff = sort(p_AvrRpt2andEV_mock_diff)[sum(q_AvrRpt2andEV_mock_diff < 0.05)]
respond_mock <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[,3] < p_cutoff & abs(x[,1]) > 1 )})
respond_mock_genes <- names(which(unlist(respond_mock) != 0))
not_respond_genes = setdiff(names(glm_fixef),respond_mock_genes)
flat_genes_profile = lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[,3] > p_cutoff)})
flat_genes = names(which(unlist(flat_genes_profile) > 19))
save(flat_genes,file = "flat_genes.Rdata")

AvrRpt2_mock_positive <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[1:10,3] < p_cutoff & x[1:10,1] > 1 )})
AvrRpt2_mock_negative <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[1:10,3] < p_cutoff & x[1:10,1] < -1 )})
EV_mock_positive <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[11:20,3] < p_cutoff & x[11:20,1] > 1 )})
EV_mock_negative <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[11:20,3] < p_cutoff & x[11:20,1] < -1 )})

AvrRpt2_mock_positive_genes1 <- names(which(unlist(AvrRpt2_mock_positive) != 0))
AvrRpt2_mock_negative_genes1 <- names(which(unlist(AvrRpt2_mock_negative) != 0))
AvrRpt2_mock_positive_genes <- setdiff(AvrRpt2_mock_positive_genes1,AvrRpt2_mock_negative_genes1)
AvrRpt2_mock_negative_genes <- setdiff(AvrRpt2_mock_negative_genes1,AvrRpt2_mock_positive_genes1)

EV_mock_positive_genes1 <- names(which(unlist(EV_mock_positive) != 0))
EV_mock_negative_genes1 <- names(which(unlist(EV_mock_negative) != 0))
EV_mock_positive_genes <- setdiff(EV_mock_positive_genes1,EV_mock_negative_genes1)

AvrRpt2EV_mock_positive_first_three_hours <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[c(1,2,3,11,12,13),3] < p_cutoff & x[c(1,2,3,11,12,13),1] > 1 )})
AvrRpt2EV_mock_positive_first_three_hours_genes <- names(which(unlist(AvrRpt2EV_mock_positive_first_three_hours) != 0))
AvrRpt2EV_mock_negative_first_three_hours <- lapply(AvrRpt2andEV_mock_diff,function(x){sum(x[c(1,2,3,11,12,13),3] < p_cutoff & x[c(1,2,3,11,12,13),1] < -1 )})
AvrRpt2EV_mock_negative_first_three_hours_genes <- names(which(unlist(AvrRpt2EV_mock_negative_first_three_hours) != 0))

AvrRpt2_mock_positive_genes_respond_before_3h <- intersect(AvrRpt2_mock_positive_genes,AvrRpt2EV_mock_positive_first_three_hours_genes)
AvrRpt2_mock_negative_genes_respond_before_3h <- intersect(AvrRpt2_mock_negative_genes,AvrRpt2EV_mock_negative_first_three_hours_genes)

AvrRpt2_mock_positive_genes <- setdiff(AvrRpt2_mock_positive_genes,AvrRpt2EV_mock_positive_first_three_hours_genes)
AvrRpt2_mock_negative_genes <- setdiff(AvrRpt2_mock_negative_genes,AvrRpt2EV_mock_negative_first_three_hours_genes)
EV_mock_positive_genes <- setdiff(EV_mock_positive_genes, AvrRpt2EV_mock_positive_first_three_hours_genes)

All.genes = names(glm_fixef)
save(AvrRpt2_mock_positive_genes,AvrRpt2_mock_positive_genes_respond_before_3h,All.genes,
     file = "AvrRpt2_genes.Rdata")
save(AvrRpt2_mock_negative_genes,AvrRpt2_mock_negative_genes_respond_before_3h,All.genes,
     file = "AvrRpt2_negative_genes.Rdata")


EV_mock_positive_first_time <- lapply(AvrRpt2andEV_mock_diff,function(x){min(which((x[11:20,3] < p_cutoff) & (x[11:20,1] > 1) ))})
EV_mock_positive_max_time <- lapply(AvrRpt2andEV_mock_diff,function(x){
  sig_indices = which((x[11:20,3] < p_cutoff) & (x[11:20,1] > 1))
  sig_indices[which.max(x[(10 + sig_indices),1])]
  })


save(EV_mock_positive_genes, EV_mock_positive_first_time, EV_mock_positive_max_time, All.genes, file = 'EV_genes.Rdata')
