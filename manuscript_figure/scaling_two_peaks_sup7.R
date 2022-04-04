# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate supplemental figure 7
# It shows a single scaling factor can account for the global peak delay

# directory of MCM project, /MCM_paper/data
setwd("/Users/liux/Documents/MCM_paper/data")

# packages
library(ComplexHeatmap)    # heatmap library
library(circlize)          # colors

# MCM fitting, long synthetic profiles, sampled at very dense time points
# optimize scaling parameter
load('complete_synthetic_profiles.Rdata')
# 1939 genes
load("select.highquality.genes.Rdata")

# first panel, sup7A begins
# sample at evenly distributed time points on linear scale
time_index = (seq(3,24,1) - 3)/0.01 + 1
selected_indices = time_index
# first and second peak profiles
mymat_1 = profiles_mat[select.highquality.genes, c(time_index, time_index + 2101)]
# color legend
col = colorRamp2(c(-2, 0, 1, 2, 6), c("cornflowerblue","white","orange",'orangered',"red4"))
dist_mat = dist(mymat_1)
mhclust = hclust(dist_mat,method = "average")
# show some but not all labels for time
mylabel_1 = c(seq(3,24,1),seq(3,24,1))
mylabel_1[seq(2,length(mylabel_1),4)] = ''
mylabel_1[seq(3,length(mylabel_1),4)] = ''
mylabel_1[seq(4,length(mylabel_1),4)] = ''
colnames(mymat_1) = mylabel_1
ha_1 = columnAnnotation(foo = anno_text(mylabel_1, gp = gpar(fontsize = 13), rot = 0, 
                                      location = unit(0.5, 'npc'), just = 'center'))
# heatmap
ht1 = Heatmap(mymat_1, col = col, cluster_rows = mhclust, name = "log(FC)",
              cluster_columns = FALSE, show_column_names = F,
              show_row_names = F,show_column_dend = FALSE, show_row_dend = T, 
              na_col = "#E0E0E0", row_title = paste(length(select.highquality.genes), "genes"),
              row_title_rot = 0, row_title_gp = gpar(fontsize=18),
              bottom_annotation = ha_1,
              row_dend_width = unit(15, "mm"), use_raster = F,height = unit(110,"mm"),
              width = unit(110,"mm"),column_title_gp = gpar(fontsize=18),
              border = "gray20",column_split = rep(c("first peaks","second peaks"),each=(dim(mymat_1)[2]/2)),
              show_heatmap_legend = F, 
              # show the peak time of the first peaks and second peaks for all genes
              # a dot in the heatmap showing the peak time
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(mymat_1[i, j] ==  max(mymat_1[i,1:(dim(mymat_1)[2]/2)]) & mymat_1[i, j] > 1){
                  if (j < 4){
                    grid.text('.', x, y, gp = gpar(fontsize = 6, col = "#0070C0"))
                  }else{
                    grid.text('.', x, y, gp = gpar(fontsize = 6, col = "cadetblue1"))
                  }
                }
                else if(mymat_1[i, j] ==  max(mymat_1[i,(dim(mymat_1)[2]/2 + 1):dim(mymat_1)[2]]) & max(mymat_1[i,1:(dim(mymat_1)[2]/2)]) > 1 & max(mymat_1[i,(dim(mymat_1)[2]/2 + 1):dim(mymat_1)[2]]) > 0.5)
                  if (j <= 33 ){
                    grid.text('.', x, y, gp = gpar(fontsize = 6, col = "#0070C0"))
                  }else{
                    grid.text('.', x, y, gp = gpar(fontsize = 6, col = "cadetblue1"))
                  }
              })
# color legend
lgd = Legend( col_fun = col, title = "log(FC)",legend_height = unit(35,"mm"),legend_width = unit(10,"mm"),
              title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14),title_gap = unit(3,"mm"))

# sup7A ends


# second panel, sup7B begins
# sample at evenly distributed time points on log scale
log_interval = (log(24) - log(3))/21
log_even_time_on_linear = round(exp(log_interval * seq(0,21) + log(3)), digits = 2)
#need to remove 3 hour above
time_index_log = (log_even_time_on_linear - 3) / 0.01 + 1
selected_indices_log = time_index_log
mymat_log = profiles_mat[select.highquality.genes, c(time_index_log, time_index_log + 2101)]

#note that we use the same clustering dengrogram, clustered based on the first heatmap panel, sup7A
dist_mat = dist(mymat_1)
mhclust = hclust(dist_mat,method = "average")
mylabel_log = round(rep(log_even_time_on_linear,2), digits = 1)
mylabel_log[seq(2,length(mylabel_log),4)] = ''
mylabel_log[seq(3,length(mylabel_log),4)] = ''
mylabel_log[seq(4,length(mylabel_log),4)] = ''
ha_log = columnAnnotation(foo = anno_text(mylabel_log, gp = gpar(fontsize = 13), rot = 0, 
                                        location = unit(0.5, 'npc'), just = 'center'))
ht_log = Heatmap(mymat_log, col = col, cluster_rows = mhclust, name = "log(FC)",
              cluster_columns = FALSE, show_column_names = F,
              show_row_names = F,show_column_dend = FALSE, show_row_dend = T, 
              na_col = "#E0E0E0", row_title = paste(length(select.highquality.genes), "genes"),
              row_title_rot = 0,bottom_annotation = ha_log,row_title_gp = gpar(fontsize=18),
              row_dend_width = unit(15, "mm"), use_raster = F,height = unit(110,"mm"),
              width = unit(110,"mm"),column_title_gp = gpar(fontsize=18),
              border = "gray20",column_split = rep(c("first peaks","second peaks"),each=dim(mymat_log)[2]/2),
              show_heatmap_legend = F, cell_fun = function(j, i, x, y, width, height, fill) {
                if(mymat_log[i, j] ==  max(mymat_log[i,1:(dim(mymat_log)[2]/2)]) & mymat_log[i, j] > 1){
                  if (j < 7){
                    grid.text('.', x, y, gp = gpar(fontsize = 6, col = "#0070C0"))
                  }else{
                    grid.text('.', x, y, gp = gpar(fontsize = 6, col = "cadetblue1"))
                  }
                }
                else if(mymat_log[i, j] ==  max(mymat_log[i,(dim(mymat_log)[2]/2+1):dim(mymat_log)[2]]) & max(mymat_log[i,1:(dim(mymat_log)[2]/2)]) > 1 & max(mymat_log[i,(dim(mymat_log)[2]/2+1):dim(mymat_log)[2]]) > 0.5){
                  if (j <= 35 ){
                    grid.text('.', x, y, gp = gpar(fontsize = 6, col = "#0070C0"))
                  }else{
                    grid.text('.', x, y, gp = gpar(fontsize = 6, col = "cadetblue1"))
                  }
                }
              })
# sup7B ends 

# the third heatmap panel, sup7C, begins
# cut some time points for the second peak, because it's so sparse
# prepare labels
mylabel_log_cut = round(rep(log_even_time_on_linear,2), digits = 1)[c(-23,-24,-25,-26,-27)]
mylabel_log_cut[seq(2,length(mylabel_log_cut),4)] = ''
mylabel_log_cut[seq(3,length(mylabel_log_cut),4)] = ''
mylabel_log_cut[seq(4,length(mylabel_log_cut),4)] = ''
ha_log_cut = columnAnnotation(foo = anno_text(mylabel_log_cut, gp = gpar(fontsize = 13), rot = 0, 
                                          location = unit(0.5, 'npc'), just = 'center'))
mymat_log_cut = mymat_log[,c(-23,-24,-25,-26,-27)]
# same clustering as previously used
ht_log_cut = Heatmap(mymat_log_cut, col = col, cluster_rows = mhclust, name = "log(FC)",
                 cluster_columns = FALSE, show_column_names = F,
                 show_row_names = F,show_column_dend = FALSE, show_row_dend = T, 
                 na_col = "#E0E0E0", row_title = paste(length(select.highquality.genes), "genes"),
                 row_title_rot = 0,bottom_annotation = ha_log_cut,row_title_gp = gpar(fontsize=18),
                 row_dend_width = unit(15, "mm"), use_raster = F,height = unit(110,"mm"),
                 width = unit(110,"mm"),column_title_gp = gpar(fontsize=18),
                 border = "gray20",column_split = rep(c("first peaks","second peaks"),times=c(22,17)),
                 show_heatmap_legend = F, cell_fun = function(j, i, x, y, width, height, fill) {
                   if(mymat_log_cut[i, j] ==  max(mymat_log_cut[i,1:22]) & mymat_log_cut[i, j] > 1){
                     if (j < 7){
                       grid.text('.', x, y, gp = gpar(fontsize = 6, col = "#0070C0"))
                     }else{
                       grid.text('.', x, y, gp = gpar(fontsize = 6, col = "cadetblue1"))
                     }
                   }
                   else if(mymat_log_cut[i, j] ==  max(mymat_log_cut[i,23:39]) & max(mymat_log_cut[i,1:22]) > 1 & max(mymat_log_cut[i,23:39]) > 0.5){
                     if (j <= 35 ){
                       grid.text('.', x, y, gp = gpar(fontsize = 6, col = "#0070C0"))
                     }else{
                       grid.text('.', x, y, gp = gpar(fontsize = 6, col = "cadetblue1"))
                     }
                 }})
#panel sup7C ends here

 
# Now we try to match the late peak profiles with early peak profiles
time_log = c(log_interval * seq(0,21) + log(3))
time_log_late = time_log[-(1:5)]
time_index_log_late = time_index_log[-(1:5)]
late_mat = profiles_mat[select.highquality.genes, time_index_log_late + 2101]
# save PCCs between early and late profiles for all the genes
all.true_pccs_time_lag = c()
# save PCCs between early and late profiles for all the genes
#all.true_pccs_time_lag.valid = c()
#try different time scaling factors
for (mytime_lag in seq(-2,0.3,0.02)){
    # based on the time lag, compute the time samples for the early peak
    early_time = round(exp(time_log_late  + mytime_lag), digits = 2)
    # determine the valid early times -- as some scaling factor will give early times out of range 3-24 hours
    valid_early_time = early_time[early_time >= 3 & early_time <= 24]
    valid_early_indices = (valid_early_time - 3) / 0.01 + 1
    # if the valid samples are so small, do not compute the PCC
    if (length(valid_early_indices) < 3){
      all.true_pccs_time_lag = c(all.true_pccs_time_lag, NA)
      next
    }
    # when the early samples are determined, need to pick the right number of late samples
    valid_late_indices = which(early_time >= 3 & early_time <= 24)
    valid_late_mat = late_mat[,valid_late_indices]
    valid_early_mat = profiles_mat[select.highquality.genes,valid_early_indices]
    # compute PCC for all genes
    true_pccs = c()
    for (i in 1:nrow(valid_early_mat)){
      x = valid_early_mat[i,]
      y = valid_late_mat[i,]
      true_pccs = c(true_pccs, cor(x,y,use="pair"))
    }
    # if many genes do not have valid PCC, should ignore this scaling factor value
    if (sum(is.na(true_pccs)) < 100){
        all.true_pccs_time_lag = c(all.true_pccs_time_lag,mean(true_pccs,na.rm = T))
    }else{
      all.true_pccs_time_lag = c(all.true_pccs_time_lag, NA)
    }
  }

# the scaling factor giving max PCC
best_ind = which.max(all.true_pccs_time_lag)
best_timelag = seq(-2,0.3,0.02)[best_ind]

#based on the optimal scaling factor, make the heatmap, the first peaks and second peaks
early_time = round(exp(time_log_late + best_timelag), digits = 2)
valid_early_time = early_time[early_time >= 3 & early_time <= 24]
valid_early_indices = (valid_early_time - 3) / 0.01 + 1
valid_late_indices = which(early_time >= 3 & early_time <= 24)
valid_late_mat = late_mat[,valid_late_indices]
valid_early_mat = profiles_mat[select.highquality.genes,valid_early_indices]
valid_late_time = exp(time_log_late[valid_late_indices])

mymat_best = cbind(valid_early_mat,valid_late_mat)
# same clustering as previously used
dist_mat = dist(mymat_1)
mhclust = hclust(dist_mat,method = "average")
mylabel_best = round(c(valid_early_time,valid_late_time), digits = 1)
mylabel_best[seq(2,length(mylabel_best),3)] = ''
mylabel_best[seq(3,length(mylabel_best),3)] = ''
ha_best = columnAnnotation(foo = anno_text(mylabel_best, gp = gpar(fontsize = 13), rot = 0, 
                                              location = unit(0.5, 'npc'), just = 'center'),
                           annotation_label = "mytime")
ht_best = Heatmap(mymat_best, col = col, cluster_rows = mhclust, name = "log(FC)",
                 cluster_columns = FALSE, show_column_names = F,
                 show_row_names = F,show_column_dend = FALSE, show_row_dend = T, 
                 na_col = "#E0E0E0", row_title = paste(length(select.highquality.genes), "genes"),
                 row_title_rot = 0, bottom_annotation = ha_best,row_title_gp = gpar(fontsize=18),
                 row_dend_width = unit(15, "mm"), use_raster = F,height = unit(110,"mm"),
                 width = unit(110,"mm"),column_title_gp = gpar(fontsize=18),
                 border = "gray20",column_split = rep(c("first peaks","second peaks"),each=dim(mymat_best)[2]/2),
                 show_heatmap_legend = F)

jpeg("/Users/liux/Documents/MCM_paper/manuscript_figure/sup.fig7.jpeg",height = 280, width = 400, units = "mm",res = 500)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht1, annotation_legend_list = lgd, ht_gap = unit(1.8, "cm"), newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht_log, annotation_legend_list = lgd, ht_gap = unit(1.8, "cm"),newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht_log_cut,annotation_legend_list = lgd, ht_gap = unit(1.8, "cm"), newpage = FALSE)

upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(ht_best,annotation_legend_list = lgd, ht_gap = unit(1.8, "cm"), newpage = FALSE)
upViewport()

dev.off()
# plot histogram showing all PCCs with different scaling factors
jpeg("/Users/liux/Documents/MCM_paper/manuscript_figure/sup.7e.jpeg",width = 130, height = 85, unit = "mm",res = 400)
par(mar = c(4.8,4.2,2,2))
plot(seq(-2,0.3,0.02), all.true_pccs_time_lag,type = "l", xlab = 'time lag (log)', 
     ylab = "mean PCC", cex.lab = 1, frame.plot = F,lwd = 1.4)
abline(v = -1.04, col = "red")
text('time lag = -1.04',x = -1.5, y = 0.25, col = "red")
dev.off()
