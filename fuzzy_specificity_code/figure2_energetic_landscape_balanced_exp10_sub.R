################################################################################################
#Source code for Figure2 (+ related supplementary material) from 
#"A complete map of specificity encoding for a partially fuzzy protein interaction"
#Zarin, T and Lehner, B. 2024
################################################################################################
#please download required data (Zenodo link) + see instructions/required packages at 
#https://github.com/lehner-lab/fuzzy_specificity/
################################################################################################
#remember to set directory 
#e.g. setwd("dir_containing_fuzzy_specificity_code_dir_AND_fuzzy_specificity_data_dir")
################################################################################################
source("./fuzzy_specificity_code/source_functions/basic_functions_fuzzy_specificity_TZ.R")
#load in variant fitness tables of cript n and cript c
cript_n <- read.delim("./fuzzy_specificity_data/criptn_norm_for_mochi_added_knyk.tsv")
cript_c <- read.delim("./fuzzy_specificity_data/criptc_norm_for_mochi.tsv")
#bind them together, delineated by "block1, block2"
cript_n$data_id <- rep("cript_n", nrow(cript_n))
cript_c$data_id <- rep("cript_c", nrow(cript_c))
n_c_combined_untrimmed <- rbind(cript_n, cript_c)
n_c_combined <- rbind(cript_n, cript_c)
#get the sequences trimmed
n_c_combined[which(n_c_combined$data_id=="cript_n"),]$aa_seq <- unlist(lapply(n_c_combined[which(n_c_combined$data_id=="cript_n"),]$aa_seq, FUN=function(x){str_sub(x, 1,4)}))
n_c_combined[which(n_c_combined$data_id=="cript_c"),]$aa_seq <- unlist(lapply(n_c_combined[which(n_c_combined$data_id=="cript_c"),]$aa_seq, FUN=function(x){str_sub(x, 5,8)}))
#want to fix the AA hamming distance for cript n (since stand-in wildtype was "SGRR" high count variant)
cript_n_diffs <- lapply(n_c_combined[which(n_c_combined$data_id=="cript_n"),]$aa_seq, list.string.diff, "KNYK")
cript_n_aa <- unlist(lapply(cript_n_diffs, nrow))
n_c_combined[which(n_c_combined$data_id=="cript_n"),]$Nham_aa <- cript_n_aa
cript_n_subset <- n_c_combined[which(n_c_combined$data_id=="cript_n"),]
cript_c_subset <- n_c_combined[which(n_c_combined$data_id=="cript_c"),]
#take out stop codons
cript_n_subset <- cript_n_subset[-which(cript_n_subset$STOP==T),]
cript_c_subset <- cript_c_subset[-which(cript_c_subset$STOP==T),]
###

#mochi linear model
pred_file_nclinear_1 <- "./fuzzy_specificity_data/mochi_results/cript_nc_balanced_exp10linear_1/task_1/predictions/predicted_phenotypes_all.txt"
cript_nc_1_linear <- combine_read_delim_heldout_pred(pred_file_nclinear_1)
mochi_obspred_linear_order_1_c <- make_gg_scatter_with_xy_line(cript_nc_1_linear$fitness[which(cript_nc_1_linear$criptc_combo==1)], cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptc_combo==1)], "Observed fitness", "Predicted fitness")
mochi_obspred_linear_order_1_n <- make_gg_scatter_with_xy_line(cript_nc_1_linear$fitness[which(cript_nc_1_linear$criptn_combo==1)], cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptn_combo==1)], "Observed fitness", "Predicted fitness")
#want to plot the residuals to these fits to quantify bias
resids_linear_order_1_c <- make_gg_scatterplot_no_corr(cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptc_combo==1)], cript_nc_1_linear$fitness[which(cript_nc_1_linear$criptc_combo==1)] - cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptc_combo==1)], "Predicted fitness", "Residual")+
  geom_hline(yintercept=0,color="darkgrey", linewidth=0.9)
resids_linear_order_1_n <- make_gg_scatterplot_no_corr(cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptn_combo==1)], cript_nc_1_linear$fitness[which(cript_nc_1_linear$criptn_combo==1)] - cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptn_combo==1)], "Predicted fitness", "Residual")+
  geom_hline(yintercept=0,color="darkgrey", linewidth=0.9)

#mochi results files
pred_file_nc_1 <- "./fuzzy_specificity_data/mochi_results/cript_nc_balanced_exp10_1/task_1/predictions/predicted_phenotypes_all.txt"
dg_file_nc_1 <- "./fuzzy_specificity_data/mochi_results/cript_nc_balanced_exp10_1/task_1/weights/weights_Folding.txt"
pred_file_nc_12 <- "./fuzzy_specificity_data/mochi_results/cript_nc_balanced_exp10_2/task_1/predictions/predicted_phenotypes_all.txt"
dg_file_nc_12 <- "./fuzzy_specificity_data/mochi_results/cript_nc_balanced_exp10_2/task_1/weights/weights_Folding.txt"
#read in prediction data from mochi
cript_nc_1 <- combine_read_delim_heldout_pred(pred_file_nc_1)
cript_nc_12 <- combine_read_delim_heldout_pred(pred_file_nc_12)
#additive trait data table from mochi
dg_criptnc_1 <- read.delim(dg_file_nc_1)
dg_cript_nc_12 <- read.delim(dg_file_nc_12)
#linear model using the subset of data for balanced cript c
first_order_c_balanced <- get_linear_model(cript_c_subset[match(unlist(lapply(cript_nc_1[which(cript_nc_1$phenotype==2),]$aa_seq, str_sub, 5, 8)), cript_c_subset$aa_seq),], "first", "QTSV")
lm_plot_firstorder_c_balanced <- make_gg_scatterplot(cript_c_subset[match(unlist(lapply(cript_nc_1[which(cript_nc_1$phenotype==2),]$aa_seq, str_sub, 5, 8)), cript_c_subset$aa_seq),]$fitness, first_order_c_balanced$fitted.values, "Observed fitness", "Predicted fitness")
#plot predicted vs. observed allowing mochi first+second order interactions
#cript_c
mochi_obspred_order1_c <- make_gg_scatter_with_xy_line(cript_nc_1$fitness[which(cript_nc_1$criptc_combo==1)], cript_nc_1$predicted_phen[which(cript_nc_1$criptc_combo==1)], "Observed fitness", "Predicted fitness")
#plot the additive trait vs. fitness
mochi_obs_dddg_order1_c <- make_gg_scatterplot_no_corr(apply(cript_nc_1[which(cript_nc_1$criptc_combo==1),50:59], 1, mean), cript_nc_1$fitness[which(cript_nc_1$criptc_combo==1)], "Additive trait",  "Observed fitness")
#first + second order mochi
mochi_obspred_order12_c <- make_gg_scatter_with_xy_line(cript_nc_12$fitness[which(cript_nc_12$criptc_combo==1)], cript_nc_12$predicted_phen[which(cript_nc_12$criptc_combo==1)], "Observed fitness", "Predicted fitness")
#plot the additive trait vs. fitness
mochi_obs_dddg_order12_c <- make_gg_scatterplot_no_corr(apply(cript_nc_12[which(cript_nc_12$criptc_combo==1), 50:59], 1, mean), cript_nc_12$fitness[which(cript_nc_12$criptc_combo==1)], "Additive trait", "Observed fitness")
#cript_n
mochi_obspred_order1_n <- make_gg_scatter_with_xy_line(cript_nc_1$fitness[which(cript_nc_1$criptn_combo==1)], cript_nc_1$predicted_phen[which(cript_nc_1$criptn_combo==1)], "Observed fitness", "Predicted fitness")
#plot the additive trait vs. fitness
mochi_obs_dddg_order1_n <- make_gg_scatterplot_no_corr(apply(cript_nc_1[which(cript_nc_1$criptn_combo==1),50:59], 1, mean), cript_nc_1$fitness[which(cript_nc_1$criptn_combo==1)], "Additive trait", "Observed fitness")
#plot predicted vs. observed allowing mochi first+second order interactions
#cript_c
mochi_obspred_order12_n <- make_gg_scatter_with_xy_line(cript_nc_12$fitness[which(cript_nc_12$criptn_combo==1)], cript_nc_12$predicted_phen[which(cript_nc_12$criptn_combo==1)], "Observed fitness", "Predicted fitness")
#plot the additive trait vs. fitness
mochi_obs_dddg_order12_n <- make_gg_scatterplot_no_corr(apply(cript_nc_12[which(cript_nc_12$criptn_combo==1),50:59], 1, mean), cript_nc_12$fitness[which(cript_nc_12$criptn_combo==1)], "Additive trait", "Observed fitness")
#resids from mochi
resids_mochi_order_1_c <- make_gg_scatterplot_no_corr(cript_nc_1$predicted_phen[which(cript_nc_1$criptc_combo==1)], cript_nc_1$fitness[which(cript_nc_1$criptc_combo==1)] - cript_nc_1$predicted_phen[which(cript_nc_1$criptc_combo==1)], "Predicted fitness", "Residual")+
  geom_hline(yintercept=0,color="darkgrey", linewidth=0.9)
resids_mochi_order_1_n <- make_gg_scatterplot_no_corr(cript_nc_1$predicted_phen[which(cript_nc_1$criptn_combo==1)], cript_nc_1$fitness[which(cript_nc_1$criptn_combo==1)] - cript_nc_1$predicted_phen[which(cript_nc_1$criptn_combo==1)], "Predicted fitness", "Residual")+
  geom_hline(yintercept=0,color="darkgrey", linewidth=0.9)




#compare first order terms first and first+2nd order model
matched_dg_12 <- match_datasets_dgs(dg_criptnc_1,dg_cript_nc_12)[-1,]
matched_dg_balanced <- make_gg_scatter_with_xy_line(matched_dg_12$x_mean_kcal.mol, matched_dg_12$y_mean_kcal.mol, "first order model balanced C", "first+2nd order model balanced C")
matched_dg_balanced_errors <- make_scatterplot_errorbars(matched_dg_12, matched_dg_12$x_mean_kcal.mol, matched_dg_12$y_mean_kcal.mol, lbl = " ", error_y = matched_dg_12$y_95ci, error_x = matched_dg_12$x_95ci, xname = "1st order model ddGs", yname="2nd order model ddGs") + 
  ggtitle(bquote(R^2 ~ .(paste(round((cor(matched_dg_12$x_mean_kcal.mol, matched_dg_12$y_mean_kcal.mol)^2), digits=3), " | N =", length(matched_dg_12$x_mean_kcal.mol), sep=" "))))
#also want to show error increases with ddG and correlation is limited to AA identity
mean_vs_error_dg_df <- dg_criptnc_1[-1,]
mean_vs_error_dg_df$id_pos <- paste(unlist(lapply(mean_vs_error_dg_df$id, FUN=function(x){str_sub(x,1,1)})),
                                    unlist(lapply(mean_vs_error_dg_df$id, get_mut_pos))-8,sep="")
mean_vs_error_dg_df$id_pos_col <- cript_cols_matrix[match(mean_vs_error_dg_df$id_pos, names(cript_cols_matrix))]
mean_vs_error_dg_df$id_pos_col_mut <- unlist(lapply(mean_vs_error_dg_df$id, change_mut_numbers, -8))
mean_vs_error_dg_df$id_pos_col_mut_conditional <- mean_vs_error_dg_df$id_pos_col_mut
mean_vs_error_dg_df$id_pos_col_mut_conditional[which(mean_vs_error_dg_df$ci95_kcal.mol<1)] <- " "
library(ggrepel)
make_scatterplot_errorbars_cols_v2 <- function(df, xaxis, yaxis, lbl, error_y, xname, yname, col_col){
  scplot <- ggplot(df, aes(x=xaxis, y=yaxis, label=lbl, fill=col_col))+
    geom_point(shape = 21, size=2, colour="black", stroke=0.5)+
    scale_fill_manual(values=cript_cols_matrix, limits=names(cript_cols_matrix), name= " ")+
    geom_errorbar(aes(ymin=yaxis-error_y, ymax=yaxis+error_y))+
    xlab(as.character(xname))+
    ylab(as.character(yname))+
    theme(plot.title = element_text(size=12, family="serif")) + 
    theme_bw()+
    #geom_text(nudge_x = -0.5, size=2.5)+
    ggrepel::geom_text_repel(segment.color = 'transparent', size=2.5, nudge_x = -0.2)+
    theme(legend.position="right")
  scplot + geom_hline(yintercept=1, color="darkgrey", size=0.8, linetype="dashed")
}
cript_dg_errors_1 <- make_scatterplot_errorbars_cols_v2(mean_vs_error_dg_df, mean_vs_error_dg_df$mean_kcal.mol, mean_vs_error_dg_df$ci95_kcal.mol, lbl = mean_vs_error_dg_df$id_pos_col_mut_conditional, error_y = 0, xname="ddG", yname="Error (95% confidence interval)", col_col = mean_vs_error_dg_df$id_pos)

#read in the mochi energies
all_dgs <- dg_cript_nc_12
all_dgs <- all_dgs[-which(all_dgs$id=="WT"),]
cript_n_dgs <- all_dgs[which(str_sub(all_dgs$Pos, 1,1) %in% c(1:4)),]
cript_c_dgs <- all_dgs[which(str_sub(all_dgs$Pos, 1,1) %in% c(5:8)),]
#sm_conf_dgs_nc <- all_dgs[intersect(which(nchar(all_dgs$id)<6), which(all_dgs$ci95_kcal.mol<1)),]
sm_conf_dgs_nc <- all_dgs[which(nchar(all_dgs$id)<6),]
sm_conf_dgs_nc$id <- unlist(lapply(sm_conf_dgs_nc$id, change_mut_numbers, -8))
single_mut_df_nc <- data.frame(mean_kcal.mol=sm_conf_dgs_nc$mean_kcal.mol, ci95_kcal.mol=sm_conf_dgs_nc$ci95_kcal.mol)
single_mut_df_nc$pos <- str_sub(sm_conf_dgs_nc$id, 1,-2)
single_mut_df_nc$mut <- str_sub(sm_conf_dgs_nc$id, -1,-1)
single_mut_df_nc$mut <- factor(single_mut_df_nc$mut, levels=new_order_aa)
single_mut_df_nc$pos <- factor(single_mut_df_nc$pos, levels=c("K-7", "N-6", "Y-5", "K-4", "Q-3", "T-2", "S-1", "V0"))
single_mut_df_nc$is_conf <- rep("± >1", nrow(single_mut_df_nc))
single_mut_df_nc$is_conf[which(single_mut_df_nc$ci95_kcal.mol<1)] <- " "
pheat_sm_nc <- ggplot(single_mut_df_nc,aes(x=pos,y=mut, fill=mean_kcal.mol))+
  geom_tile()+
  scale_fill_gradient2(low="orangered", mid="white", high="midnightblue", na.value = "darkgrey", midpoint=0, name="ddG", limits=c(-1.3,2.1))+
  theme(axis.text.x = element_text(size = 10,
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 90),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("Position")+ylab("Mutation")+
  theme_bw()+theme(legend.key.width = unit(0.3,"cm"))
#make a row of average effects also
avg_effects_sm_nc_ddg <- lapply(split.data.frame(single_mut_df_nc, single_mut_df_nc$pos), FUN=function(x){median(x$mean_kcal.mol)})
avg_effects_sm_nc_ddg_df <- data.frame(pos=names(avg_effects_sm_nc_ddg), avg_dddg=as.vector(unlist(avg_effects_sm_nc_ddg)), mut=rep("Mean", length(avg_effects_sm_nc_ddg)))
avg_effects_sm_nc_ddg_df$pos <- factor(avg_effects_sm_nc_ddg_df$pos, levels=c("K-7", "N-6", "Y-5", "K-4", "Q-3", "T-2", "S-1", "V0"))
pheat_avg_sm_nc <- ggplot(avg_effects_sm_nc_ddg_df, aes(x=pos, y=mut, fill=avg_dddg))+
  geom_tile()+
  scale_fill_gradient2(low="orangered", mid="white", high="midnightblue", na.value = "darkgrey", midpoint=0, name="ddG", limits=c(-1.3,2.1))+
  theme(axis.text.x = element_text(size = 10,
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 90),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #geom_text(size=3)+
  xlab("Position")+
  theme_bw()+theme(legend.key.width = unit(0.3,"cm"))
  
#2nd order interactions (supplementary)
#want to make matrices for clustering (position-position) -- 2nd order coefficients (+ 1st order)
#first load the data and filter all the matrices to make more amenable to clustering
#second order ddgs
change_mut_numbers_dm <- function(x){paste(unlist(lapply(strsplit(x, "_"), FUN=function(x){unlist(lapply(x, change_mut_numbers, -8))})), collapse="_")}
#dm_conf_dgs <- all_dgs[intersect(which(nchar(all_dgs$id)>6), which(all_dgs$ci95_kcal.mol<1)),]
dm_conf_dgs <- all_dgs[which(nchar(all_dgs$id)>6),]
dm_conf_dgs$id <- unlist(lapply(dm_conf_dgs$id, change_mut_numbers_dm))
split_dm_conf_dgs_pos <- split.data.frame(dm_conf_dgs, dm_conf_dgs$Pos_ref)
make_submatrix <- function(x){
  cript2xa3_dm_df <- data.frame(x$mean)
  cript2xa3_dm_df$mut1 <- unlist(lapply(x$id, FUN=function(y){strsplit(y, "_")[[1]][1]}))
  cript2xa3_dm_df$mut2 <- unlist(lapply(x$id, FUN=function(y){strsplit(y, "_")[[1]][2]}))
  double_mut_m <- (dcast(cript2xa3_dm_df, mut1~mut2, value.var="x.mean"))
  rownames(double_mut_m) <- double_mut_m[,1]
  double_mut_m$mut1 <- NULL
  return(double_mut_m)
}
split_cast_dm_matrices <- lapply(split_dm_conf_dgs_pos, make_submatrix)
#for each matrix, pre-clustering, need to remove NA rows/cols ideally
filter_matrix_rows_na <- function(x){
  lengths_of_nas <- apply(x, 1, FUN=function(x){length(which(is.na(x)))/length(x)*100})
  if(length(which(lengths_of_nas >= 80))>0){
    dm_matrix_filt <- x[-which(lengths_of_nas >= 80),]
    return(dm_matrix_filt)
  }else{return(x)}
}
filter_matrix_cols_na <- function(x){
  lengths_of_nas <- apply(x, 2, FUN=function(x){length(which(is.na(x)))/length(x)*100})
  if(length(which(lengths_of_nas >= 80))>0){
    dm_matrix_filt <- x[,-which(lengths_of_nas >= 80)]
    return(dm_matrix_filt)
  }else{return(x)}
}
filter_matrix_cols_zero <- function(x){
  if(length(which(colSums(x, na.rm = T)==0))>0){
    dm_matrix_nozero <- x[,-which(colSums(x, na.rm = T)==0)]
    return(dm_matrix_nozero)}else{return(x)}
}
filter_matrix_rows_zero <- function(x){
  if(length(which(rowSums(x, na.rm = T)==0))>0){
    dm_matrix_nozero <- x[-which(rowSums(x, na.rm = T)==0),]
    return(dm_matrix_nozero)}else{return(x)}
}
#perform all filtering functions on each submatrix
filtered_cast_dm_matrices <- lapply(split_cast_dm_matrices, FUN=function(x){
  y <- filter_matrix_cols_na(x)
  y <- filter_matrix_rows_na(y)
  y <- filter_matrix_cols_zero(y)
  y <- filter_matrix_rows_zero(y)
  return(y)
})
#set up the color palette
paletteLength <- 50
myColor <- colorRampPalette(c("orangered", "white", "midnightblue"))(paletteLength)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
get_breaks_matrices <- function(x){
  myBreaks <- c(seq(min(x, na.rm=T), 0, length.out=ceiling(paletteLength/2)),
                seq(max(x, na.rm=T)/paletteLength, max(x, na.rm=T), length.out=floor(paletteLength/2)))
  return(myBreaks)
}

all_submatrix_breaks <- lapply(filtered_cast_dm_matrices, get_breaks_matrices)
#cluster all the individual matrices
###combine everything together
library(plyr)
all_mats <- rbind.fill.matrix(filtered_cast_dm_matrices)
breaks_for_all_mats <- get_breaks_matrices(all_mats)
#all_submatrix_clustered_maps_together <- list()
#for (i in 1:length(filtered_cast_dm_matrices)){
#  all_submatrix_clustered_maps_together[[i]] <- pheatmap(filtered_cast_dm_matrices[[i]], color=myColor, breaks = breaks_for_all_mats)
#}
###threshold the scale
all_mats_thresholded <- all_mats
max_neg_value <- max(all_mats_thresholded[which(-1.5 >=all_mats_thresholded)])
max_pos_value <- max(all_mats_thresholded[which(all_mats_thresholded >= 1.5)])
all_mats_thresholded[which(max_neg_value >=all_mats_thresholded)] <- max_neg_value
all_mats_thresholded[which(all_mats_thresholded >= max_pos_value)] <- max_pos_value
breaks_for_all_mats_thresh <- get_breaks_matrices(all_mats_thresholded)

#also putting the single mutant values on each axis in order to reference each double mut
sm_conf_dgs_all <- all_dgs[which(nchar(all_dgs$id)<6 & nchar(all_dgs$id)>2),]
sm_conf_dgs_all$id <- unlist(lapply(sm_conf_dgs_all$id, change_mut_numbers, -8))
#make annotation row/col for each submatrix
make_annot_col <- function(x){
  annot_col_df <- data.frame(dgs=sm_conf_dgs_all$mean[which(sm_conf_dgs_all$id %in% colnames(x))])
  rownames(annot_col_df) <- sm_conf_dgs_all$id[which(sm_conf_dgs_all$id %in% colnames(x))]
  return(annot_col_df)
}
make_annot_row <- function(x){
  annot_row_df <- data.frame(dgs=sm_conf_dgs_all$mean[which(sm_conf_dgs_all$id %in% rownames(x))])
  rownames(annot_row_df) <- sm_conf_dgs_all$id[which(sm_conf_dgs_all$id %in% rownames(x))]
  return(annot_row_df)
}
library(ComplexHeatmap)
library(circlize)
#clustering where everything is thresholded 
all_submatrix_clustered_maps_together_thresh <- list()
for (i in 1:length(filtered_cast_dm_matrices)){
  all_submatrix_clustered_maps_together_thresh[[i]] <- Heatmap(filtered_cast_dm_matrices[[i]], 
                                                               col=colorRamp2(breaks_for_all_mats_thresh, myColor),
                                                               row_names_gp = grid::gpar(fontsize=8),
                                                               column_names_gp = grid::gpar(fontsize = 8),
                                                               width = ncol(filtered_cast_dm_matrices[[i]])*unit(3, "mm"), 
                                                               height = nrow(filtered_cast_dm_matrices[[i]])*unit(3, "mm"),
                                                               name = "dddG", 
                                                               top_annotation = HeatmapAnnotation(sm=make_annot_col(filtered_cast_dm_matrices[[i]])$dgs, col=list(sm=colorRamp2(breaks_for_all_mats_thresh, myColor)), show_legend=F, show_annotation_name = F),
                                                               right_annotation = rowAnnotation(sm=make_annot_row(filtered_cast_dm_matrices[[i]])$dgs, col=list(sm=colorRamp2(breaks_for_all_mats_thresh, myColor)), show_legend=F, show_annotation_name=F))
}
