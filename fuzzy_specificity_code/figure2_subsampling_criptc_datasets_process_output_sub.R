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
#processing results from 10 subsampled cript c combinatorial mochi runs
#get the files
home_dir <- getwd()
source("./fuzzy_specificity_code/source_functions/basic_functions_fuzzy_specificity_TZ.R")
setwd("./fuzzy_specificity_data/mochi_results/")
files_w_mochi_results <- list.files(recursive=T)
files_w_mochi_results_subsamp_linear <- files_w_mochi_results[grep("cript_nc_balanced_exp10_linear_[0-9]*", files_w_mochi_results)]
files_w_mochi_results_subsamp_1 <- files_w_mochi_results[grep("cript_nc_balanced_exp10_1_[0-9]*", files_w_mochi_results)]
files_w_mochi_results_subsamp_2 <- files_w_mochi_results[grep("cript_nc_balanced_exp10_2_[0-9]*", files_w_mochi_results)]
#get predictions files
subsamp_pred_files_linear <- files_w_mochi_results_subsamp_linear[grep("predicted_phenotypes_all.txt", files_w_mochi_results_subsamp_linear)]
subsamp_pred_files_nc_1 <- files_w_mochi_results_subsamp_1[grep("predicted_phenotypes_all.txt", files_w_mochi_results_subsamp_1)]
subsamp_pred_files_nc_12 <- files_w_mochi_results_subsamp_2[grep("predicted_phenotypes_all.txt", files_w_mochi_results_subsamp_2)]
subsamp_dg_files_linear <- files_w_mochi_results_subsamp_linear[grep("weights_Folding.txt", files_w_mochi_results_subsamp_linear)]
subsamp_dg_files_nc_1 <- files_w_mochi_results_subsamp_1[grep("weights_Folding.txt", files_w_mochi_results_subsamp_1)]
subsamp_dg_files_nc_12 <- files_w_mochi_results_subsamp_2[grep("weights_Folding.txt", files_w_mochi_results_subsamp_2)]
#read in files
subsamp_pred_files_linear_read <- lapply(subsamp_pred_files_linear, combine_read_delim_heldout_pred)
subsamp_pred_files_nc_1_read <- lapply(subsamp_pred_files_nc_1, combine_read_delim_heldout_pred)
subsamp_pred_files_nc_12_read <- lapply(subsamp_pred_files_nc_12, combine_read_delim_heldout_pred)
subsamp_dg_files_linear_read <- lapply(subsamp_dg_files_linear, read.delim)
subsamp_dg_files_nc_1_read <- lapply(subsamp_dg_files_nc_1, read.delim)
subsamp_dg_files_nc_12_read <- lapply(subsamp_dg_files_nc_12, read.delim)
#function to plot everything
make_plots_compare_performance <- function(cript_nc_1_linear, cript_nc_1, cript_nc_12, main_supp){
  mochi_obspred_linear_order_1_c <- make_gg_scatter_with_xy_line(cript_nc_1_linear$fitness[which(cript_nc_1_linear$criptc_combo==1)], cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptc_combo==1)], "Observed fitness", "Predicted fitness")
  mochi_obspred_linear_order_1_n <- make_gg_scatter_with_xy_line(cript_nc_1_linear$fitness[which(cript_nc_1_linear$criptn_combo==1)], cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptn_combo==1)], "Observed fitness", "Predicted fitness")
  #want to plot the residuals to these fits to quantify bias
  resids_linear_order_1_c <- make_gg_scatterplot_no_corr(cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptc_combo==1)], cript_nc_1_linear$fitness[which(cript_nc_1_linear$criptc_combo==1)] - cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptc_combo==1)], "Predicted fitness", "Residual")+
    geom_hline(yintercept=0,color="darkgrey", linewidth=0.9)
  resids_linear_order_1_n <- make_gg_scatterplot_no_corr(cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptn_combo==1)], cript_nc_1_linear$fitness[which(cript_nc_1_linear$criptn_combo==1)] - cript_nc_1_linear$predicted_phen[which(cript_nc_1_linear$criptn_combo==1)], "Predicted fitness", "Residual")+
    geom_hline(yintercept=0,color="darkgrey", linewidth=0.9)
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
  predictions_main <- ggarrange(mochi_obspred_order1_n+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),#+ ggtitle("CRIPT N MoCHI 1st order"),
                                mochi_obspred_order1_c+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),#+ggtitle("CRIPT C MoCHI 1st order"), 
                                mochi_obs_dddg_order1_n+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),#+ggtitle("CRIPT N MoCHI 2nd order terms"),
                                mochi_obs_dddg_order1_c+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),#+ggtitle("CRIPT C MoCHI 2nd order terms"),
                                ncol=2, nrow=2, legend="right")
  predictions_supp <- ggarrange(mochi_obspred_linear_order_1_n+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),
                                mochi_obspred_linear_order_1_c+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),
                                resids_linear_order_1_n+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),
                                resids_linear_order_1_c+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),
                                resids_mochi_order_1_n+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),
                                resids_mochi_order_1_c+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),
                                mochi_obspred_order12_n+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),#+ggtitle("CRIPT N 2nd order"),
                                mochi_obspred_order12_c+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),#+ggtitle("CRIPT C 2nd order"),
                                mochi_obs_dddg_order12_n+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),#+ggtitle("CRIPT N MoCHI 2nd order terms"),
                                mochi_obs_dddg_order12_c+theme(legend.key.width = unit(0.3,"cm"), legend.key.height = unit(0.4,"cm")),#+ggtitle("CRIPT C MoCHI 2nd order terms"), 
                                ncol=2, nrow=5, legend="right")
  if(main_supp=="print_main"){
    return(predictions_main)
  }else{
    return(predictions_supp)
    }
}
ten_subsamp_plots_main <- NULL
for (i in 1:length(subsamp_pred_files_nc_1_read)){
  ten_subsamp_plots_main[[i]] <- make_plots_compare_performance(subsamp_pred_files_linear_read[[i]],
                                                                subsamp_pred_files_nc_1_read[[i]],
                                                                subsamp_pred_files_nc_12_read[[i]],
                                                                "print_main")
}
ten_subsamp_plots_supp <- NULL
for (i in 1:length(subsamp_pred_files_nc_1_read)){
  ten_subsamp_plots_supp[[i]] <- make_plots_compare_performance(subsamp_pred_files_linear_read[[i]],
                                                                subsamp_pred_files_nc_1_read[[i]],
                                                                subsamp_pred_files_nc_12_read[[i]],
                                                                "print_supp")
}

get_corrs_compare_performance <- function(cript_nc_1, n_or_c){
  if(n_or_c==1){
    mochi_obspred_order1_c <- cor.test(cript_nc_1$fitness[which(cript_nc_1$criptc_combo==1)], cript_nc_1$predicted_phen[which(cript_nc_1$criptc_combo==1)])
  }else{if(n_or_c==2){
    mochi_obspred_order1_n <- cor.test(cript_nc_1$fitness[which(cript_nc_1$criptn_combo==1)], cript_nc_1$predicted_phen[which(cript_nc_1$criptn_combo==1)])
  }}
}
ten_corrs_c_linear <- NULL
for (i in 1:length(subsamp_pred_files_linear_read)){
  ten_corrs_c_linear[[i]] <- (get_corrs_compare_performance(subsamp_pred_files_linear_read[[i]],1)$estimate)^2
}
ten_corrs_c_1 <- NULL
for (i in 1:length(subsamp_pred_files_nc_1_read)){
  ten_corrs_c_1[[i]] <- (get_corrs_compare_performance(subsamp_pred_files_nc_1_read[[i]],1)$estimate)^2
}
ten_corrs_c_2 <- NULL
for (i in 1:length(subsamp_pred_files_nc_12_read)){
  ten_corrs_c_2[[i]] <- (get_corrs_compare_performance(subsamp_pred_files_nc_12_read[[i]],1)$estimate)^2
}
ten_corrs_n_linear <- NULL
for (i in 1:length(subsamp_pred_files_linear_read)){
  ten_corrs_n_linear[[i]] <- (get_corrs_compare_performance(subsamp_pred_files_linear_read[[i]],2)$estimate)^2
}
ten_corrs_n_1 <- NULL
for (i in 1:length(subsamp_pred_files_nc_1_read)){
  ten_corrs_n_1[[i]] <- (get_corrs_compare_performance(subsamp_pred_files_nc_1_read[[i]],2)$estimate)^2
}
ten_corrs_n_2 <- NULL
for (i in 1:length(subsamp_pred_files_nc_12_read)){
  ten_corrs_n_2[[i]] <- (get_corrs_compare_performance(subsamp_pred_files_nc_12_read[[i]],2)$estimate)^2
}
#put it into a df
summary_corrs_all <- data.frame(cript_n_linear = unlist(ten_corrs_n_linear), 
                                cript_n_1 = unlist(ten_corrs_n_1), 
                                cript_n_2 = unlist(ten_corrs_n_2),
                                cript_c_linear = unlist(ten_corrs_c_linear), 
                                cript_c_1 = unlist(ten_corrs_c_1), 
                                cript_c_2 = unlist(ten_corrs_c_2))
summary_corrs_all_melted <- reshape2::melt(summary_corrs_all)
summary_plot_corrs <- ggplot(summary_corrs_all_melted, aes(variable, value))+
  xlab("MoCHI model")+
  ylab(bquote(R^2))
summary_plot_corrs <- summary_plot_corrs + geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.2)+
  theme_bw()
setwd(home_dir)





