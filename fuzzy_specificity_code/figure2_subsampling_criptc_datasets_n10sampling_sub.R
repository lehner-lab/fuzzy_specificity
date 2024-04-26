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
#this script creates 10 iterations of sampling the distribution of cript combinatorial c variants such that alive variants are heavily favoured
#read in data, find dead mode threshold
criptc_for_mochi <- read.delim("C:/Users/tzarin/OneDrive - CRG - Centre de Regulacio Genomica/2024-03_pdzcript_figs_scripts/criptc_norm_for_mochi.tsv")
criptc_dead_mode_threshold <- mean(criptc_for_mochi$fitness[which(criptc_for_mochi$STOP==T)])+2*sd(criptc_for_mochi$fitness[which(criptc_for_mochi$STOP==T)])
aliven_downsample <- 2*(length(which(criptc_for_mochi$fitness>criptc_dead_mode_threshold)))
criptc_dead_mode_peak <- mean(criptc_for_mochi$fitness[which(criptc_for_mochi$STOP==T)])
distance_to_dead_mode_peak <- 10^10^(1 - criptc_for_mochi$fitness/criptc_dead_mode_peak)
#want to scale this between 0 and 1
distance_to_dead_mode_peak_rescaled <- distance_to_dead_mode_peak - min(distance_to_dead_mode_peak) / max(distance_to_dead_mode_peak) - min(distance_to_dead_mode_peak)
prob_distance_to_dead_mode_peak <- distance_to_dead_mode_peak_rescaled/sum(distance_to_dead_mode_peak_rescaled)
#make the new dataframes
set.seed(905)
subsampled_variants_aliven <- NULL
subsampled_variants_aliven <- replicate(10,sample(1:nrow(criptc_for_mochi), size = aliven_downsample, replace=F, prob = distance_to_dead_mode_peak))
###############################
###initially only sampled once -- to reproduce results in figures from one sampling iteration: 
#set.seed(905)
#subsampled_variants_aliven <- criptc_for_mochi[sample(1:nrow(criptc_for_mochi), size = aliven_downsample, replace=F, prob = distance_to_dead_mode_peak),]
#subsampled_variants_aliven <- make_wt_col_na_free(subsampled_variants_aliven)
#write.table(subsampled_variants_aliven, "./fuzzy_specificity_data/figure2_energetic_landscape/subsampled_variants_aliven_10.tsv", quote=F, sep="\t", row.names=F)
###############################
criptc_subset_dfs <- NULL
for (i in 1:ncol(subsampled_variants_aliven)){
  criptc_subset_dfs[[i]] <- criptc_for_mochi[subsampled_variants_aliven[,i],]
}
#check how the sampling affects distribution of variants
#make a graph of pre- and post- subsampled/balanced datasets
hist_pre_subsampling <- ggplot(criptc_for_mochi, aes(x=fitness))+
  geom_histogram(fill="white", colour="black", bins=100)+
  theme_bw()+
  xlab("CRIPT C normalized fitness (all variants)")+
  ggtitle(paste("n = ", nrow(criptc_for_mochi), " | alive = ", round(length(which(criptc_for_mochi$fitness>criptc_dead_mode_threshold))/nrow(criptc_for_mochi)*100,2), "%",sep=""))
hist_pre_subsampling <- hist_pre_subsampling+geom_vline(xintercept=criptc_dead_mode_threshold, col=custom_scale_discrete_1[4])

make_plot_subsampled_fitness <- function(x){
  hist_post_subsampling <- ggplot(x, aes(x=fitness))+
    geom_histogram(fill="white", colour="black", bins=100)+
    theme_bw()+
    xlab("CRIPT C normalized fitness (subsampled)")+
    ggtitle(paste("n = ", nrow(x), " | alive = ", round(length(which(x$fitness>criptc_dead_mode_threshold))/nrow(x)*100,2), "%", sep=""))
  hist_post_subsampling <- hist_post_subsampling+geom_vline(xintercept=criptc_dead_mode_threshold, col=custom_scale_discrete_1[4])
  hist_post_subsampling
}
plots_subsampled <- lapply(criptc_subset_dfs, make_plot_subsampled_fitness)
subsampled_data_criptc <- cowplot::plot_grid(plotlist=plots_subsampled, nrow=5, ncol=2, align = "hv")


subsampled_variants_aliven_dfs <- NULL
for (i in 1:10){
  subsampled_variants_aliven_dfs[[i]] <- make_wt_col_na_free(criptc_subset_dfs[[i]])
  write.table(subsampled_variants_aliven_dfs[[i]], paste("./fuzzy_specificity_data/figures/figure2_energetic_landscape/subsampled_variants_aliven_", i, ".tsv", sep=""), quote=F, sep="\t", row.names=F)
}



