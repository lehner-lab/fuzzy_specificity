################################################################################################
#Source code for Figure3-5 (+ related supplementary material) from 
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
#this script does all data analysis related to pdz3-cript double mutant library
wt_seq <- "GEEDIPREPRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQYKPEEYSRFEAKNYKQTSV"
pdz_dimsum <- read.delim("./fuzzy_specificity_data/normalized_block1_block2_pdz_df_block2wt.tsv")
#for plots where block1 and block2 fitness shown together, make fitness col = normalized fitness
pdz_dimsum_fitness2norm <- pdz_dimsum
pdz_dimsum_fitness2norm$fitness <- pdz_dimsum_fitness2norm$nor_fitness
repcorr_mat <- get_rep_corr_matrix(pdz_dimsum_fitness2norm)
pdz_cript_fitcorr_mat <- ggpairs(pdz_dimsum_fitness2norm[,grep("fitness[0-9]",colnames(pdz_dimsum_fitness2norm))],
                         lower = list(continuous = "hexbin", combo = "facethist", discrete = "facetbar", na =  "na"))+theme_bw()
#looking at nucleotide variants (for supp)
b1_variant_counts <- read.delim("./fuzzy_specificity_data/stage_4_dimsum_pdz3_input/pdz3b1_short_long_consolidated_variant_counts_wtfiltered.txt")
b2_variant_counts <- read.delim("./fuzzy_specificity_data/stage_4_dimsum_pdz3_input/pdz3b2_short_long_consolidated_variant_counts_wtfiltered.txt")
#counting the number of nuc substitutions in ones that have at least one rep > 10 reads
b1_design <- read.table("./fuzzy_specificity_data/bc_files/pdz3b1_bcid_long.txt", header=T)
b2_design <- read.table("./fuzzy_specificity_data/bc_files/pdz3b2_bcid_long.txt", header=T)
##
wt_b1_full <- paste(pdz_b1_nuc_wt_full, cript_nuc_full, sep="")
wt_b2_full <- paste(pdz_b2_nuc_wt_full, cript_nuc_full, sep="")
##
#for each variant with at least 1 count in each replicate, how many substitutions away is it from wt?
b1_nuc_hamming <- lapply(b1_variant_counts$nt_seq[index_count_cols_input(b1_variant_counts)], list.string.diff, wt_b1_full)
b1_nuc_hamming_num_diffs <- unlist(lapply(b1_nuc_hamming, nrow))
b2_nuc_hamming <- lapply(b2_variant_counts$nt_seq[index_count_cols_input(b2_variant_counts)], list.string.diff, wt_b2_full)
b2_nuc_hamming_num_diffs <- unlist(lapply(b2_nuc_hamming, nrow))
b1_design_nuc_hamming <- lapply(tolower(b1_design$barcode), list.string.diff, wt_b1_full)
b1_design_nuc_hamming_diffs <- unlist(lapply(b1_design_nuc_hamming, nrow))
b2_design_nuc_hamming <- lapply(tolower(b2_design$barcode), list.string.diff, wt_b2_full)
b2_design_nuc_hamming_diffs <- unlist(lapply(b2_design_nuc_hamming, nrow))
#add in "1" distance for b1/b2 nuc hamming 
b1_nuc_hamming_num_diffs_table <- table(b1_nuc_hamming_num_diffs)
b1_nuc_hamming_num_diffs_table <- c(b1_nuc_hamming_num_diffs_table, 0)
names(b1_nuc_hamming_num_diffs_table)[7] <- "1"
#same for b2
b2_nuc_hamming_num_diffs_table <- table(b2_nuc_hamming_num_diffs)
b2_nuc_hamming_num_diffs_table <- c(b2_nuc_hamming_num_diffs_table, 0)
names(b2_nuc_hamming_num_diffs_table)[7] <- "1"
#actual libraries >10 in input:
b1_b2_nuc_hamming <- rbind(b1_nuc_hamming_num_diffs_table, b2_nuc_hamming_num_diffs_table)[,c("0", "1", "2", "3", "4", "5", "6")]
#designed libraries:
b1_designed_table <- table(b1_design_nuc_hamming_diffs)
b2_designed_table <- table(b2_design_nuc_hamming_diffs)
designed_b1b2_hamming <- rbind(b1_designed_table, b2_designed_table)
all_nuc_hamming <- rbind(b1_b2_nuc_hamming, designed_b1b2_hamming)
all_nuc_hamming_melted <- reshape2::melt(all_nuc_hamming)
###
nuc_ham_libraries <- ggplot(data=all_nuc_hamming_melted, aes(x=Var2, y=value, fill=Var1))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_x_continuous(breaks=seq(0,6,1))+
  scale_fill_manual(name=" ", values=custom_scale_discrete_1, breaks=c("b1_designed_table", "b1_nuc_hamming_num_diffs_table", "b2_designed_table", "b2_nuc_hamming_num_diffs_table"), labels=c("B1_designed", "B1_seq_>10", "B2_designed", "B2_seq_>10"))+
  theme_bw()+
  xlab("Number of nucleotide substitutions")+
  ylab("Number of variants")
#############
densityplot_pdz_cript <- make_gg_density_nor_fit_general(pdz_dimsum_fitness2norm)
aaham_pdz_cript <- make_aaham_plot(pdz_dimsum_fitness2norm)
#dead variant cutoff and justification
mean_stops <- mean(pdz_dimsum_fitness2norm$fitness[which(pdz_dimsum_fitness2norm$STOP==T)])
mode_stops <- as.numeric(names(table(pdz_dimsum_fitness2norm$fitness[which(pdz_dimsum_fitness2norm$STOP==T)])[order(table(pdz_dimsum_fitness2norm$fitness[which(pdz_dimsum_fitness2norm$STOP==T)]), decreasing=T)][1]))
median_stops <- median(pdz_dimsum_fitness2norm$fitness[which(pdz_dimsum_fitness2norm$STOP==T)])
dead_var_cutoff <- mean(pdz_dimsum_fitness2norm$fitness[which(pdz_dimsum_fitness2norm$STOP==T)])+2*sd(pdz_dimsum_fitness2norm$fitness[which(pdz_dimsum_fitness2norm$STOP==T)])
dead_threshold_hist <- ggplot(pdz_dimsum_fitness2norm[which(pdz_dimsum_fitness2norm$STOP==T),], aes(x=fitness))+
  geom_histogram(fill="white", colour="black", bins=50)+
  theme_bw()+
  xlab("Normalized fitness STOP variants")
dead_threshold_hist <- dead_threshold_hist+
  geom_vline(xintercept=dead_var_cutoff, col=custom_scale_discrete_1[4])+
  geom_vline(xintercept=mean_stops, col=custom_scale_discrete_1[3], linetype="dashed")
#single mutant fitness heatmap with adjusted position coordinates
#separating the pdz from cript
pdz_dimsum_wtcript <- pdz_dimsum_fitness2norm[which(unlist(lapply(pdz_dimsum_fitness2norm$aa_seq, FUN=str_sub, -8,-1)) %in% str_sub(wt_seq,-8,-1)),]
pdz_dimsum_wtpdz <- pdz_dimsum_fitness2norm[which(unlist(lapply(pdz_dimsum_fitness2norm$aa_seq, FUN=str_sub, 1,100)) %in% str_sub(wt_seq,1,100)),]
#because cript single mutant effects is across 2 blocks and variants are repeated, need to collapse
#look at block1 and block2 separately - different?
pheat_sm_fitness_pdz3cript_pdz <- make_pheat_fitness_sm(pdz_dimsum_wtcript, 302)+
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.key.height = unit(0.3, "cm"), legend.key.width=unit(1, "cm"))+
  scale_fill_gradientn(breaks=c(-1,0,1), colors=c("midnightblue", "white", "orangered"), limits=c(-1.3,1.3), name="Normalized fitness")
#before consolidating make sure wt are correlated
matched_wt_pdz_block1_block2 <- match_datasets(pdz_dimsum_wtpdz[which(pdz_dimsum_wtpdz$block=="block1"),], pdz_dimsum_wtpdz[which(pdz_dimsum_wtpdz$block=="block2"),])
block1_block2_wt_scatter <- make_gg_scatter_with_xy_line(matched_wt_pdz_block1_block2$x_fitness, matched_wt_pdz_block1_block2$y_fitness, "Block 1 WT PDZ", "Block 2 WT PDZ")
pdz_dimsum_wtpdz_consolidated <- pdz_dimsum_wtpdz %>%
  group_by(aa_seq) %>%
  dplyr::summarise(WT=unique(WT), Nham_aa=unique(Nham_aa), fitness_med=median(fitness)) %>%
  as.data.frame
colnames(pdz_dimsum_wtpdz_consolidated)[which(colnames(pdz_dimsum_wtpdz_consolidated)=="fitness_med")] <- "fitness"
make_pheat_fitness_sm(pdz_dimsum_wtpdz_consolidated, -108)
#make the cript plot
pheat_sm_fitness_pdz3cript_cript <- make_pheat_fitness_sm(pdz_dimsum_wtpdz_consolidated, -108)+
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.key.height = unit(0.3, "cm"), legend.key.width = unit(1, "cm"))+
  scale_fill_gradientn(breaks=c(-1,0,1), colors=c("midnightblue", "white", "orangered"), limits=c(-1.3,1.3), name="Normalized fitness")
pdz_cript_sm_together <- cowplot::plot_grid(pheat_sm_fitness_pdz3cript_cript, pheat_sm_fitness_pdz3cript_pdz, nrow=1, rel_widths = c(8, 50), rel_heights = c(1,1), align = "hv")
##################################################
#plot how the average effect of mutation is at each position (binding fitness on structure)
matrix_sm_fitness_pdz3cript_cript <- make_matrix_fitness_sm(pdz_dimsum_wtpdz_consolidated, -108)
matrix_sm_fitness_pdz3cript_cript_minus_stops <- matrix_sm_fitness_pdz3cript_cript[-which(matrix_sm_fitness_pdz3cript_cript$mut=="*"),]
matrix_sm_fitness_pdz3cript_cript_minus_stops_average <- matrix_sm_fitness_pdz3cript_cript_minus_stops %>%
  group_by(pos) %>%
  dplyr::summarise(avg_fitness=mean(fitness), avg_pos=unique(pos)) %>%
  as.data.frame
eg_setatt_file <-read.delim("./fuzzy_specificity_data/eg_files/percentExposed.txt")
home_dir <- getwd()
setwd("./fuzzy_specificity_data/figures/figure3_pdzcript_setup/")
setatt_df_cript <- data.frame(empty_col=rep("", nrow(matrix_sm_fitness_pdz3cript_cript_minus_stops_average)), pos=paste("/B:", 2:9, sep=""), att=matrix_sm_fitness_pdz3cript_cript_minus_stops_average$avg_fitness)
setatt_df_cript_table <- write.table(setatt_df_cript, "criptavgsmfitness.defattr", quote=F, sep="\t",row.names=F, col.names=F, eol = "\r")
content_to_change <- readLines("criptavgsmfitness.defattr")
new_content <- c("attribute: criptavgsmfitness","match mode: 1-to-1", "recipient: residues", content_to_change)
writeLines(new_content, "criptavgsmfitness.defattr", sep="\n")
#same as above, now for pdz
#plot how the average effect of mutation is at each position (binding fitness on structure)
matrix_sm_fitness_pdz3cript_pdz3 <- make_matrix_fitness_sm(pdz_dimsum_wtcript, 302)
matrix_sm_fitness_pdz3cript_pdz3_minus_stops <- matrix_sm_fitness_pdz3cript_pdz3[-which(matrix_sm_fitness_pdz3cript_pdz3$mut=="*"),]
matrix_sm_fitness_pdz3cript_pdz3_minus_stops_average <- matrix_sm_fitness_pdz3cript_pdz3_minus_stops %>%
  group_by(pos) %>%
  dplyr::summarise(avg_fitness=mean(fitness), avg_pos=unique(pos)) %>%
  as.data.frame
min_pos <- min(as.numeric(str_sub(unique(matrix_sm_fitness_pdz3cript_pdz3$pos), 2, 4)))
max_pos <- max(as.numeric(str_sub(unique(matrix_sm_fitness_pdz3cript_pdz3$pos), 2, 4)))
setatt_df_pdz3 <- data.frame(empty_col=rep("", nrow(matrix_sm_fitness_pdz3cript_pdz3_minus_stops_average)), pos=paste("/A:", min_pos:max_pos , sep=""), att=matrix_sm_fitness_pdz3cript_pdz3_minus_stops_average$avg_fitness)
setatt_df_pdz3_table <- write.table(setatt_df_pdz3, "pdz3avgsmfitness.defattr", quote=F, sep="\t",row.names=F, col.names=F, eol = "\r")
content_to_change <- readLines("pdz3avgsmfitness.defattr")
new_content <- c("attribute: pdz3avgsmfitness","match mode: 1-to-1", "recipient: residues", content_to_change)
writeLines(new_content, "pdz3avgsmfitness.defattr", sep="\n")
setwd(home_dir)
#color byattribute criptavgsmfitness palette -0.8,#5555ff:0,white:0.8,orange
#look at effect of mutations with hamming distance but also cript vs pdz if applicable
pdz_dimsum_fitness2norm$Nham_aa <- as.factor(pdz_dimsum_fitness2norm$Nham_aa)
pdz_dimsum_aaham <- ggplot(pdz_dimsum_fitness2norm[which(pdz_dimsum_fitness2norm$Nham_aa %in% c(0,1,2,3,4)),], aes(x=Nham_aa, y=fitness)) + 
  geom_violin(draw_quantiles = 0.5) + 
  theme_bw() +
  xlab("AA Hamming distance")+
  ylab("Fitness")+
  theme(legend.title = element_blank())
########################################################################################################################################################
###evaluate how model is performing
#take out the wt and use only held-out fold prediction rather than mean
data_path <- "./fuzzy_specificity_data/mochi_results/pdz_b1b2/task_1/predictions/predicted_phenotypes_all.txt"
pdz_merged <- combine_read_delim_heldout_pred(data_path)
mod_perf_pdzcript <- make_gg_scatter_with_xy_line(pdz_merged$fitness, pdz_merged$predicted_phen, "observed pdz_merged", "predicted")
addtrait_pdzcript <- make_gg_scatterplot(apply(pdz_merged[,49:58], 1, mean), pdz_merged$fitness, "additive trait pdz merged", "fitness")
###read in weights for first order model, compare to n/c first order weights 
pdz_merged_weights <- read.delim("./fuzzy_specificity_data/mochi_results/pdz_b1b2/task_1/weights/weights_Folding.txt")
pdz_merged_weights_idchange <- pdz_merged_weights
pdz_merged_weights_idchange$id <- unlist(lapply(pdz_merged_weights_idchange$id, change_mut_numbers, -100))
criptnc_weights <- read.delim("./fuzzy_specificity_data/mochi_results/cript_nc_balanced_exp10_1/task_1/weights/weights_Folding.txt")
matched_ddgs_pdzcript_criptcombo <- match_datasets_dgs(pdz_merged_weights_idchange, criptnc_weights)
matched_ddgs_pdzcript_criptcombo_errors <- make_scatterplot_errorbars(matched_ddgs_pdzcript_criptcombo, matched_ddgs_pdzcript_criptcombo$x_mean_kcal.mol, matched_ddgs_pdzcript_criptcombo$y_mean_kcal.mol, lbl = " ", error_y = matched_ddgs_pdzcript_criptcombo$y_95ci, error_x = matched_ddgs_pdzcript_criptcombo$x_95ci, xname = "pdz3cript trans ddGs", yname="combinatorial ddGs") + 
  ggtitle(bquote(R^2 ~ .(paste(round((cor(matched_ddgs_pdzcript_criptcombo$x_mean_kcal.mol, matched_ddgs_pdzcript_criptcombo$y_mean_kcal.mol)^2), digits=3), " | N =", length(matched_ddgs_pdzcript_criptcombo$x_mean_kcal.mol), sep=" "))))
#pdz annotations
pdz3cript_mochi <- read.delim(data_path)
#wanna make a heatmap with the residuals
#get residuals
pdz3cript_resids <- pdz3cript_mochi$fitness - pdz3cript_mochi$mean
#make a map of residuals vs fitness
resids_pdzcript <- make_gg_scatterplot(pdz3cript_mochi$mean, pdz3cript_resids, "Predicted fitness", "Residual")
pdz3cript_double_mut_pos <- which(pdz3cript_mochi$Nham_aa==2)
pdz3cript_dm <- cbind(pdz3cript_mochi[pdz3cript_double_mut_pos,], pdz3cript_resids[pdz3cript_double_mut_pos])
#compute double mutant positions/ids
#note: this part takes a few mins
pdz3cript_dm_diffs <- lapply(pdz3cript_dm$aa_seq, FUN=list.string.diff, wt_seq)
pdz3cript_dm_df_mut1 <- NULL
for (i in 1:length(pdz3cript_dm_diffs)){
  pdz3cript_dm_df_mut1[i] <- paste(pdz3cript_dm_diffs[[i]][1,]$poly.seq.b, pdz3cript_dm_diffs[[i]][1,]$position, pdz3cript_dm_diffs[[i]][1,]$poly.seq.a, sep="")
}
pdz3cript_dm_df_mut2 <- NULL
for (i in 1:length(pdz3cript_dm_diffs)){
  pdz3cript_dm_df_mut2[i] <- paste(pdz3cript_dm_diffs[[i]][2,]$poly.seq.b, pdz3cript_dm_diffs[[i]][2,]$position, pdz3cript_dm_diffs[[i]][2,]$poly.seq.a, sep="")
}
#analysis of double mutants
pdz3cript_dm$pdz_mut1 <- pdz3cript_dm_df_mut1
pdz3cript_dm$pdz_mut1_realprot_coord <- unlist(lapply(pdz3cript_dm$pdz_mut1, change_mut_numbers, 302))
pdz3cript_dm$cript_mut2 <- pdz3cript_dm_df_mut2
pdz3cript_dm$cript_mut2_realprot_coord <- unlist(lapply(pdz3cript_dm$cript_mut2, change_mut_numbers, -108))
pdz3cript_dm$pos_pairs <- paste(unlist(lapply(pdz3cript_dm_df_mut1, get_mut_pos)), unlist(lapply(pdz3cript_dm_df_mut2, get_mut_pos)), sep="_")
pdz3cript_dm$pos_pairs_realprot_coord <- paste(as.numeric(unlist(lapply(pdz3cript_dm_df_mut1, get_mut_pos))+302), as.numeric(unlist(lapply(pdz3cript_dm_df_mut2, get_mut_pos))-108), sep="_")
pdz3cript_dm$mut1_mut2 <- paste(pdz3cript_dm$pdz_mut1, pdz3cript_dm$cript_mut2, sep="_")
pdz3cript_dm$mut1_mut2_realprot_coord <- paste(pdz3cript_dm$pdz_mut1_realprot_coord, pdz3cript_dm$cript_mut2_realprot_coord, sep="_")
pdz3cript_dm$type_cript <- rep("", nrow(pdz3cript_dm))
pdz3cript_dm$type_cript[which(unlist(lapply(pdz3cript_dm$cript_mut2, get_mut_pos)) %in% c(101:104))] <- "N"
pdz3cript_dm$type_cript[which(unlist(lapply(pdz3cript_dm$cript_mut2, get_mut_pos)) %in% c(105:108))] <- "C"
#stats to see which residuals are significant
pdz3cript_dm$zscore_dm <- (pdz3cript_dm$`pdz3cript_resids[pdz3cript_double_mut_pos]`)/sqrt((pdz3cript_dm$sigma)^2)
pdz3cript_dm$ztest_dm <- 2*pnorm(-abs(pdz3cript_dm$zscore_dm)) 
pdz3cript_dm$fdr_dm <- p.adjust(pdz3cript_dm$ztest_dm, method="fdr", n=length(pdz3cript_dm$ztest_dm))
pdz3cript_dm$passes_fdr <- rep(0, nrow(pdz3cript_dm))
pdz3cript_dm$passes_fdr[which(pdz3cript_dm$fdr_dm<0.1)] <- 1
#plot the residuals
sc_plot_resids <- ggplot(pdz3cript_dm, aes(x=fitness, y=mean, fill=`pdz3cript_resids[pdz3cript_double_mut_pos]`))+
  geom_point(size=1.5, shape=21)+
  scale_fill_viridis()+
  theme(element_text(size=12, family="serif"))+
  xlab("Observed fitness")+
  ylab("Predicted fitness")+
  labs(fill="Residual")+
  theme_bw()+theme(legend.key.width = unit(0.3, "cm"))+geom_density2d(color="white")
#plot the FDR of residuals
sc_plot_fdr <- ggplot(pdz3cript_dm, aes(x=fitness, y=mean, fill=fdr_dm*100))+
  geom_point(size=1.5, shape=21)+
  scale_fill_viridis_c(direction=-1, option="B", begin= 0.1, end = 0.9)+
  theme(element_text(size=12, family="serif"))+
  xlab("Observed fitness")+
  ylab("Predicted fitness")+
  labs(fill="FDR of residual (%)")+
  theme_bw()+theme(legend.key.width = unit(0.3, "cm"))+geom_density2d(color="white")+
  geom_vline(xintercept=dead_var_cutoff, color="grey", linetype="dashed", linewidth=1.5)
###
sc_plots_resid_fdr <- plot_grid(sc_plot_resids, sc_plot_fdr, nrow=2, ncol=1, align = "hv") 
#also want the z-score distribution
hist_zscore_residuals <- ggplot(pdz3cript_dm, aes(x=zscore_dm))+
  geom_histogram(fill="white", colour="black")+
  theme_bw()+
  xlab("Z-score of residuals")+
  ylab("Count")
#another way of showing fdr vs residual
hist_fdr_residuals <- ggplot(pdz3cript_dm, aes(x=fdr_dm*100))+
  geom_histogram(fill="white", colour="black", breaks=seq(0,100, by = 1))+
  theme_bw()+
  xlab("FDR of residuals")+
  ylab("Count (log)")+scale_y_log10()+geom_vline(xintercept = 10, col="red",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")
  #ylab("Count")+geom_vline(xintercept = 10, col="red",linetype="dashed")
#another way using hexbins/showing residuals so as not to plot all 200k datapoints
fdr_resid_pred_plot <- function(df_x){
  dist_plot_nofdr <- ggplot(df_x, aes(x=mean, y=`pdz3cript_resids[pdz3cript_double_mut_pos]`))+
    geom_hex(bins=60)+
    stat_summary_hex(aes(z=fdr_dm, color=after_stat(as.character(value))), bins=60, fun = ~+all(.x<0.1), fill = NA, linewidth=0.25)+
    scale_color_manual(values=c("0" = "transparent", "1" = "#F9CB64"), guide="none")+
    #scale_fill_gradient(low="#5E5E5E", high="lightgrey")+
    scale_fill_gradientn(colours = viridis(256, option = "B"))+
    xlab("Predicted fitness")+
    ylab("Residual")+
    scale_y_continuous(limits=c(min(df_x$`pdz3cript_resids[pdz3cript_double_mut_pos]`), max(df_x$`pdz3cript_resids[pdz3cript_double_mut_pos]`)))+
    theme_bw()+theme(legend.key.width = unit(0.3, "cm"))
  dist_plot_nofdr + geom_hline(yintercept = 0, linetype="dashed", color = "lightgrey", linewidth=0.25)
}
predfit_resid_pdz3cript <- fdr_resid_pred_plot(pdz3cript_dm)
###
pdz3cript_dm$positive_epi <- rep(0, nrow(pdz3cript_dm))
pdz3cript_dm$positive_epi[which(pdz3cript_dm$fitness>dead_var_cutoff & pdz3cript_dm$`pdz3cript_resids[pdz3cript_double_mut_pos]`>0 & pdz3cript_dm$passes_fdr==1)] <- 1
#creating an enrichment test for each pair of positions as compared to the rest of the dataset (with FDR calculated at level of whole dataset)
pdz3cript_dm_split <- split.data.frame(pdz3cript_dm, f=pdz3cript_dm$pos_pairs_realprot_coord)
#for each pair i also want to see if it's enriched in specificity changing mutations as compared to the whole dataset
hit_in_samp_per_pair <- lapply(pdz3cript_dm_split, FUN=function(x){sum(x$positive_epi)})
sample_size_per_pair <- lapply(pdz3cript_dm_split, FUN=function(x){nrow(x)})
hit_in_pop <- sum(do.call(rbind, pdz3cript_dm_split)$positive_epi)
fail_in_pop <- nrow(do.call(rbind, pdz3cript_dm_split))-hit_in_pop
hypergeom_test_per_pair <- phyper(unlist(hit_in_samp_per_pair)-1, hit_in_pop, fail_in_pop, unlist(sample_size_per_pair), lower.tail = F)
df_hypergeom_test <- data.frame(pair_name=names(hypergeom_test_per_pair),
                                #mut_name=unlist(lapply(pdz3cript_dm_split, FUN=function(x){unique(x$mut1_mut2_realprot_coord)})),
                                type_intxn=unlist(lapply(pdz3cript_dm_split, FUN=function(x){unique(x$type_cript)})),
                                pvalue_phyper_test=as.numeric(hypergeom_test_per_pair),
                                pvalue_phyper_test_fdr=p.adjust(as.numeric(hypergeom_test_per_pair), method="fdr", n=length(pdz3cript_dm_split)),
                                hit_in_samp=as.numeric(unlist(hit_in_samp_per_pair)),
                                position_a=unlist(lapply(names(pdz3cript_dm_split), FUN=function(x){strsplit(x, "_")[[1]][1]})),
                                position_b=unlist(lapply(names(pdz3cript_dm_split), FUN=function(x){strsplit(x, "_")[[1]][2]})),
                                aa_pos_a=unlist(lapply(pdz3cript_dm_split, FUN=function(x){str_sub(unique(strsplit(x$mut1_mut2_realprot_coord, "_")[[1]][1]),1,-2)})),
                                aa_pos_b=unlist(lapply(pdz3cript_dm_split, FUN=function(x){str_sub(unique(strsplit(x$mut1_mut2_realprot_coord, "_")[[1]][2]), 1, -2)})))
bf_p_value_cutoff <- (0.05/length(pdz3cript_dm_split))
df_hypergeom_test$aa_pos_a <- factor(df_hypergeom_test$aa_pos_a, levels=paste(as.vector(s2c(wt_seq)), (1+302):(nchar(wt_seq)+302), sep=""))
df_hypergeom_test$aa_pos_b <- factor(df_hypergeom_test$aa_pos_b, levels=c("K-7", "N-6", "Y-5", "K-4", "Q-3", "T-2", "S-1", "V0"))
hypergeom_plot <- ggplot(data=df_hypergeom_test)+
  geom_point(aes(x=aa_pos_a, y=pvalue_phyper_test_fdr, fill=type_intxn), size=2, shape=21)+
  #if you only want to colour here significant stuff
  #geom_point(data = subset(df_hypergeom_test, -log10(pvalue_phyper_test) > -log10(p_value_cutoff)), aes(x=aa_pos_a, y=-log10(pvalue_phyper_test), fill=type_intxn), size=2, shape=21)+
  #if you want to colour everything by type of intxn
  #geom_point(data = subset(df_hypergeom_test, -log10(pvalue_phyper_test) > -log10(p_value_cutoff)), aes(x=aa_pos_a, y=-log10(pvalue_phyper_test), fill=type_intxn), size=2, shape=21)+
  scale_fill_manual(values=custom_scale_discrete_4)+
  xlab("PDZ position")+
  ylab("-log10 pvalue")+
  #for only p-value cutoff passing stuff
  #labs(fill=paste("Pairs w p<", round(p_value_cutoff, 4), sep=""))+
  labs(fill=" ")+
  #scale_y_continuous(trans="log10")+
  #labeling it so that only the top p values are labelled, between pairs
  #geom_text(data = subset(df_hypergeom_test, -log10(pvalue_phyper_test) > -log10(0.0001)), aes(x=aa_pos_a, y=-log10(pvalue_phyper_test), label=pair_name),vjust=0, nudge_y=1, size=2.5)+
  #label only the cript residue bc pdz is already on the plot
  geom_text(data = subset(df_hypergeom_test, 0.01 > df_hypergeom_test$pvalue_phyper_test), aes(x=aa_pos_a, y=pvalue_phyper_test, label=paste(aa_pos_b, sep="")),vjust=1, hjust=1.5,nudge_x=0.5, size=2.5)+
  #geom_text(aes(x=aa_pos_a, y=-log10(pvalue_phyper_test), label=paste(aa_pos_b, sep="")),vjust=0, nudge_x=0.15, size=2.5, log10="y")+
  theme_bw()
hypergeom_plot+  
  theme(axis.text.x = element_text(angle = 90))
  #geom_hline(yintercept=-log10(p_value_cutoff), color="grey", linetype="dashed")
#make a heatmap version of this
df_hypergeom_test$hit_in_samp_thresholded <- df_hypergeom_test$hit_in_samp
df_hypergeom_test$hit_in_samp_thresholded[which(df_hypergeom_test$hit_in_samp>10)] <- 10
heat_hitinsamp_matrix <- ggplot(df_hypergeom_test, aes(x=position_a, y=aa_pos_b))+
  geom_tile(aes(fill=hit_in_samp_thresholded))+
  scale_fill_gradient2(low="white", high="#5555ff", midpoint=0, name="No. of spec. changing muts", 
                       labels=c("0", "5","10+"), breaks=c(0,5,10))+
  #scale_fill_gradient2(low="midnightblue", mid="white", high="orangered", midpoint=0, name="No. of spec. changing muts")+
  xlab("PDZ3 position")+
  ylab("CRIPT position")+
  #stuff that passes fdr<0.1
  geom_tile(data=df_hypergeom_test[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1),],aes(x=position_a, y=aa_pos_b), color="black", fill="transparent", linewidth=0.1, linetype="dashed")+
  geom_tile(data=df_hypergeom_test[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.01),],aes(x=position_a, y=aa_pos_b), color="black", fill="transparent", linewidth=0.2)+
  scale_x_discrete(breaks=seq(303,402,10))+
  #theme_bw()+
  theme(axis.text.x = element_text(size = 8,
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 90),
        axis.text.y = element_text(size=10),
        legend.key.height = unit(0.3, "cm"),
        panel.border = element_rect(fill=NA,colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.direction = "horizontal", legend.position = "bottom")
heat_hitinsamp_matrix
#general histogram to show number of pairs with significant residuals
hist_num_spec_changing_per_pair <- ggplot(df_hypergeom_test, aes(x=hit_in_samp))+
  geom_histogram(fill="white", colour="black", bins=100)+
  theme_bw()+
  xlab("Number of spec. changing mutations per position pair")+
  ylab("Count (log)")
#how many major spec. changers at each pdz residue?
num_spec_changers_cript_res <- reshape2::melt(table(df_hypergeom_test$position_b[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1)]))
barplot_num_spec_changers_cript_res <- ggplot(num_spec_changers_cript_res, aes(x=Var1, y=value))+
  geom_bar(stat="identity", position=position_dodge(), fill="white", color="black")+
  theme_bw()+
  xlab("CRIPT position")+ylab("Number of major spec. encoding residues")+
  scale_x_continuous(breaks=seq(-5,0,1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#also breaking down the number of spec changing per pdz position, cript position etc
df_spec_changing_pdz <- reshape2::melt(table(unlist(lapply(unique(pdz3cript_dm$pdz_mut1_realprot_coord[which(pdz3cript_dm$positive_epi==1)]), get_mut_pos))))
barplot_num_spec_changing_pdz <- ggplot(data=df_spec_changing_pdz, aes(x=Var1, y=value))+
  geom_bar(stat="identity", fill="black", color="black", width=0.5, position = position_dodge(width=0, preserve = "total"))+
  theme_bw()+
  xlab("PDZ position")+ylab("Number of spec. changing mutations")+
  scale_x_continuous(breaks=seq(303,402,10))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.title.y = element_blank())#+ggtitle("Number of spec. changing mutations unweighted")
df_spec_changing_cript <- reshape2::melt(table(-(unlist(lapply(unique(pdz3cript_dm$cript_mut2_realprot_coord[which(pdz3cript_dm$positive_epi==1)]), get_mut_pos)))))
barplot_num_spec_changing_cript <- ggplot(data=df_spec_changing_cript, aes(x=Var1, y=value))+
  geom_bar(stat="identity", fill="black", color="black", width=0.5, position = position_dodge(width=0, preserve = "total"))+
  theme_bw()+
  xlab("CRIPT position")+ylab("Number of spec. changing mutations")+
  scale_x_continuous(breaks=seq(-7,0,1))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#plot pdz and cript together
pdz_cript_number_spec_changing <- plot_grid(barplot_num_spec_changing_cript, barplot_num_spec_changing_pdz, rel_widths = c(2.2,7), align = "hv")
#is the number of spec. changing mutations at each site based on coverage? --for supp (a: no)
#plot the number of reads per pdz position
pdz3cript_dm_split_pdzpos <- split.data.frame(pdz3cript_dm, f=str_sub(pdz3cript_dm$pdz_mut1_realprot_coord, 1,-2))
sum_all_seqs <- sum(unlist(lapply(pdz3cript_dm_split_pdzpos, FUN=function(x){sum(x$mean_count)})))
order_pdz <- unique(str_sub(unique(pdz3cript_dm$pdz_mut1_realprot_coord),1,-2))[order(as.numeric(str_sub(unique(str_sub(unique(pdz3cript_dm$pdz_mut1_realprot_coord),1,-2)),2,-1)), decreasing=F)]
df_reads_pdz_pos <- reshape2::melt(unlist(lapply(pdz3cript_dm_split_pdzpos, FUN=function(x){sum(x$mean_count)}))[order_pdz]/sum_all_seqs)
df_reads_pdz_pos$pos <- rownames(df_reads_pdz_pos)
df_reads_pdz_pos$pos <- factor(df_reads_pdz_pos$pos, levels=order_pdz)
df_reads_pdz_pos$pos_num <- as.numeric(str_sub(df_reads_pdz_pos$pos,2,4))
barplot_reads_pdz_pos <- ggplot(data=df_reads_pdz_pos, aes(x=pos_num, y=value*100))+
  geom_bar(stat="identity", position=position_dodge(), fill="white", color="black", width = 0.5)+
  theme_bw()+
  xlab("PDZ position")+ylab("Percentage of reads / total reads in input")+
  scale_x_continuous(breaks=seq(303,402,20))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_hline(yintercept=1/100*100, col="red", linetype="dashed")+ggtitle("Percentage of reads/total [red dash = expected for totally uniform distribution]")
#does not really change when you weight the number of spec changing mutations in pdz by the % reads
number_mut_spec_changing_per_pdz <- reshape2::melt(table(unlist(lapply(unique(pdz3cript_dm$pdz_mut1_realprot_coord[which(pdz3cript_dm$positive_epi==1)]), get_mut_pos))))
number_mut_spec_changing_per_pdz_weighted <- number_mut_spec_changing_per_pdz
#some positions that don't have any spec. changing muts
missing_pos_pdz_spec_changing <- c(303:402)[-match(number_mut_spec_changing_per_pdz_weighted$Var1, 303:402)]
missing_pos_pdz_spec_changing_num <- rep(0, length(missing_pos_pdz_spec_changing))
missing_pos_pdz_df <- data.frame(Var1=missing_pos_pdz_spec_changing, value=missing_pos_pdz_spec_changing_num)
number_mut_spec_changing_per_pdz_weighted <- rbind(number_mut_spec_changing_per_pdz_weighted, missing_pos_pdz_df)
number_mut_spec_changing_per_pdz_weighted <- number_mut_spec_changing_per_pdz_weighted[order(number_mut_spec_changing_per_pdz_weighted$Var1),]
number_mut_spec_changing_per_pdz_weighted$value <- number_mut_spec_changing_per_pdz_weighted$value*df_reads_pdz_pos$value
barplot_num_spec_changing_pdz_weighted <- ggplot(data=number_mut_spec_changing_per_pdz_weighted, aes(x=Var1, y=value))+
  geom_bar(stat="identity", position=position_dodge(), fill="white", color="black")+
  theme_bw()+
  xlab("PDZ position")+ylab("Number of spec. changing mutations")+
  scale_x_continuous(breaks=seq(303,402,1))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Number of spec. changing mutations weighted")

weighted_vs_not_vs_reads_pdz <- plot_grid(barplot_num_spec_changing_pdz, barplot_num_spec_changing_pdz_weighted, barplot_reads_pdz_pos, nrow=3, ncol=1)

#contact analysis
########################################################################################start here:
# plot vs. distance to ligand 
#source script from Faure et al. 2022 to get minimum sidechain heavy atom distance to ligand
source("./fuzzy_specificity_code/source_functions/doubledeepms__minimum_interchain_distances_from_PDB_v2.R")
#the above script calculates distances between sidechains of heavy atoms per residue in chain A vs. B of pdb structure
#in the case of 5heb (pdb struc that i'm using) -- cript residues are labeled 2-9 (and K-7 is missing)
#thus want to re-name colnames accordingly to match the indexing for cript and add a column for K-7 which is basically identical to N-6 (the next residue over)
colnames(pdb5heb_dists_perres)[which(colnames(pdb5heb_dists_perres) %in% c(3,4,5,6,7,8,9))] <- c("-6", "-5", "-4", "-3", "-2", "-1", "0")
pdb5heb_dists_perres <- cbind(pdb5heb_dists_perres[,1:3], pdb5heb_dists_perres[,4], pdb5heb_dists_perres[,4:10])
colnames(pdb5heb_dists_perres)[4] <- "-7"
df_with_dist_resid <- data.frame(pdz3_cript_mut=paste(pdz3cript_dm$pdz_mut1, pdz3cript_dm$cript_mut2_realprot_coord, sep="|"), 
                                 dist_to_lig_min=pdb5heb_dists_perres$scHAmin_ligand[match(unlist(lapply(pdz3cript_dm$pdz_mut1_realprot_coord, get_mut_pos)), pdb5heb_dists_perres$Pos)],
                                 pdz_pos=unlist(lapply(pdz3cript_dm$pdz_mut1_realprot_coord, get_mut_pos)),
                                 binding_resid=pdz3cript_dm$`pdz3cript_resids[pdz3cript_double_mut_pos]`,
                                 fdr_dm=pdz3cript_dm$fdr_dm, 
                                 observed_fitness=pdz3cript_dm$fitness,
                                 type_cript=pdz3cript_dm$type_cript,
                                 pos_cript=unlist(lapply(strsplit(pdz3cript_dm$pos_pairs_realprot_coord, "_"), FUN=function(x){x[2]})),
                                 pdz_cript_pair=pdz3cript_dm$pos_pairs_realprot_coord)
#add a column for inter-residue distance (looks similar but more specific)
pdb5heb_dists_perres <- as.data.frame(pdb5heb_dists_perres)
matched_rows_dist <- match(df_with_dist_resid$pdz_pos,pdb5heb_dists_perres$Pos)
matched_cols_dist <- match(df_with_dist_resid$pos_cript, colnames(pdb5heb_dists_perres))
dists_matched <- NULL
for(i in 1:nrow(df_with_dist_resid)){
  dists_matched[i] <- (pdb5heb_dists_perres[matched_rows_dist[i], matched_cols_dist[i]])
  }
df_with_dist_resid$dist_to_lig <- dists_matched
#separate this further by cript residue
df_with_dist_resid_split_cript <- split.data.frame(df_with_dist_resid, f=df_with_dist_resid$pos_cript)
#original version with FDR coloured continuously
make_dist_plot_fdr_plot <- function(df_x){
  dist_plot_fdr <- ggplot(df_x, aes(x=dist_to_lig, y=binding_resid, label=pdz3_cript_mut))+
    geom_point(size=1.5, aes(col=fdr_dm*100))+
    scale_color_viridis_c(direction=-1, option="B", begin= 0.1, end = 0.9)+
    xlab("Dist. PDZ3 residue to lig. (A)")+
    ylab("Residual")+
    labs(color="FDR of residual (%)")+
    scale_y_continuous(limits=c(min(df_with_dist_resid$binding_resid), max(df_with_dist_resid$binding_resid)))+
    theme_bw()+theme(legend.key.height = unit(0.3, "cm"))
  dist_plot_fdr + geom_vline(xintercept = 5, linetype="dashed", color = "darkgrey", linewidth=1.5) + ggtitle(unique(df_x$pos_cript))
}
#contour coloured only
make_dist_plot_nofdr_plot <- function(df_x){
  dist_plot_nofdr <- ggplot(df_x, aes(x=dist_to_lig, y=binding_resid))+
    geom_hex(bins=20)+
    stat_summary_hex(aes(z=fdr_dm, color=after_stat(as.character(value))), bins=20, fun = ~+all(.x<0.1), fill = NA)+
    scale_color_manual(values=c("0" = "transparent", "1" = "#F9CB64"), guide="none")+
    scale_fill_gradient(low="#5E5E5E", high="lightgrey")+
    xlab("Dist. PDZ3 residue to lig. (A)")+
    ylab("Residual")+
    scale_y_continuous(limits=c(min(df_with_dist_resid$binding_resid), max(df_with_dist_resid$binding_resid)))+
    theme_bw()+theme(legend.key.height = unit(0.3, "cm"))
  dist_plot_nofdr + geom_vline(xintercept = 5, linetype="dashed", color = "orangered", linewidth=1) + ggtitle(unique(df_x$pos_cript))
}
dist_plot_fdr_cript_residue <- lapply(df_with_dist_resid_split_cript, make_dist_plot_nofdr_plot)
dist_plot_fdr_cript_residue_ordered <- dist_plot_fdr_cript_residue[c("0", "-1", "-2", "-3", "-4", "-5", "-6", "-7")]
grid_dist_cript_res_plots <- ggarrange(plotlist=dist_plot_fdr_cript_residue_ordered, nrow=2, ncol=4, common.legend = T, legend="bottom")

pdz_pos_posepi_by_cript <- lapply(df_with_dist_resid_split_cript, FUN=function(x){x$pdz_pos[which(x$fdr_dm<0.1)]})
summary_pdz_posepi_by_cript <- reshape2::melt(lapply(pdz_pos_posepi_by_cript, table))
cript_cols_matrix_mod <- cript_cols_matrix
names(cript_cols_matrix_mod) <-  -(as.numeric(unlist(lapply(names(cript_cols_matrix), str_sub, -1,-1))))
pdz_pos_posepi_by_cript_plot <- ggplot(summary_pdz_posepi_by_cript, aes(fill=L1, y=value, x=Var1))+
      geom_bar(position="stack", stat="identity")+
      scale_x_continuous(breaks=seq(303,402,5))+
      scale_fill_manual(values=cript_cols_matrix_mod, limits=names(cript_cols_matrix_mod), name= " ")+
      xlab("PDZ3 position")+
      ylab("Number of spec. changing mutations")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      theme_bw()
#separate by pdz position
df_with_dist_resid_split_pdz <- split.data.frame(df_with_dist_resid, f=df_with_dist_resid$pdz_cript_pair)
df_with_dist_resid_split_pdz_absres <- lapply(df_with_dist_resid_split_pdz, FUN=function(x){
  data.frame(mean_pos_binding_res = if(length(x$binding_resid)>10){mean(abs(x$binding_resid[which(x$binding_resid>0)]))}else{NA},
             pos_pdz = unique(x$pdz_pos),
             pos_cript= unique(x$pos_cript),
             dist_to_lig = unique(x$dist_to_lig))
             })
df_with_dist_resid_split_pdz_absres_bound <- do.call(rbind,df_with_dist_resid_split_pdz_absres)
df_with_dist_resid_split_pdz_absres_bound_nona <- na.omit(df_with_dist_resid_split_pdz_absres_bound)
split_df_with_dist_resid_pdz_absres_nona <- split.data.frame(df_with_dist_resid_split_pdz_absres_bound_nona, f=df_with_dist_resid_split_pdz_absres_bound_nona$pos_cript)
#do the distance calculation residue<->residue
make_dist_plot_mean_resid_per_pos <- function(df_to_plot){
  dist_plot_fdr_intxn_pair <- ggplot(df_to_plot, aes(x=dist_to_lig, y=mean_pos_binding_res))+
    geom_point(size=1.5)+
    #geom_errorbar(aes(ymin=`N_term_a3_min`, ymax=`N_term_a3_max`, xmin=`N_term_min`, xmax=`N_term_max`))+
    xlab("distance to ligand")+
    ylab("mean abs. binding fitness pos. residual")+
    scale_x_continuous(limits=c(0,25))+
    scale_y_continuous(limits=c(0,0.3))+
    theme_bw()
  dist_plot_fdr_intxn_pair + geom_vline(xintercept = 5, linetype="dashed", color = "darkgrey", linewidth=1.5)+
    ggtitle(unique(df_to_plot$pos_cript))
}
split_df_with_dist_resid_pdz_absres_nona_ordered <- split_df_with_dist_resid_pdz_absres_nona[c("0", "-1", "-2", "-3", "-4", "-5", "-6", "-7")]
plots_abs_resid_per_pair <- lapply(split_df_with_dist_resid_pdz_absres_nona_ordered, make_dist_plot_mean_resid_per_pos)
plots_mean_pos_resid <- plot_grid(plotlist=plots_abs_resid_per_pair, nrow=2, ncol=4)

#read in the contacts from getcontacts
get_contacts_results <- read.delim("./fuzzy_specificity_data/getcontacts_analysis/5heb_contacts_apr324_reformatcol_headers.tsv")
res1_chain <- unlist(lapply(get_contacts_results$atom_1, FUN=function(x){strsplit(x, ":")[[1]][1]}))
res1_index <- unlist(lapply(get_contacts_results$atom_1, FUN=function(x){strsplit(x, ":")[[1]][3]}))
res2_chain <- unlist(lapply(get_contacts_results$atom_2, FUN=function(x){strsplit(x, ":")[[1]][1]}))
res2_index <- unlist(lapply(get_contacts_results$atom_2, FUN=function(x){strsplit(x, ":")[[1]][3]}))
res1_res2_chain <- paste(res1_chain, res2_chain, sep="_")
#fix the indices for cript residues
res1_index_realcoord <- res1_index
res1_index_realcoord[which(res1_chain=="B")] <- as.numeric(res1_index[which(res1_chain=="B")])-9
res2_index_realcoord <- res2_index
res2_index_realcoord[which(res2_chain=="B")] <- as.numeric(res2_index[which(res2_chain=="B")])-9
res1res2_pairs <- paste(res1_index_realcoord, res2_index_realcoord, sep="_")
#make a df
get_contacts_results_df <- cbind(get_contacts_results, res1_chain, res1_index, res2_chain, res2_index, res1_res2_chain,
                                 res1_index_realcoord, res2_index_realcoord, res1res2_pairs)

sigpairs_direct_contact <- get_contacts_results_df[na.omit(match(df_hypergeom_test$pair_name[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1)], get_contacts_results_df$res1res2_pairs)),]
sigpairs_pdz_contact_a <- get_contacts_results_df[na.omit(match(df_hypergeom_test$position_a[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1)], get_contacts_results_df$res1_index)),]
sigpairs_pdz_contact_a_nodup <- sigpairs_pdz_contact_a[-grep("\\.", rownames(sigpairs_pdz_contact_a)),]
sigpairs_pdz_contact_b <- get_contacts_results_df[na.omit(match(df_hypergeom_test$position_a[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1)], get_contacts_results_df$res2_index)),]
sigpairs_pdz_contact_b_nodup <- sigpairs_pdz_contact_b[-grep("\\.", rownames(sigpairs_pdz_contact_b)),]
sigpairs_pdz_contact_ab_nodup <- rbind(sigpairs_pdz_contact_a_nodup, sigpairs_pdz_contact_b_nodup)


#for plotting the structures
get_contacts_results_df_melted_intxn_type <- get_contacts_results_df
get_contacts_results_df_melted_intxn_type <- get_contacts_results_df_melted_intxn_type[,-c(3,4)]
aggregated_df_contact_map <- aggregate(get_contacts_results_df_melted_intxn_type$interaction_type, list(get_contacts_results_df_melted_intxn_type$res1res2_pairs, get_contacts_results_df_melted_intxn_type$res1_res2_chain), paste, collapse="_")
colnames(aggregated_df_contact_map) <- c("interaction_pairs", "res1_res2_chain", "interaction_type")
aggregated_df_contact_map$interaction_type <- unlist(lapply(aggregated_df_contact_map$interaction_type, FUN=function(x){paste(unique(unlist(strsplit(x, "_"))), collapse="_")}))

#want to represent in a new df
df_contacts_heatmap <- df_hypergeom_test
df_contacts_heatmap$contact_type <- aggregated_df_contact_map$interaction_type[match(df_hypergeom_test$pair_name, aggregated_df_contact_map$interaction_pairs)]
df_contacts_heatmap$contact_type[which(is.na(df_contacts_heatmap$contact_type)==T)] <- "none"
df_contacts_heatmap$contact_type <- factor(df_contacts_heatmap$contact_type, levels=unique(df_contacts_heatmap$contact_type))
#for supplementary -- types of contacts via getcontacts for each pair of positions
df_contacts_matrix <- ggplot(df_contacts_heatmap, aes(x=aa_pos_a, y=aa_pos_b))+
    geom_tile(aes(fill=contact_type))+
    scale_fill_manual(values=c("white", RColorBrewer::brewer.pal(10, "Set3")), name="contact type", labels=unique(as.vector(df_contacts_heatmap$contact_type)))+
    xlab("PDZ3 position")+
    ylab("CRIPT position")+
    #stuff that passes fdr<0.1
    geom_tile(data=df_contacts_heatmap[which(df_contacts_heatmap$pvalue_phyper_test_fdr<0.1),],aes(x=aa_pos_a, y=aa_pos_b), color="black", fill="transparent", linewidth=0.1, linetype="dashed")+
    geom_tile(data=df_contacts_heatmap[which(df_contacts_heatmap$pvalue_phyper_test_fdr<0.01),],aes(x=aa_pos_a, y=aa_pos_b), color="black", fill="transparent", linewidth=0.2)+
    theme(axis.text.x = element_text(size = 8,
                                     vjust = 0.5,
                                     hjust = 0.5,
                                     angle = 90),
          axis.text.y = element_text(size=10),
          legend.key.width = unit(0.3, "cm"),
          panel.border = element_rect(fill=NA,colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
  theme(legend.key.height = unit(0.3, "cm"), legend.direction = "horizontal", legend.position = "bottom")
  df_contacts_matrix

df_hypergeom_sig_pairs <- df_hypergeom_test[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1),]
df_hypergeom_sig_pairs_split  <- split.data.frame(df_hypergeom_sig_pairs, f=df_hypergeom_sig_pairs$position_b)
#want to see how many interactions there are between significant pdz stuff:
all_poss_sig_pos_interactions_pdz <- unique(apply(expand.grid(as.vector(df_hypergeom_sig_pairs$position_a), as.vector(df_hypergeom_sig_pairs$position_a)), 1, paste, collapse="_"))
#make the same contact map but for pdz-pdz
aggregated_df_contact_map$aa_pos_a <- unlist(lapply(aggregated_df_contact_map$interaction_pairs, FUN=function(x){strsplit(x, "_")[[1]][1]}))
aggregated_df_contact_map$aa_pos_b <- unlist(lapply(aggregated_df_contact_map$interaction_pairs, FUN=function(x){strsplit(x, "_")[[1]][2]}))
aggregated_df_contact_map_only_pdz <- aggregated_df_contact_map[which(aggregated_df_contact_map$res1_res2_chain=="A_A"),]
aggregated_df_contact_map_only_pdz$interacting_with_other_interactors <- rep(0, nrow(aggregated_df_contact_map_only_pdz))
aggregated_df_contact_map_only_pdz$interacting_with_other_interactors[which(aggregated_df_contact_map_only_pdz$interaction_pairs %in% all_poss_sig_pos_interactions_pdz)] <- 1
df_contacts_matrix_pdz <- ggplot(aggregated_df_contact_map_only_pdz, aes(x=aa_pos_a, y=aa_pos_b))+
  geom_tile(aes(fill=interaction_type))+
  scale_fill_manual(values=viridis(24), name="contact type", labels=unique(as.vector(aggregated_df_contact_map_only_pdz$interaction_type)))+
  xlab("PDZ3 position")+
  ylab("PDZ3 position")+
  geom_tile(data=aggregated_df_contact_map[which(aggregated_df_contact_map_only_pdz$interacting_with_other_interactors==1),],aes(x=aa_pos_a, y=aa_pos_b), color="black", fill="transparent", linewidth=0.5)+
  theme(axis.text.x = element_text(size = 8,
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 90),
        axis.text.y = element_text(size=10),
        legend.key.width = unit(0.3, "cm"),
        panel.border = element_rect(fill=NA,colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.key.height = unit(0.3, "cm"), legend.direction = "horizontal", legend.position = "bottom")

#looking at each significant pair/network of contacts in pdz and otherwise
make_pdz_contact_matrix_plot <- function(df_x){
  df_x_positions_of_interest <- c(which(aggregated_df_contact_map_only_pdz$aa_pos_a %in% df_x$position_a),
                                  which(aggregated_df_contact_map_only_pdz$aa_pos_a %in% df_x$position_b))
  df_for_geomtile <- aggregated_df_contact_map_only_pdz[df_x_positions_of_interest,]
  df_contacts_matrix_pdz <- ggplot(df_for_geomtile, aes(x=(aa_pos_b), y=(aa_pos_a)))+
    geom_tile(aes(fill=interaction_type))+
    scale_fill_manual(values=viridis(7), name="contact type", labels=unique(as.vector(df_for_geomtile$interaction_type)))+
    xlab("PDZ3 position")+
    ylab("PDZ3 position")+
    theme_bw()+
    theme(axis.text.x = element_text(size = 8,                                 
                                     vjust = 0.5,
                                     hjust = 0.5,
                                     angle = 90),
          axis.text.y = element_text(size=8),
          legend.key.width = unit(0.3, "cm"),
          panel.border = element_rect(fill=NA,colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme(legend.key.width=unit(0.3,"cm"))
    #theme(legend.key.height = unit(0.3, "cm"), legend.direction = "horizontal", legend.position = "bottom")
  df_contacts_matrix_pdz
}
pdz_pos_matrix_contacts <- make_pdz_contact_matrix_plot(df_hypergeom_sig_pairs)


#distance between pairs of all the significant things
dist_sig_spec_changing <- (lapply(df_with_dist_resid_split_pdz, FUN=function(x){unique(x$dist_to_lig)}))[which(names(lapply(df_with_dist_resid_split_pdz, FUN=function(x){unique(x$dist_to_lig)})) %in% df_hypergeom_test$pair_name[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1)])]

#make a scatterplot of the distance vs. # of spec. changing mutations
df_dist_fdr_spec_change <- df_hypergeom_test
if(identical(names(lapply(df_with_dist_resid_split_pdz, FUN=function(x){unique(x$dist_to_lig)})), df_hypergeom_test$pair_name)==T){
  df_dist_fdr_spec_change$dist_to_lig <- unlist(lapply(df_with_dist_resid_split_pdz, FUN=function(x){unique(x$dist_to_lig)}))
}else{print("check id of pairs")}

dist_vs_spec_changing <- ggplot(df_dist_fdr_spec_change, aes(x=dist_to_lig, y=hit_in_samp, fill=pvalue_phyper_test_fdr*100))+
  geom_point(size=1.5, shape=21)+
  scale_fill_viridis_c(direction=-1, option="B", begin= 0.1, end = 0.9)+
  theme(element_text(size=8, family="serif"))+
  xlab("Distance between PDZ3-CRIPT pair")+
  ylab("Number of spec. changing mutations")+
  labs(fill="FDR of enrichment")+
  theme_bw()+theme(legend.key.height = unit(0.3, "cm"), legend.direction = "horizontal", legend.position = "bottom")

#look at the significant spec. changers and categorize:
cast_significant <- reshape2::dcast(df_hypergeom_test, position_a~position_b, value.var="pvalue_phyper_test_fdr")
list_of_sig_positions_pdz <- NULL
for (i in colnames(cast_significant)[2:ncol(cast_significant)]){
  list_of_sig_positions_pdz[[i]] <- (cast_significant$position_a[which(cast_significant[i]<0.1)])
  }

#plot on structure
cast_hit_in_samp <- reshape2::dcast(df_hypergeom_test, position_a~position_b, value.var="hit_in_samp")
#set the NAs to zero for chimerax
cast_hit_in_samp[which(is.na(cast_hit_in_samp), arr.ind=T)] <- 0
#trying to make some vectors to colour structure -- need to change example set attribute file from chimerax
eg_setatt_file <-read.delim("./fuzzy_specificity_data/eg_files/percentExposed.txt")
home_dir <- getwd()
setwd("./fuzzy_specificity_data/figures/figure5_pdz_specificity/")
setatt_df <- NULL
for (i in 2:ncol(cast_hit_in_samp)){
  setatt_df[[i]] <- data.frame(empty_col=rep("", nrow(cast_hit_in_samp)), pos=paste(":", cast_hit_in_samp$position_a, sep=""), att=cast_hit_in_samp[,i])
  write.table(setatt_df[[i]], paste("hitinsamp_", colnames(cast_hit_in_samp)[i], ".defattr", sep=""), quote=F, sep="\t",row.names=F, col.names=F, eol = "\r")
  content_to_change <- (readLines(paste("hitinsamp_", colnames(cast_hit_in_samp)[i], ".defattr", sep="")))
  new_content <- c("attribute: x","match mode: 1-to-1", "recipient: residues", content_to_change)
  writeLines(new_content, paste("hitinsamp_", colnames(cast_hit_in_samp)[i], ".defattr", sep=""), sep="\n")
}
setwd(home_dir)
#to select residues with 5+ spec. changing mutations in PDZ (to then show sidechains in chimerax)
commands_chimera_sidechains <- NULL
for (i in 2:ncol(cast_hit_in_samp)){
  commands_chimera_sidechains[[i]] <- paste("sel #1/A:",paste(cast_hit_in_samp$position_a[which(cast_hit_in_samp[,i]>=5)], collapse = ","), sep=" ")
  }
#function to put into chimerax after loading the basic structure (custom-coloured already and outlining the relevant cript res.) 
#color byattribute x palette white:#5555ff range 0,10 key true; key white:0 #5555ff:10+
######################################################################################################################################################################
#19x19 maps of the positions that have >>>spec. changing mutations
#for each double mutant pair, want to plot the 19x19 heatmap
make_matrix_for_clustering <- function(df_x, value_to_cluster, filtering){
  potential_matrix <- reshape2::dcast(df_x, pdz_mut1_realprot_coord~cript_mut2_realprot_coord, value.var=as.character(value_to_cluster))
  new_matrix <- potential_matrix[,2:ncol(potential_matrix)]
  rownames(new_matrix) <- potential_matrix[,1]
  #----heatmapping in r
  #filtering matrix pre-clustering using custom funcs (see sourced basic funcs file)
  if(filtering==TRUE){
    filter_matrix_rows_na <- function(x){
      lengths_of_nas <- apply(x, 1, FUN=function(x){length(which(is.na(x)))/length(x)*100})
      if(length(which(lengths_of_nas >= 50))>0){
        dm_matrix_filt <- x[-which(lengths_of_nas >= 50),]
        return(dm_matrix_filt)
      }else{return(x)}
    }
    filter_matrix_cols_na <- function(x){
      lengths_of_nas <- apply(x, 2, FUN=function(x){length(which(is.na(x)))/length(x)*100})
      if(length(which(lengths_of_nas >= 50))>0){
        dm_matrix_filt <- x[,-which(lengths_of_nas >= 50)]
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
    
    filter_matrices_func <- function(x){
      y <- filter_matrix_cols_na(x)
      y <- filter_matrix_rows_na(y)
      y <- filter_matrix_cols_zero(y)
      y <- filter_matrix_rows_zero(y)
      return(y)
    }
    filtered_matrix <- filter_matrices_func(new_matrix)
    filtered_matrix <- as.matrix(filtered_matrix)
    return(filtered_matrix)
  }else{
    return(as.matrix(new_matrix))
  }
}

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
get_breaks_matrices <- function(x){
  #  myBreaks <- c(seq(min(x, na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), ###USE THIS for pheatmap -- it doesn't care about break length being 1+colour length
  #trying something else
  myBreaks <- c(seq(min(x, na.rm=T), 0, length.out=ceiling(paletteLength/2)),
                seq(max(x, na.rm=T)/paletteLength, max(x, na.rm=T), length.out=floor(paletteLength/2)))
  return(myBreaks)
}

#add an annotation layer in here where i can colour stuff by amino acid type
zappo_cols_matrix <- zappo_col_scheme$col
names(zappo_cols_matrix) <- zappo_col_scheme$letter
make_mat_to_clust_w_fdr_layer <- function(x){
  mat_to_cluster <- make_matrix_for_clustering(x, "pdz3cript_resids[pdz3cript_double_mut_pos]", filtering=T)
  mat_fdr_values <- make_matrix_for_clustering(x, "fdr_dm", filtering=F)
  paletteLength <- 50
  myColor <- colorRampPalette(c("midnightblue", "white", "orangered"))(paletteLength)
  #diff color palette
  myColor_v2 <- colorRamp2(seq(-0.6, 0.6, length=3), c("midnightblue", "white", "orangered"))
  heatmap_clust <- Heatmap(mat_to_cluster)
  heatmap_clust_ordered <- heatmap_clust@matrix[row_order(heatmap_clust), column_order(heatmap_clust)]
  fdr_matrix_ordered <- mat_fdr_values[rownames(heatmap_clust_ordered), colnames(heatmap_clust_ordered)]
  mut_physicochem_col <- unlist(lapply(colnames(mat_to_cluster), str_sub, -1, -1))
  mut_physicochem_row <- unlist(lapply(rownames(mat_to_cluster), str_sub, -1,-1))
  col_annot_mat <- HeatmapAnnotation(cn=anno_simple(mut_physicochem_col, col=zappo_cols_matrix),
                                     annotation_name_gp= gpar(fontsize = 10),
                                     annotation_name_side = "left",
                                     show_annotation_name = F)
  row_annot_mat <- rowAnnotation(cn=anno_simple(mut_physicochem_row, col=zappo_cols_matrix),
                                     annotation_name_gp= gpar(fontsize = 10),
                                     show_annotation_name = F)
  Heatmap(mat_to_cluster, col=myColor_v2, name = "binding residual", column_names_gp = grid::gpar(fontsize=8), row_names_gp = grid::gpar(fontsize=8), row_names_side = "left",
          show_column_dend = F, show_row_dend = F,
          bottom_annotation = col_annot_mat,
          left_annotation = row_annot_mat,
          width = ncol(mat_to_cluster)*unit(3, "mm"), 
          height = nrow(mat_to_cluster)*unit(3, "mm"),
          layer_fun = function(j, i, x, y, width, height, fill) {
            #v= pindex(fdr_matrix_ordered, i,j)
            v=fdr_matrix_ordered
            l= v > .05 & v < .10   
            l2= v > .01 & v < .05
            l3= v > 0 & v < .01
            grid.text("*", x[l], y[l], gp = gpar(fontsize = 8))
            grid.text("**", x[l2], y[l2], gp = gpar(fontsize = 8))
            grid.text("***", x[l3], y[l3], gp = gpar(fontsize = 8))
          })
}
heatmaps_sig_pos_w_fdr <- lapply(pdz3cript_dm_split[df_hypergeom_test$pair_name[order(df_hypergeom_test$hit_in_samp, decreasing=T)]][1:20], make_mat_to_clust_w_fdr_layer)
grobbed_hm_cols <- lapply(heatmaps_sig_pos_w_fdr, ggplotify::as.grob)
#all the specificity-changing pairs that pass fdr 10%
grid_hm_top_spec_changing <- cowplot::plot_grid(plotlist=grobbed_hm_cols, nrow=5, ncol=4)
##########################################################################################################################################################
#plot how pleiotropic pdz positions are -- eg how often is one pdz position interacting with more than one cript residue
no_spec_change_per_sig_pair_pdz <- table(unlist(lapply(df_hypergeom_test$pair_name[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1)], FUN=function(x){str_sub(x, 1,3)})))
no_spec_change_per_sig_pair_pdz_melted <- reshape2::melt(no_spec_change_per_sig_pair_pdz)
pleiotropicity_pdz_barplot <- ggplot(no_spec_change_per_sig_pair_pdz_melted, aes(x=as.factor(value)))+
  geom_bar(stat="count", position=position_dodge(), fill="white", color="black")+
  xlab("Number of CRIPT residues that PDZ residue interacts with")+
  ylab("Count (out of all sig. specificity-changing PDZ-CRIPT pairs)")+
  theme_bw()
#max number of cript residues that each pdz position interacts with is 2 -- and there are only 3 (for FDR 10%), 2 (for FDR 1%)
pairs_pleiotropic_pdz <- names(which(no_spec_change_per_sig_pair_pdz>1))
pairs_to_compare_pleiotropic_pdz <- df_hypergeom_test$pair_name[which(df_hypergeom_test$position_a %in% pairs_pleiotropic_pdz & df_hypergeom_test$pvalue_phyper_test_fdr<0.1)]
############################################################################################################################################################
#also wanted to make the pdz matrices (i.e. one matrix per pdz mutation)
pdz3cript_dm_split_pdzpos <- split.data.frame(pdz3cript_dm, f=str_sub(pdz3cript_dm$pdz_mut1_realprot_coord, 1,-2))
#look at the pdz mutations that are implicated in specificity changes 
sig_change_implicated_mut <- which(unlist(lapply(names(pdz3cript_dm_split_pdzpos), FUN=function(x){str_sub(x,2,-1)})) %in% names(no_spec_change_per_sig_pair_pdz))
heatmaps_pdzmut <- lapply(pdz3cript_dm_split_pdzpos[sig_change_implicated_mut], make_mat_to_clust_w_fdr_layer)
grobbed_hm_cols_pdzmut <- lapply(heatmaps_pdzmut, ggplotify::as.grob)
grid_hm_top_spec_changing_pdzmut_1 <- cowplot::plot_grid(plotlist=grobbed_hm_cols_pdzmut[1:8], nrow=8)
grid_hm_top_spec_changing_pdzmut_2 <- cowplot::plot_grid(plotlist=grobbed_hm_cols_pdzmut[9:17], nrow=9)
####write plots to file
#for figure 3, add annotations to the double mutant heatmap
ordering_aa_muts_cript <- NULL 
for (i in 1:length(order_cript)){
  ordering_aa_muts_cript[[i]] <- paste(order_cript[i], new_order_aa, sep="")
}
ordering_aa_muts_cript <- unlist(ordering_aa_muts_cript)
#ordering of pdz positions
order_pdz <- unique(str_sub(unique(pdz3cript_dm$pdz_mut1_realprot_coord),1,-2))[order(as.numeric(str_sub(unique(str_sub(unique(pdz3cript_dm$pdz_mut1_realprot_coord),1,-2)),2,-1)), decreasing=F)]
ordering_aa_muts_pdz <- NULL
for (i in 1:length(order_pdz)){
  ordering_aa_muts_pdz[[i]] <- paste(order_pdz[i], new_order_aa, sep="")
}
ordering_aa_muts_pdz <- unlist(ordering_aa_muts_pdz)
cript_cols_matrix <- viridis(length(order_cript), option = "B")
names(cript_cols_matrix) <- order_cript
pdz_cols_matrix <- viridis(length(order_pdz), option="B")
names(pdz_cols_matrix) <- order_pdz
make_mat_to_clust_mutcolor_fitness <- function(x){
  mat_to_cluster <- make_matrix_for_clustering(x, "fitness", filtering=F)
  #reorder the matrix since it won't be clustered, just want to show annotations
  #order by cript residue/pdz residue and by physicochemical AA properties
  mat_to_cluster <- mat_to_cluster[,index_x_to_y(colnames(mat_to_cluster), ordering_aa_muts_cript)]
  mat_to_cluster <- mat_to_cluster[index_x_to_y(rownames(mat_to_cluster), ordering_aa_muts_pdz),]
  #diff color palette
  myColor_v2 <- colorRamp2(c(-1.3,0,1.3), c("midnightblue", "white", "orangered"))
  mut_physicochem_col <- unlist(lapply(colnames(mat_to_cluster), str_sub, -1, -1))
  mut_physicochem_row <- unlist(lapply(rownames(mat_to_cluster), str_sub, -1,-1))
  mut_position_col <- unlist(lapply(colnames(mat_to_cluster), str_sub,1,-2))
  mut_position_row <- unlist(lapply(rownames(mat_to_cluster), str_sub,1,-2))
  col_annot_mat <- HeatmapAnnotation(cn=anno_simple(mut_physicochem_col, col=zappo_cols_matrix),
                                     cn2=anno_simple(mut_position_col, col=cript_cols_matrix),
                                     annotation_name_gp= gpar(fontsize = 10),
                                     annotation_name_side = "left",
                                     show_annotation_name = F,
                                     simple_anno_size = unit(0.6,"cm"))
  row_annot_mat <- rowAnnotation(cn=anno_simple(mut_physicochem_row, col=zappo_cols_matrix),
                                 cn2=anno_simple(mut_position_row, col=pdz_cols_matrix),
                                 annotation_name_gp= gpar(fontsize = 10),
                                 show_annotation_name = F, simple_anno_size=unit(0.6, "cm"))
  Heatmap(mat_to_cluster, col=myColor_v2, name = "binding fitness", 
          show_column_names = F, show_row_names = F,
          show_column_dend = F, show_row_dend = F,
          bottom_annotation = col_annot_mat,
          left_annotation = row_annot_mat,
          cluster_rows = F,
          cluster_columns = F, na_col = "#EBEBEB")
}
#huge plot of all double mutants (heatmap of binding fitness)
#pheat_dm_annotated <- make_mat_to_clust_mutcolor_fitness(pdz3cript_dm)
##############################################################################################################################
#making some additional plots where i compare the binding of a wt pdz to pdz with mutation to all the cript residues
#compute double mutant positions/ids
sm_diffs_pdzwt <- lapply(pdz_dimsum_wtpdz_consolidated$aa_seq, FUN=list.string.diff, wt_seq)
pdz3cript_sm_df_mut1 <- NULL
for (i in 1:length(sm_diffs_pdzwt)){
  pdz3cript_sm_df_mut1[i] <- paste(sm_diffs_pdzwt[[i]][1,]$poly.seq.b, sm_diffs_pdzwt[[i]][1,]$position, sm_diffs_pdzwt[[i]][1,]$poly.seq.a, sep="")
}
sm_diffs_df <- data.frame(pdz_dimsum_wtpdz_consolidated, cript_mut=pdz3cript_sm_df_mut1, cript_mut_real=unlist(lapply(pdz3cript_sm_df_mut1,change_mut_numbers, -108)), 
                          observed_fitness=pdz_dimsum_wtpdz_consolidated$fitness)
sm_diffs_df <- sm_diffs_df[-which(pdz_dimsum_wtpdz_consolidated$Nham_aa==0),]
#compare, e.g. for each position pair, what the binding fitness of that variant looks like compared to that variant with wt pdz:
df_only_sig_pairs <- df_hypergeom_test[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1),]
df_only_sig_pairs_with_fitness <- pdz3cript_dm_split[df_only_sig_pairs$pair_name]
#for each pair in the list, add in the binding fitness for same cript mut in wt pdz background
df_only_sig_pairs_with_fitness_plus_wtpdz_fitness <- lapply(df_only_sig_pairs_with_fitness, FUN=function(x){
  x$sm_fitness_wt_pdz <- sm_diffs_df$fitness[match(x$cript_mut2_realprot_coord, sm_diffs_df$cript_mut_real)]
  return(x)
})
split_by_mut_sig_pairs_plus_wtpdz_fitness <- lapply(df_only_sig_pairs_with_fitness_plus_wtpdz_fitness, FUN=function(x){
  split.data.frame(x, f=x$pdz_mut1_realprot_coord)
})

#scatterplot for all wt vs mutant fitness per double mutant:
mut_sig_pairs_scatterplots <- lapply(split_by_mut_sig_pairs_plus_wtpdz_fitness, FUN=function(x){
lapply(x, FUN=function(y){make_scatterplot_errorbars(y, y$sm_fitness_wt_pdz, y$fitness, lbl = y$cript_mut2_realprot_coord, error_y=0, yname="Mutant PDZ fitness", xname="Wildtype PDZ fitness")+ggtitle(unique(y$pdz_mut1_realprot_coord))})
  })

#keep only the pdz mutations that have a positive epi signal and split per pdz mutation
#also add wt binding info
pdz3cript_dm$wt_pdz_fitness <- sm_diffs_df$fitness[match(pdz3cript_dm$cript_mut2_realprot_coord, sm_diffs_df$cript_mut_real)]
pdz3cript_dm$wt_pdz_mut_pdz_diff <- pdz3cript_dm$fitness - pdz3cript_dm$wt_pdz_fitness
pdz3cript_dm$cript_pos_realcoord <- -unlist(lapply(pdz3cript_dm$cript_mut2_realprot_coord, get_mut_pos))
#filter for specificity-changing mutations on the cript<->pdz row level (might need to change this if visualizing differently)
pdz3cript_dm_split_criptmut <- split.data.frame(pdz3cript_dm, f=pdz3cript_dm$cript_pos_realcoord)
pdz3cript_dm_split_criptmut_posepi <- lapply(pdz3cript_dm_split_criptmut, FUN=function(x){
  posi_epi_muts <- unique(x$pdz_mut1_realprot_coord[which(x$positive_epi==1)])
  x[which(x$pdz_mut1_realprot_coord %in% posi_epi_muts),]
  })
#alternative to filter spec-changing at top level
pdz3cript_dm_split_criptmut_only_pos_spec <- split.data.frame(pdz3cript_dm[which(pdz3cript_dm$positive_epi==1),], f=pdz3cript_dm[which(pdz3cript_dm$positive_epi==1),]$cript_pos_realcoord)
#plot as heatmap, one for every cript position
ordering_aa_muts_cript <- NULL 
for (i in 1:length(order_cript)){
  ordering_aa_muts_cript[[i]] <- paste(order_cript[i], new_order_aa, sep="")
}
ordering_aa_muts_cript <- unlist(ordering_aa_muts_cript)
#ordering of pdz positions
order_pdz <- unique(str_sub(unique(pdz3cript_dm$pdz_mut1_realprot_coord),1,-2))[order(as.numeric(str_sub(unique(str_sub(unique(pdz3cript_dm$pdz_mut1_realprot_coord),1,-2)),2,-1)), decreasing=F)]
ordering_aa_muts_pdz <- NULL
for (i in 1:length(order_pdz)){
  ordering_aa_muts_pdz[[i]] <- paste(order_pdz[i], new_order_aa, sep="")
}
ordering_aa_muts_pdz <- unlist(ordering_aa_muts_pdz)
cript_cols_matrix <- viridis(length(order_cript), option = "B")
names(cript_cols_matrix) <- order_cript
pdz_cols_matrix <- viridis(length(order_pdz), option="B")
names(pdz_cols_matrix) <- order_pdz
make_mat_to_clust_posepi_fitness <- function(x){
  mat_to_cluster <- make_matrix_for_clustering(x, "fitness", filtering=F)
  #reorder the matrix since it won't be clustered, just want to show annotations
  #order by cript residue/pdz residue and by physicochemical AA properties
  mat_to_cluster <- mat_to_cluster[,index_x_to_y(colnames(mat_to_cluster), ordering_aa_muts_cript)]
  mat_to_cluster <- mat_to_cluster[index_x_to_y(rownames(mat_to_cluster), ordering_aa_muts_pdz),]
  #diff color palette
  myColor_v2 <- colorRamp2(c(-0.6,0.2), c("#EBEBEB", "orangered"))
  mut_physicochem_col <- unlist(lapply(colnames(mat_to_cluster), str_sub, -1, -1))
  mut_physicochem_row <- unlist(lapply(rownames(mat_to_cluster), str_sub, -1,-1))
  mut_position_col <- unlist(lapply(colnames(mat_to_cluster), str_sub,1,-2))
  mut_position_row <- unlist(lapply(rownames(mat_to_cluster), str_sub,1,-2))
  col_annot_mat <- HeatmapAnnotation(cn=anno_simple(mut_physicochem_col, col=zappo_cols_matrix),
                                     #cn2=anno_simple(mut_position_col, col=cript_cols_matrix),
                                     annotation_name_gp= gpar(fontsize = 10),
                                     annotation_name_side = "left",
                                     show_annotation_name = F,
                                     simple_anno_size = unit(0.6,"cm"))
  row_annot_mat <- rowAnnotation(cn=anno_simple(mut_physicochem_row, col=zappo_cols_matrix),
                                 cn2=anno_simple(mut_position_row, col=pdz_cols_matrix),
                                 annotation_name_gp= gpar(fontsize = 10),
                                 show_annotation_name = F, simple_anno_size=unit(0.6, "cm"))
  hm_to_draw <- Heatmap(mat_to_cluster, col=myColor_v2, name = "binding fitness", 
          show_column_names = T, show_row_names = T,
          show_column_dend = F, show_row_dend = F,
          bottom_annotation = col_annot_mat,
          left_annotation = row_annot_mat,
          cluster_rows = F, row_names_gp = gpar(fontsize=7), 
          cluster_columns = F, na_col = "white")
  draw(hm_to_draw)
}
#for general info 
#pheat_fitness_posepi_only <- lapply(pdz3cript_dm_split_criptmut_only_pos_spec, make_mat_to_clust_posepi_fitness)
#################################################################################################################################################################
#checking all the heatmaps without filtering; only take out the pdz mutations that have no spec. changing mut across any  cript positions
#split pdz3cript dm by pdz mutation, then filter, then re-partition
split_pdz3criptdm_pos_epi_mut <- split.data.frame(pdz3cript_dm, f=pdz3cript_dm$pdz_mut1_realprot_coord)
number_of_spec_changing_muts_per_pdz_mut <- lapply(split_pdz3criptdm_pos_epi_mut, FUN=function(x){sum(x$positive_epi)})
#there are 340 non-zero mutations
split_pdz3criptdm_pos_epi_mut_filtered_one_spec_changing_per_row <- split_pdz3criptdm_pos_epi_mut[-which(number_of_spec_changing_muts_per_pdz_mut==0)]
#re-bind and re-split
pdz3cript_dm_pos_epi_atleast_one_row <- do.call(rbind, split_pdz3criptdm_pos_epi_mut_filtered_one_spec_changing_per_row)
pdz3cript_dm_pos_epi_atleast_one_row_cript_split <- split.data.frame(pdz3cript_dm_pos_epi_atleast_one_row, f=pdz3cript_dm_pos_epi_atleast_one_row$cript_pos_realcoord)
#split further to keep only the rows that have one pos epi at least for that specific mutation with cript
filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split <- lapply(pdz3cript_dm_pos_epi_atleast_one_row_cript_split, FUN=function(x){
  split_df <- split.data.frame(x, f=x$pdz_mut1_realprot_coord)
  split_df_filtered <- split_df[which(lapply(split_df, FUN=function(y){sum(y$positive_epi)})>0)]
  do.call(rbind, split_df_filtered)
})
####doing it differently: where i only keep rows from spec-changing enriched pairs to see how it looks
filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split <- lapply(filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split, FUN=function(x){
  x <- x[which(x$pos_pairs_realprot_coord %in% df_hypergeom_test$pair_name[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1)]),]
  return(x)
})
filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split <- filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split[which(unlist(lapply(filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split, nrow))>0)]

filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split <- lapply(filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split, FUN=function(x){
  x <- x[which(x$pos_pairs_realprot_coord %in% df_hypergeom_test$pair_name[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1)]),]
  return(x)
})
filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split <- filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split[which(unlist(lapply(filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split, nrow))>0)]
##################################################################################################################################################################
#bind data together again
filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split_bound <- do.call(rbind, filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split)
#split by pdz mutation
spec_change_by_pdz_mut <- split.data.frame(filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split_bound, f=filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split_bound$pdz_mut1_realprot_coord)
#write df for cluster 3.0
#this is all spec changing cript/pdz
write_df_clustering <- function(x){
  x <- as.data.frame(x)
  intermed_df <- cbind(rownames(x), rownames(x),x)
  colnames(intermed_df)[1:2] <- c("GID","NAME")
  return(intermed_df)
  }
pdz3cript_dm_all_pos_epi <- pdz3cript_dm[(which(pdz3cript_dm$pdz_mut1_realprot_coord %in% pdz3cript_dm$pdz_mut1_realprot_coord[which(pdz3cript_dm$positive_epi==1)])),]
matrices_resid_for_export <- make_matrix_for_clustering(pdz3cript_dm_all_pos_epi, "pdz3cript_resids[pdz3cript_double_mut_pos]",filtering=F)
matrices_resid_for_export_df <- write_df_clustering(matrices_resid_for_export)
#write.table(matrices_resid_for_export_df, quote=F, sep="\t", row.names=F, "./fuzzy_specificity_data/figures/figure5_pdz_specificity/pdz3cript_dm_all_posepi.txt")
########################################################################################################################################################################
#heatmaps for fig 5 -- pos sig 
bfactor_vals_ligand <- read.delim("./fuzzy_specificity_data/figures/figure1_genetic_landscape/bfactor_5heb_ligand.txt")
unique_spec_changes_each_pos <- lapply(filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split, FUN=function(x){
  unique(x$cript_mut2_realprot_coord[which(x$positive_epi==1)])
  })
no_unique_spec_changes_each_pos <- lapply(unique_spec_changes_each_pos, length)
melted_no_unique_spec_changes_each_pos <- reshape2::melt(no_unique_spec_changes_each_pos)
melted_no_unique_spec_changes_each_pos$L1 <- factor(melted_no_unique_spec_changes_each_pos$L1, levels=as.character(seq(-5,0,1)))
melted_no_unique_spec_changes_each_pos$bfactor <- bfactor_vals_ligand$B_factor[which(bfactor_vals_ligand$residue_NtoC %in% melted_no_unique_spec_changes_each_pos$L1)]
no_unique_spec_changes_cript_bfactor <- ggplot(data=melted_no_unique_spec_changes_each_pos, aes(x=bfactor, y=value, fill=L1))+
  #geom_bar(stat="identity", position=position_dodge(), fill="white", color="black")+
  geom_point(shape=21, size=2)+
  scale_fill_manual(values=as.vector(cript_cols_matrix[3:length(cript_cols_matrix)]), limits=factor(-as.numeric(unlist(lapply(names(cript_cols_matrix[3:length(cript_cols_matrix)]), str_sub,-1,-1)))))+
  theme_bw()+
  xlab("b-factor")+
  ylab("Number of unique specificity changes (/19)")+
  scale_y_continuous(breaks=seq(0,20,5))
#pdz3cript_dm with any spec changing pdz mutation by all cript mutations to put into cluster 3.0
pdz3cript_dm_split_criptmut_posepi_bound <- do.call(rbind, pdz3cript_dm_split_criptmut_posepi)
matrices_resid_for_export_bound <- make_matrix_for_clustering(pdz3cript_dm_split_criptmut_posepi_bound, "pdz3cript_resids[pdz3cript_double_mut_pos]", filtering=F)
matrices_resid_for_export_bound_df <- write_df_clustering(matrices_resid_for_export_df)
matrices_resid_for_export_bound_df <- matrices_resid_for_export_bound_df[,-c(1:2)]
write.table(matrices_resid_for_export_bound_df, quote=F, sep ="\t", row.names = F,
            "./fuzzy_specificity_data/figures/figure5_pdz_specificity/all_spec_changers_resid_matrix.txt")

#output the number of specificity-changing mutations per pdz-cript residue in a table
write.table(df_hypergeom_test, "./fuzzy_specificity_data/figures/figure5_pdz_specificity/num_specificity_muts_per_pair_pdz_cript.txt", sep="\t", row.names = F, quote=F)
#######################################################################################################################
#add the wildtype to mutant-mutant heatmaps
#pdz_mut1_realprot_coord~cript_mut2_realprot_coord
#for the consolidated wt pdz; 
pdz_wt_row_for_dm_matrices <- pdz_dimsum_wtpdz_consolidated
pdz_wt_row_for_dm_matrices$pdz_mut1_realprot_coord <- "WT"
pdz_wt_row_for_dm_matrices$cript_mut2_realprot_coord <- " "
diffs_cript_wt_pdz <- lapply(str_sub(pdz_wt_row_for_dm_matrices$aa_seq,-8,-1), list.string.diff, "KNYKQTSV")
pdz_wt_row_cript_mut2_realprot_coord <- NULL
for (i in 1:length(diffs_cript_wt_pdz)){
  pdz_wt_row_cript_mut2_realprot_coord[[i]] <- paste(diffs_cript_wt_pdz[[i]][1,]$poly.seq.b, diffs_cript_wt_pdz[[i]][1,]$position, diffs_cript_wt_pdz[[i]][1,]$poly.seq.a, sep="")
}
pdz_wt_row_for_dm_matrices$cript_mut2_realprot_coord <- unlist(lapply(pdz_wt_row_cript_mut2_realprot_coord, change_mut_numbers,-8))
pdz_wt_row_for_dm_matrices$cript_mut2_realprot_coord[which(pdz_wt_row_for_dm_matrices$WT==T)] <- "WT"
pdz_wt_row_for_dm_matrices <- pdz_wt_row_for_dm_matrices[-which(pdz_wt_row_for_dm_matrices$cript_mut2_realprot_coord=="NNAA"),]
#do the same for the wt cript;
cript_wt_col_for_dm_matrices <- pdz_dimsum_wtcript
cript_wt_col_for_dm_matrices$pdz_mut1_realprot_coord <- " "
cript_wt_col_for_dm_matrices$cript_mut2_realprot_coord <- "WT"
diffs_pdz_wtcript <- lapply(str_sub(cript_wt_col_for_dm_matrices$aa_seq,1,100), list.string.diff, pdz_aa_seq_wt)
cript_wt_col_pdz_mut1_realprot_coord <- NULL
for (i in 1:length(diffs_pdz_wtcript)){
  cript_wt_col_pdz_mut1_realprot_coord[[i]] <- paste(diffs_pdz_wtcript[[i]][1,]$poly.seq.b, diffs_pdz_wtcript[[i]][1,]$position, diffs_pdz_wtcript[[i]][1,]$poly.seq.a, sep="")
}
cript_wt_col_for_dm_matrices$pdz_mut1_realprot_coord <- unlist(lapply(cript_wt_col_pdz_mut1_realprot_coord, change_mut_numbers,+302))
cript_wt_col_for_dm_matrices$pdz_mut1_realprot_coord[which(cript_wt_col_for_dm_matrices$WT==T)] <- "WT"
cript_wt_col_for_dm_matrices <- cript_wt_col_for_dm_matrices[-which(cript_wt_col_for_dm_matrices$pdz_mut1_realprot_coord=="NNAA"),]
###################################################
#also add white wt row for the other two heatmaps
#Figure 5 heatmaps
make_resid_binding_mats_to_clust_w_fdr_layer_wwt_blank_row <- function(x){
  resids_to_cluster <- make_matrix_for_clustering(x, "pdz3cript_resids[pdz3cript_double_mut_pos]", filtering=T)
  mat_fdr_values <- make_matrix_for_clustering(x, "fdr_dm", filtering=F)
  paletteLength <- 50
  myColor <- colorRampPalette(c("midnightblue", "white", "orangered"))(paletteLength)
  #diff color palette
  myColor_v2 <- colorRamp2(seq(-0.6, 0.6, length=3), c("midnightblue", "white", "orangered"))
  #resids_to_cluster <- Heatmap(resids_to_cluster)
  resids_to_cluster <- resids_to_cluster[,index_x_to_y(colnames(resids_to_cluster), ordering_aa_muts_cript)]
  resids_to_cluster <- resids_to_cluster[index_x_to_y(rownames(resids_to_cluster), ordering_aa_muts_pdz),]
  resids_to_cluster_ordered <- resids_to_cluster
  resids_to_cluster_ordered <- rbind(rep(0, ncol(resids_to_cluster_ordered)), resids_to_cluster_ordered)
  rownames(resids_to_cluster_ordered)[1] <- "WT"
  #resids_to_cluster_ordered <- resids_to_cluster@matrix[row_order(resids_to_cluster), column_order(resids_to_cluster)]
  fdr_matrix_ordered <- mat_fdr_values[rownames(resids_to_cluster_ordered)[-1], colnames(resids_to_cluster_ordered)]
  fdr_matrix_ordered <- rbind(rep(100, ncol(fdr_matrix_ordered)), fdr_matrix_ordered)
  rownames(fdr_matrix_ordered)[1] <- "WT"
  mut_physicochem_col <- unlist(lapply(colnames(resids_to_cluster_ordered), str_sub, -1, -1))
  mut_physicochem_row <- unlist(lapply(rownames(resids_to_cluster_ordered), str_sub, -1,-1))
  mut_position_col <- unlist(lapply(colnames(resids_to_cluster_ordered), str_sub,1,-2))
  mut_position_row <- unlist(lapply(rownames(resids_to_cluster_ordered), str_sub,1,-2))
  col_annot_mat <- HeatmapAnnotation(cn=anno_simple(mut_physicochem_col, col=zappo_cols_matrix),
                                     annotation_name_gp= gpar(fontsize = 10),
                                     annotation_name_side = "left",
                                     show_annotation_name = F)
  pdz_cols_matrix_expanded <- c("white", pdz_cols_matrix)
  names(pdz_cols_matrix_expanded)[1] <- "W"
  row_annot_mat <- rowAnnotation(cn=anno_simple(mut_physicochem_row, col=zappo_cols_matrix),
                                 #cn2=anno_simple(mut_position_row, col=pdz_cols_matrix_expanded),
                                 annotation_name_gp= gpar(fontsize = 10),
                                 show_annotation_name = F, simple_anno_size=unit(0.4, "cm"))
  resid_heatmap <- Heatmap(resids_to_cluster_ordered, col=myColor_v2, name = "binding residual", column_names_gp = grid::gpar(fontsize=8), row_names_gp = grid::gpar(fontsize=8), row_names_side = "left",
                           show_column_dend = F, show_row_dend = F,
                           cluster_rows = F, cluster_columns = F,
                           bottom_annotation = col_annot_mat,
                           left_annotation = row_annot_mat,
                           width = ncol(resids_to_cluster)*unit(3, "mm"), 
                           height = nrow(resids_to_cluster)*unit(3, "mm"),
                           show_row_names = T,
                           layer_fun = function(j, i, x, y, width, height, fill) {
                             #v= pindex(fdr_matrix_ordered, i,j)
                             v=fdr_matrix_ordered
                             l= v > .05 & v < .10   
                             l2= v > .01 & v < .05
                             l3= v > 0 & v < .01
                             grid.text("*", x[l], y[l], gp = gpar(fontsize = 8))
                             grid.text("**", x[l2], y[l2], gp = gpar(fontsize = 8))
                             grid.text("***", x[l3], y[l3], gp = gpar(fontsize = 8))
                           })
  bindingfit_to_cluster <- make_matrix_for_clustering(x, "wt_pdz_mut_pdz_diff", filtering=T)
  #reorder the matrix since it won't be clustered, just want to show annotations
  #order by cript residue/pdz residue and by physicochemical AA properties
  bindingfit_to_cluster <- bindingfit_to_cluster[,index_x_to_y(colnames(bindingfit_to_cluster), ordering_aa_muts_cript)]
  bindingfit_to_cluster <- bindingfit_to_cluster[index_x_to_y(rownames(bindingfit_to_cluster), ordering_aa_muts_pdz),]
  bindingfit_to_cluster <- rbind(rep(0, ncol(bindingfit_to_cluster)), bindingfit_to_cluster)
  rownames(bindingfit_to_cluster)[1] <- "WT"
  #diff color palette
  myColor_v2 <- colorRamp2(c(-1,0, 1.05), c("midnightblue","white", "orangered"))
  mut_physicochem_col <- unlist(lapply(colnames(bindingfit_to_cluster), str_sub, -1, -1))
  mut_physicochem_row <- unlist(lapply(rownames(bindingfit_to_cluster), str_sub, -1,-1))
  mut_position_col <- unlist(lapply(colnames(bindingfit_to_cluster), str_sub,1,-2))
  mut_position_row <- unlist(lapply(rownames(bindingfit_to_cluster), str_sub,1,-2))
  col_annot_mat <- HeatmapAnnotation(cn=anno_simple(mut_physicochem_col, col=zappo_cols_matrix),
                                     #cn2=anno_simple(mut_position_col, col=cript_cols_matrix),
                                     annotation_name_gp= gpar(fontsize = 10),
                                     annotation_name_side = "left",
                                     show_annotation_name = F,
                                     simple_anno_size = unit(0.6,"cm"))
  wt_binding_heatmap <- Heatmap(bindingfit_to_cluster, col=myColor_v2, name = "MUT-WT", 
                                column_names_gp = grid::gpar(fontsize=8), 
                                width = ncol(bindingfit_to_cluster)*unit(3, "mm"), 
                                height = nrow(bindingfit_to_cluster)*unit(3, "mm"),
                                layer_fun = function(j, i, x, y, width, height, fill) {
                                  #v= pindex(fdr_matrix_ordered, i,j)
                                  v=fdr_matrix_ordered
                                  l= v < .10   
                                  #l2= v > .01 & v < .05
                                 # l3= v > 0 & v < .01
                                 # grid.text("*", x[l], y[l], gp = gpar(fontsize = 8))
                                #  grid.text("**", x[l2], y[l2], gp = gpar(fontsize = 8))
                                 # grid.text("***", x[l3], y[l3], gp = gpar(fontsize = 8))
                                  grid.rect(x = x[l], y = y[l], width = width, height = height, 
                                            gp = gpar(col = "black", fill = NA, lwd=0.25))
                                },
                                  show_column_names = T, show_row_names = F,
                                show_column_dend = F, show_row_dend = F,
                                bottom_annotation = col_annot_mat,
                                #right_annotation = row_annot_mat,
                                clustering_distance_rows = "pearson",
                                clustering_distance_columns = "pearson",
                                cluster_rows = F, row_names_gp = gpar(fontsize=8), 
                                cluster_columns = F, na_col = "#EBEBEB")
  #################optional/for supp:
  #bring in the wildtype clustering info
  x_subset_cols <- x[,colnames(pdz_wt_row_for_dm_matrices)]
  z <- rbind(x_subset_cols,pdz_wt_row_for_dm_matrices)
  criptwts_to_add <- cript_wt_col_for_dm_matrices[which(cript_wt_col_for_dm_matrices$pdz_mut1_realprot_coord %in% x_subset_cols$pdz_mut1_realprot_coord),]
  #also add in the WT
  criptwts_to_add <- rbind(criptwts_to_add,cript_wt_col_for_dm_matrices[which(cript_wt_col_for_dm_matrices$WT==T),])
  criptwts_to_add$cript_mut2_realprot_coord <- paste(str_sub(z$cript_mut2_realprot_coord[1],1,1), -get_mut_pos(z$cript_mut2_realprot_coord[1]), str_sub(z$cript_mut2_realprot_coord[1],1,1), sep ="")
  criptwts_to_add_subset_cols <- criptwts_to_add[,colnames(pdz_wt_row_for_dm_matrices)]
  z <- rbind(z, criptwts_to_add_subset_cols)
  rawfitness_to_cluster <- make_matrix_for_clustering(z, "fitness", filtering=T)
  #reorder the matrix since it won't be clustered, just want to show annotations
  #order by cript residue/pdz residue and by physicochemical AA properties
  rawfitness_to_cluster <- rawfitness_to_cluster[,index_x_to_y(colnames(rawfitness_to_cluster), ordering_aa_muts_cript)]
  #for one position (T-2) in a couple of mutations/rows, there are just under 50% NAs which don't get filtered in raw binding but do get filtered in residuals
  #this is to make sure the rownames match across the 3 heatmaps
  rawfitness_to_cluster <- rawfitness_to_cluster[rownames(bindingfit_to_cluster),]
  rawfitness_to_cluster <- rawfitness_to_cluster[index_x_to_y(rownames(rawfitness_to_cluster), c("WT",ordering_aa_muts_pdz)),]
  #need a diff fdr palette for this heatmap since it has additional columns:
  #rownames(rawfitness_to_cluster)[-which(rownames(rawfitness_to_cluster)=="WT")]
#find WT col in cript
which_wt_col <- function(cript_colnames){str_sub(cript_colnames,1,1)==str_sub(cript_colnames,-1,-1)}
cript_wt_col <- colnames(rawfitness_to_cluster)[which(unlist(lapply(colnames(rawfitness_to_cluster), which_wt_col))==T)]
fdr_matrix_ordered_rawbinding <- mat_fdr_values[rownames(rawfitness_to_cluster)[-which(rownames(rawfitness_to_cluster)=="WT")], colnames(rawfitness_to_cluster)[-which(colnames(rawfitness_to_cluster)==cript_wt_col)]]
#add in wt row (blank)
fdr_matrix_ordered_rawbinding <- rbind(rep(100, ncol(fdr_matrix_ordered_rawbinding)), fdr_matrix_ordered_rawbinding)
rownames(fdr_matrix_ordered_rawbinding)[1] <- "WT"
#add in wt col (blank)
fdr_matrix_ordered_rawbinding <- cbind(rep(100, nrow(fdr_matrix_ordered_rawbinding)), fdr_matrix_ordered_rawbinding)
colnames(fdr_matrix_ordered_rawbinding)[1] <- cript_wt_col
fdr_matrix_ordered_rawbinding <- fdr_matrix_ordered_rawbinding[,colnames(rawfitness_to_cluster)]
  #diff color palette
  myColor_v2 <- colorRamp2(c(-1.2,dead_var_cutoff, 0.2), c("midnightblue","white", "orangered"))
  mut_physicochem_col <- unlist(lapply(colnames(rawfitness_to_cluster), str_sub, -1, -1))
  mut_physicochem_row <- unlist(lapply(rownames(rawfitness_to_cluster), str_sub, -1,-1))
  mut_position_col <- unlist(lapply(colnames(rawfitness_to_cluster), str_sub,1,-2))
  mut_position_row <- unlist(lapply(rownames(rawfitness_to_cluster), str_sub,1,-2))
  col_annot_mat <- HeatmapAnnotation(cn=anno_simple(mut_physicochem_col, col=zappo_cols_matrix),
                                     #cn2=anno_simple(mut_position_col, col=cript_cols_matrix),
                                     annotation_name_gp= gpar(fontsize = 10),
                                     annotation_name_side = "left",
                                     show_annotation_name = F,
                                     simple_anno_size = unit(0.4,"cm"))
  row_annot_mat_2 <- rowAnnotation(#cn=anno_simple(mut_physicochem_row, col=zappo_cols_matrix),
                                 cn2=anno_simple(mut_position_row, col=pdz_cols_matrix_expanded),
                                 annotation_name_gp= gpar(fontsize = 10),
                                 show_annotation_name = F, simple_anno_size=unit(0.4, "cm"))
  
  #expand the pdz annotation matrix
  raw_binding_heatmap <- Heatmap(rawfitness_to_cluster, col=myColor_v2, name = "binding", 
                                 column_names_gp = grid::gpar(fontsize=8), 
                                 width = ncol(rawfitness_to_cluster)*unit(3, "mm"), 
                                 height = nrow(rawfitness_to_cluster)*unit(3, "mm"),
                                 layer_fun = function(j, i, x, y, width, height, fill) {
                                   #v= pindex(fdr_matrix_ordered, i,j)
                                   v=fdr_matrix_ordered_rawbinding
                                   l= v < .10   
                                   grid.rect(x = x[l], y = y[l], width = width, height = height, 
                                             gp = gpar(col = "black", fill = NA, lwd=0.25))
                                 },
                                 show_column_names = T, show_row_names = T,
                                 show_column_dend = F, show_row_dend = F,
                                 bottom_annotation = col_annot_mat,
                                 right_annotation = row_annot_mat_2,
                                 clustering_distance_rows = "pearson",
                                 clustering_distance_columns = "pearson",
                                 cluster_rows = F, row_names_gp = gpar(fontsize=8), 
                                 cluster_columns = F, na_col = "#EBEBEB")
  resid_heatmap+wt_binding_heatmap+raw_binding_heatmap
}
list_of_pos_epi_mutwt_maps_wwt_blank <- lapply(filtered_pdz3cript_dm_pos_epi_atleast_one_row_cript_split,make_resid_binding_mats_to_clust_w_fdr_layer_wwt_blank_row)
