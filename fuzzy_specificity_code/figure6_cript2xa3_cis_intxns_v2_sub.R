################################################################################################
#Source code for Figure6 (+ related supplementary material) from 
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
#cript cis double mutant library analyses
#data quality etc
cript_2xa3 <- read.delim("./fuzzy_specificity_data/cript2xa3_norm_for_mochi.tsv")
cript_2xa3_with_syn <- read.delim("./fuzzy_specificity_data/cript_2xa3_normalized_pdz3_cript_sm.txt")
cript2xa3_aa_ham <- make_aaham_plot(cript_2xa3_with_syn)
cript_2xa3_repcorr <- get_rep_corr_matrix(cript_2xa3_with_syn)
cript_2xa3_fitcorr_mat <- ggpairs(cript_2xa3_with_syn[,grep("fitness[0-9]",colnames(cript_2xa3_with_syn))],
                                  lower = list(continuous = "hexbin", combo = "facethist", discrete = "facetbar", na =  "na"))+theme_bw()
densityplot_cript2xa3 <- make_gg_density_fit_general(cript_2xa3_with_syn)
#dead variant cutoff and justification
mean_stops <- mean(cript_2xa3$fitness[which(cript_2xa3$STOP==T)])
mode_stops <- as.numeric(names(table(cript_2xa3$fitness[which(cript_2xa3$STOP==T)])[order(table(cript_2xa3$fitness[which(cript_2xa3$STOP==T)]), decreasing=T)][1]))
median_stops <- median(cript_2xa3$fitness[which(cript_2xa3$STOP==T)])
dead_var_cutoff <- mean(cript_2xa3$fitness[which(cript_2xa3$STOP==T)])+2*sd(cript_2xa3$fitness[which(cript_2xa3$STOP==T)])
hist_dead_cutoff_cript2xa3 <- ggplot(cript_2xa3[which(cript_2xa3$STOP==T),], aes(x=fitness))+
  geom_histogram(fill="white", colour="black", bins=100)+
  xlab("Normalized fitness STOP variants")+
  geom_vline(xintercept = dead_var_cutoff, col=custom_scale_discrete_1[4], linetype="dashed")+
  geom_vline(xintercept = mean_stops, col=custom_scale_discrete_1[3], linetype="dashed")+theme_bw()
#read in mochi results
cript2xa3_mochi_path <- "./fuzzy_specificity_data/mochi_results/cript2xa3/task_1/predictions/predicted_phenotypes_all.txt"
cript_2xa3_mochi <- combine_read_delim_heldout_pred(cript2xa3_mochi_path)
#linear model
#lm_cript2xa3 <- get_linear_model(cript_2xa3_mochi, "first", "KNYKQTSV")
#cript_2xa3_lm <- make_gg_scatterplot(cript_2xa3_mochi$fitness,lm_cript2xa3$fitted.values, "Observed fitness", "Predicted fitness")
cript_2xa3_mochi_pred_perf <- make_gg_scatter_with_xy_line(cript_2xa3_mochi$fitness, cript_2xa3_mochi$predicted_phen, "Observed cript2xa3", "Predicted")
cript_2xa3_addtrait_fitness <- make_gg_scatterplot(apply(cript_2xa3_mochi[,49:58], 1, mean), cript_2xa3_mochi$fitness, "Additive trait cript2xa3", "Fitness")
#get residuals
cript2xa3_resids <- cript_2xa3_mochi$fitness - cript_2xa3_mochi$mean
cript2xa3_double_mut_pos <- which(cript_2xa3_mochi$Nham_aa==2)
cript2xa3_dm <- cbind(cript_2xa3_mochi[cript2xa3_double_mut_pos,], cript2xa3_resids[cript2xa3_double_mut_pos])
#compute double mutant positions/ids
cript2xa3_dm_diffs <- lapply(cript2xa3_dm$aa_seq, FUN=list.string.diff, "KNYKQTSV")
cript2xa3_dm_df_mut1 <- NULL
for (i in 1:length(cript2xa3_dm_diffs)){
  cript2xa3_dm_df_mut1[i] <- paste(cript2xa3_dm_diffs[[i]][1,]$poly.seq.b, cript2xa3_dm_diffs[[i]][1,]$position, cript2xa3_dm_diffs[[i]][1,]$poly.seq.a, sep="")
}
cript2xa3_dm_df_mut2 <- NULL
for (i in 1:length(cript2xa3_dm_diffs)){
  cript2xa3_dm_df_mut2[i] <- paste(cript2xa3_dm_diffs[[i]][2,]$poly.seq.b, cript2xa3_dm_diffs[[i]][2,]$position, cript2xa3_dm_diffs[[i]][2,]$poly.seq.a, sep="")
}
#double mut heatmap
cript2xa3_dm_df <- data.frame(cript2xa3_dm$`cript2xa3_resids[cript2xa3_double_mut_pos]`)
colnames(cript2xa3_dm_df) <- "residuals"
cript2xa3_dm_df$mut1 <- factor(cript2xa3_dm_df_mut1, levels=unique(cript2xa3_dm_df_mut1[order(unlist(lapply(cript2xa3_dm_df_mut1, get_mut_pos)))]))
cript2xa3_dm_df$mut2 <- factor(cript2xa3_dm_df_mut2, levels=unique(cript2xa3_dm_df_mut2[order(unlist(lapply(cript2xa3_dm_df_mut2, get_mut_pos)))]))
pheat_dm <- ggplot(cript2xa3_dm_df, aes(x=mut1, y=mut2))+
  geom_tile(aes(fill=residuals))+
  scale_fill_gradient2(low="midnightblue", mid="white", high="orangered", midpoint=0, name="residuals")+
  xlab("Mutation 1")+
  ylab("Mutation 2")+
  theme(axis.text.x = element_text(size = 6,
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 90),
        axis.text.y = element_text(size=6),
        legend.key.width = unit(0.3, "cm"),
        panel.border = element_rect(fill=NA,colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#double mut heatmap fitness
cript2xa3_dm_df_fit <- data.frame(cript2xa3_dm$fitness)
colnames(cript2xa3_dm_df_fit) <- "fitness"
cript2xa3_dm_df_fit$mut1 <- factor(cript2xa3_dm_df_mut1, levels=unique(cript2xa3_dm_df_mut1[order(unlist(lapply(cript2xa3_dm_df_mut1, get_mut_pos)))]))
cript2xa3_dm_df_fit$mut2 <- factor(cript2xa3_dm_df_mut2, levels=unique(cript2xa3_dm_df_mut2[order(unlist(lapply(cript2xa3_dm_df_mut2, get_mut_pos)))]))
pheat_dm_fit <- ggplot(cript2xa3_dm_df_fit, aes(x=mut1, y=mut2))+
  geom_tile(aes(fill=fitness))+
  scale_fill_gradient2(low="midnightblue", mid="white", high="orangered", midpoint=0, name="fitness")+
  xlab("CRIPT Mutation 1")+
  ylab("CRIPT Mutation 2")+
  theme(axis.text.x = element_text(size = 6,
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 90),
        axis.text.y = element_text(size=6),
        legend.key.width = unit(0.3, "cm"),
        panel.border = element_rect(fill=NA,colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#double mutant analysis -- looking for interactions
###################################################################
cript2xa3_dm$mut1 <- cript2xa3_dm_df_mut1
cript2xa3_dm$mut1_realprot_coord <- unlist(lapply(cript2xa3_dm$mut1, change_mut_numbers, -8))
cript2xa3_dm$mut2 <- cript2xa3_dm_df_mut2
cript2xa3_dm$mut2_realprot_coord <- unlist(lapply(cript2xa3_dm$mut2, change_mut_numbers, -8))
cript2xa3_dm$pos_pairs <- paste(unlist(lapply(cript2xa3_dm_df_mut1, get_mut_pos)), unlist(lapply(cript2xa3_dm_df_mut2, get_mut_pos)), sep="_")
cript2xa3_dm$pos_pairs_realprot_coord <- paste(as.numeric(unlist(lapply(cript2xa3_dm_df_mut1, get_mut_pos))-8), as.numeric(unlist(lapply(cript2xa3_dm_df_mut2, get_mut_pos))-8), sep="_")
cript2xa3_dm$mut1_mut2 <- paste(cript2xa3_dm$mut1, cript2xa3_dm$mut2, sep="_")
cript2xa3_dm$mut1_mut2_realprot_coord<- unlist(lapply(cript2xa3_dm$mut1_mut2, change_mut_numbers_double_mut_pos, -8))
cript2xa3_dm$mut1_mut2_realprotmut_coord<- unlist(lapply(cript2xa3_dm$mut1_mut2, change_mut_numbers_double_mut, -8))
#which are n-n, n-c, c-c (defined in basic funcs code file)
cript2xa3_dm$which_nc <- rep("", nrow(cript2xa3_dm))
cript2xa3_dm$which_nc[which(cript2xa3_dm$pos_pairs %in% cript_n_pos_pairs)] <- "n_n"
cript2xa3_dm$which_nc[which(cript2xa3_dm$pos_pairs %in% cript_c_pos_pairs)] <- "c_c"
cript2xa3_dm$which_nc[which(cript2xa3_dm$pos_pairs %in% cript_cross_pos_pairs)] <- "n_c"
#stats ! to see which residuals are significant
cript2xa3_dm$zscore_dm <- (cript2xa3_dm$`cript2xa3_resids[cript2xa3_double_mut_pos]`)/sqrt((cript2xa3_dm$sigma)^2)
cript2xa3_dm$ztest_dm <- 2*pnorm(-abs(cript2xa3_dm$zscore_dm)) #p-values 
cript2xa3_dm$fdr_dm <- p.adjust(cript2xa3_dm$ztest_dm, method="fdr", n=length(cript2xa3_dm$ztest_dm))
cript2xa3_dm$passes_fdr <- rep(0, nrow(cript2xa3_dm))
cript2xa3_dm$passes_fdr[which(cript2xa3_dm$fdr_dm<0.1)] <- 1
#plot the residuals
sc_plot_resids <- ggplot(cript2xa3_dm, aes(x=fitness, y=mean, fill=`cript2xa3_resids[cript2xa3_double_mut_pos]`))+
  geom_point(size=1.5, shape=21)+
  scale_fill_viridis()+
  theme(element_text(size=12, family="serif"))+
  xlab("Observed fitness")+
  ylab("Predicted fitness")+
  labs(fill="Residual")+
  theme_bw()+theme(legend.key.width = unit(0.3, "cm"))+geom_density2d(color="white")
#plot the FDR of residuals
sc_plot_fdr <- ggplot(cript2xa3_dm, aes(x=fitness, y=mean, fill=fdr_dm*100))+
  geom_point(size=1.5, shape=21)+
  scale_fill_viridis_c(direction=-1, option="B", begin= 0.1, end = 0.9)+
  theme(element_text(size=12, family="serif"))+
  xlab("Observed fitness")+
  ylab("Predicted fitness")+
  labs(fill="FDR of residual (%)")+
  theme_bw()+theme(legend.key.width = unit(0.3, "cm"))+geom_density2d(color="white")+
  geom_vline(xintercept=dead_var_cutoff, color="grey", linetype="dashed")
###
#plot them together
sc_plots_resid_fdr <- plot_grid(sc_plot_resids, sc_plot_fdr, nrow=2, ncol=1, align = "hv") 
#make a map of residuals vs fitness
resids_cript2xa3 <- make_gg_scatterplot(cript_2xa3_mochi$mean, cript2xa3_resids, "Predicted fitness", "Residual")
###
#testing for enrichment of spec. changing mutations (pos. epistasis)
cript2xa3_dm$positive_epi <- rep(0, nrow(cript2xa3_dm))
cript2xa3_dm$positive_epi[which(cript2xa3_dm$fitness>dead_var_cutoff & cript2xa3_dm$`cript2xa3_resids[cript2xa3_double_mut_pos]`>0 & cript2xa3_dm$passes_fdr==1)] <- 1
#creating an enrichment test for each pair of positions as compared to the rest of the dataset (with FDR calculated at level of whole dataset)
cript2xa3_dm_split <- split.data.frame(cript2xa3_dm, f=cript2xa3_dm$pos_pairs_realprot_coord)
#for each pair i also want to see if it's enriched in specificity changing mutations as compared to the whole dataset
hit_in_samp_per_pair <- lapply(cript2xa3_dm_split, FUN=function(x){sum(x$positive_epi)})
sample_size_per_pair <- lapply(cript2xa3_dm_split, FUN=function(x){nrow(x)})
hit_in_pop <- sum(do.call(rbind, cript2xa3_dm_split)$positive_epi)
fail_in_pop <- nrow(do.call(rbind, cript2xa3_dm_split))-hit_in_pop
hypergeom_test_per_pair <- phyper(unlist(hit_in_samp_per_pair)-1, hit_in_pop, fail_in_pop, unlist(sample_size_per_pair), lower.tail = F)
df_hypergeom_test <- data.frame(pair_name=names(hypergeom_test_per_pair),
                                type_intxn=unlist(lapply(cript2xa3_dm_split, FUN=function(x){unique(x$which_nc)})),
                                pvalue_phyper_test=as.numeric(hypergeom_test_per_pair),
                                pvalue_phyper_test_fdr=p.adjust(as.numeric(hypergeom_test_per_pair), method="fdr", n=length(cript2xa3_dm_split)),
                                hit_in_samp=as.numeric(unlist(hit_in_samp_per_pair)),
                                position_a=unlist(lapply(names(cript2xa3_dm_split), FUN=function(x){strsplit(x, "_")[[1]][1]})),
                                position_b=unlist(lapply(names(cript2xa3_dm_split), FUN=function(x){strsplit(x, "_")[[1]][2]})),
                                aa_pos_a=unlist(lapply(cript2xa3_dm_split, FUN=function(x){str_sub(unique(strsplit(x$mut1_mut2_realprot_coord, "_")[[1]][1]),1,-1)})),
                                aa_pos_b=unlist(lapply(cript2xa3_dm_split, FUN=function(x){str_sub(unique(strsplit(x$mut1_mut2_realprot_coord, "_")[[1]][2]), 1, -1)})))
bf_p_value_cutoff <- (0.05/length(cript2xa3_dm_split))
df_hypergeom_test$aa_pos_a <- factor(df_hypergeom_test$aa_pos_a, levels=c("K-7", "N-6", "Y-5", "K-4", "Q-3", "T-2", "S-1", "V0"))
df_hypergeom_test$aa_pos_b <- factor(df_hypergeom_test$aa_pos_b, levels=c("K-7", "N-6", "Y-5", "K-4", "Q-3", "T-2", "S-1", "V0"))
hypergeom_plot <- ggplot(data=df_hypergeom_test)+
  geom_point(aes(x=aa_pos_a, y=pvalue_phyper_test_fdr, fill=type_intxn), size=2, shape=21)+
  scale_fill_manual(values=custom_scale_discrete_5)+
  xlab("PDZ position")+
  ylab("-log10 pvalue")+
  #for only p-value cutoff passing stuff
  #labs(fill=paste("Pairs w p<", round(p_value_cutoff, 4), sep=""))+
  labs(fill=" ")+
  geom_text(data = subset(df_hypergeom_test, 0.01 > df_hypergeom_test$pvalue_phyper_test), aes(x=aa_pos_a, y=pvalue_phyper_test, label=paste(aa_pos_b, sep="")),vjust=1, hjust=1.5,nudge_x=0.5, size=2.5)+
  theme_bw()
hypergeom_plot+  
  theme(axis.text.x = element_text(angle = 90))
#geom_hline(yintercept=-log10(p_value_cutoff), color="grey", linetype="dashed")
#make a heatmap version of this
df_hypergeom_test$hit_in_samp_thresholded <- df_hypergeom_test$hit_in_samp
df_hypergeom_test$hit_in_samp_thresholded[which(df_hypergeom_test$hit_in_samp>10)] <- 10
heat_hitinsamp_matrix <- ggplot(df_hypergeom_test, aes(x=aa_pos_a, y=aa_pos_b))+
  geom_tile(aes(fill=hit_in_samp_thresholded))+
  scale_fill_gradient2(low="white", high="#5555ff", midpoint=0, name="No. of spec. changing muts", 
                       labels=c("0", "5","10+"), breaks=c(0,5,10))+
  #scale_fill_gradient2(low="midnightblue", mid="white", high="orangered", midpoint=0, name="No. of spec. changing muts")+
  xlab("CRIPT position")+
  ylab("CRIPT position")+
  #outline BF-corrected significant stuff
  #geom_tile(data=df_hypergeom_test[which(df_hypergeom_test$pvalue_phyper_test<p_value_cutoff),],aes(x=aa_pos_a, y=aa_pos_b), color="black", fill="transparent", linewidth=0.1)+
  #stuff that passes fdr<0.1
  geom_tile(data=df_hypergeom_test[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.1),],aes(x=aa_pos_a, y=aa_pos_b), color="black", fill="transparent", linewidth=0.1, linetype="dashed")+
  geom_tile(data=df_hypergeom_test[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.05),],aes(x=aa_pos_a, y=aa_pos_b), color="black", fill="transparent", linewidth=0.1, linetype="dashed")+
  geom_tile(data=df_hypergeom_test[which(df_hypergeom_test$pvalue_phyper_test_fdr<0.01),],aes(x=aa_pos_a, y=aa_pos_b), color="black", fill="transparent", linewidth=0.1)+
  #theme_bw()+
  theme(axis.text.x = element_text(size = 8,
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 90),
        axis.text.y = element_text(size=8),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "bottom",
        panel.border = element_rect(fill=NA,colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
###########figure 6
heat_hitinsamp_matrix
#general histogram to show number of pairs with significant residuals
no_sig_resids_cript2xa3 <- ggplot(df_hypergeom_test, aes(x=hit_in_samp))+
  geom_histogram(fill="white", colour="black", bins=10)+
  theme_bw()+
  xlab("Number of spec. changing mutations per position pair")+
  ylab("Count")
########################################################
#19x19 maps of the positions that have >>>spec. changing mutations
#for each double mutant pair, want to plot the 19x19 heatmap
make_matrix_for_clustering <- function(df_x, value_to_cluster, filtering){
  potential_matrix <- reshape2::dcast(df_x, mut1_realprot_coord~mut2_realprot_coord, value.var=as.character(value_to_cluster))
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
  myBreaks <- c(seq(min(x, na.rm=T), 0, length.out=ceiling(paletteLength/2)),
                seq(max(x, na.rm=T)/paletteLength, max(x, na.rm=T), length.out=floor(paletteLength/2)))
  return(myBreaks)
}
#want to quantify the number of significant interactions per pair
number_of_sig_pairs_10 <- lapply(cript2xa3_dm_split, FUN=function(x){length(which(x$fdr_dm < 0.10))})
ordered_sig_pairs_by_number_10 <- number_of_sig_pairs_10[order(unlist(number_of_sig_pairs_10), decreasing = T)]
#want to add an annotation layer in here where i can colour stuff by amino acid type
zappo_cols_matrix <- zappo_col_scheme$col
names(zappo_cols_matrix) <- zappo_col_scheme$letter
#alternatively, want to split up the matrix making and clustering so i can move between the two/plot multiple things on the same clustered heatmap:
make_mat_to_clust_w_fdr_layer <- function(x){
  mat_to_cluster <- make_matrix_for_clustering(x, "cript2xa3_resids[cript2xa3_double_mut_pos]", filtering=T)
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
  Heatmap(mat_to_cluster, col=myColor_v2, name = "binding residual", column_names_gp = grid::gpar(fontsize=8), row_names_gp = grid::gpar(fontsize=8),
          show_column_dend = F, show_row_dend = F,
          bottom_annotation = col_annot_mat,
          right_annotation = row_annot_mat,
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
heatmaps_sig_pos_w_fdr <- lapply(cript2xa3_dm_split[df_hypergeom_test$pair_name[order(df_hypergeom_test$hit_in_samp, decreasing=T)]][1:15], make_mat_to_clust_w_fdr_layer)
grobbed_hm_cols <- lapply(heatmaps_sig_pos_w_fdr, ggplotify::as.grob)
grid_hm_top_spec_changing <- cowplot::plot_grid(plotlist=grobbed_hm_cols[1:4], nrow=2, ncol=2)
#write out the number of spec. changing muts per pair in cript-cis interactions
write.table(df_hypergeom_test, "./fuzzy_specificity_data/figures/figure6_cript_cis_intxns/num_specificity_muts_per_pair_cript_cript.txt", sep="\t", row.names = F, quote=F)

#supplementary
#plot 1st order terms in this data vs. combo data 
dgs_cript2xa3 <- read.delim("./fuzzy_specificity_data/mochi_results/cript2xa3/task_1/weights/weights_Folding.txt")
dgs_combo <- read.delim("./fuzzy_specificity_data/mochi_results/cript_nc_balanced_exp10_1/task_1/weights/weights_Folding.txt")
matched_cript2xa3_combo_dgs <- match_datasets_dgs(dgs_cript2xa3, dgs_combo)
make_gg_scatter_with_xy_line(matched_cript2xa3_combo_dgs$x_mean_kcal.mol, matched_cript2xa3_combo_dgs$y_mean_kcal.mol, "cript2xa3 binding ddgs", "combinatorial binding ddgs")