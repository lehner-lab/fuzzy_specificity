library(ComplexHeatmap)
library(ggplot2)
library(stringr)
library(ggseqlogo)
library(GGally)
library(ggpubr)
library(viridis)
library(viridisLite)
library(seqinr)
library(cowplot)
library(data.table)
####################################################################################################################################
###############aesthetics###########################################################################################################
####################################################################################################################################
#some custom colour scales
custom_scale_discrete_1 <- c("midnightblue", "orangered", "darkolivegreen", "lightsalmon1")
custom_scale_discrete_2 <- c("lightslateblue", "darkseagreen1", "sienna1")
custom_scale_discrete_3 <- c("orangered", "midnightblue", "#F9CB64", "darkseagreen1")
custom_scale_discrete_4 <- c("#8C8CB7", "#F49F80")
custom_scale_discrete_5 <- c("lightsalmon", "#D0CEEF", "#F9CB64", "lightblue")
custom_scale_discrete_6 <- c(custom_scale_discrete_5, "#C3FF40", "#4072FF", "#81EFCA", "#492BF9", "#FFC55F", "orangered", "#FFAE55")
#for ggseqlogo
clustal_col_scheme = make_col_scheme(chars=c('A', 'I', 'L', 'M', 'F', 'W', 'V', 
                                             'D', 'E', 
                                             'N', 'Q', 'S', 'T',
                                             'R', 'K',
                                             'C',
                                             'G', 
                                             'H', 'Y',
                                             'P'),
                                     groups=c('hydrophobic', 'hydrophobic', 'hydrophobic', 'hydrophobic', 'hydrophobic', 'hydrophobic', 'hydrophobic',
                                              'negative charge', 'negative charge',
                                              'polar', 'polar', 'polar', 'polar',
                                              'positive charge', 'positive charge',
                                              'cysteine',
                                              'glycine',
                                              'aromatic', 'aromatic',
                                              'proline'), 
                                     cols=c('#383C82', '#383C82', '#383C82', '#383C82', '#383C82', '#383C82', '#383C82',
                                            '#FF6B3E', '#FF6B3E',
                                            '#2B6141', '#2B6141', '#2B6141', '#2B6141',
                                            '#C7D3F9', '#C7D3F9',
                                            '#FFCCFF',
                                            '#F9CB64',
                                            '#74B6B7', '#74B6B7',
                                            '#C1FF69'))
zappo_col_scheme = make_col_scheme(chars=c('A', 'I', 'L', 'M', 'V', 
                                           'D', 'E', 
                                           'N', 'Q', 'S', 'T',
                                           'R', 'K','H',
                                           'C',
                                           'G','P', 
                                           'Y', 'W', 'F'),
                                   groups=c('aliphatic/hydrophobic', 'aliphatic/hydrophobic', 'aliphatic/hydrophobic', 'aliphatic/hydrophobic', 'aliphatic/hydrophobic',
                                            'negative', 'negative',
                                            'hydrophilic', 'hydrophilic', 'hydrophilic', 'hydrophilic',
                                            'positive', 'positive','positive',
                                            'cysteine',
                                            'special conform.','special conform.',
                                            'aromatic', 'aromatic', 'aromatic'), 
                                   cols=c('#383C82', '#383C82', '#383C82', '#383C82', '#383C82',
                                          '#FF6B3E', '#FF6B3E',
                                          '#2B6141', '#2B6141', '#2B6141', '#2B6141',
                                          '#C7D3F9', '#C7D3F9','#C7D3F9',
                                          '#FFCCFF',
                                          '#F9CB64','#F9CB64',
                                          '#74B6B7', '#74B6B7','#74B6B7'))

####################################################################################################################################
###############sequence stuff#######################################################################################################
####################################################################################################################################
#compute variant diffs:
list.string.diff <- function(a, b, exclude = c("-", "?"), ignore.case = TRUE, show.excluded = FALSE)
{
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  if(ignore.case)
  {
    a <- toupper(a)
    b <- toupper(b)
  }
  split_seqs <- strsplit(c(a, b), split = "")
  only.diff <- (split_seqs[[1]] != split_seqs[[2]])
  only.diff[
    (split_seqs[[1]] %in% exclude) |
      (split_seqs[[2]] %in% exclude)
  ] <- NA
  diff.info<-data.frame(which(is.na(only.diff)|only.diff),
                        split_seqs[[1]][only.diff],split_seqs[[2]][only.diff])
  names(diff.info)<-c("position","poly.seq.a","poly.seq.b")
  if(!show.excluded) diff.info<-na.omit(diff.info)
  diff.info
}
#get position of mutation
get_mut_pos <- function(x){
  as.numeric(c2s(s2c(x)[grep("[0-9]", s2c(x))]))
}
#change mutation index, e.g. to match that of the real protein
change_mut_numbers <- function(x,y){
  mut_number <- as.numeric(c2s(s2c(x)[grep("[0-9]", s2c(x))]))
  aa_one <- str_sub(x,1,1)
  aa_last <- str_sub(x,-1,-1)
  paste(aa_one, mut_number+y, aa_last, sep="")
}
change_mut_numbers_pos <- function(x,y){
  mut_number <- as.numeric(c2s(s2c(x)[grep("[0-9]", s2c(x))]))
  aa_one <- str_sub(x,1,1)
  paste(aa_one, mut_number+y, sep="")
}
#change mutation index of a double mutant in e.g. "T6A_V8K" format; x is the mutant id, y is index to sum
change_mut_numbers_double_mut <- function(x,y){
  mut_number_1 <- strsplit(x, "_")[[1]][1]
  mut_number_2 <- strsplit(x, "_")[[1]][2]
  mut_number_1_new <- change_mut_numbers(mut_number_1, y)
  mut_number_2_new <- change_mut_numbers(mut_number_2, y)
  paste(mut_number_1_new, mut_number_2_new, sep = "_")
}
#change mutation index of double mutant position e.g. "T6_V8" format; x is the mutant id, y is index to sum
change_mut_numbers_double_mut_pos <- function(x,y){
  mut_number_1 <- strsplit(x, "_")[[1]][1]
  mut_number_2 <- strsplit(x, "_")[[1]][2]
  mut_number_1_new <- change_mut_numbers_pos(mut_number_1, y)
  mut_number_2_new <- change_mut_numbers_pos(mut_number_2, y)
  paste(mut_number_1_new, mut_number_2_new, sep = "_")
}
pdz_aa_seq_wt <- "GEEDIPREPRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQYKPEEYSRFEA"
pdz_b1_nuc_wt_full <- tolower("GGGGAGGAAGACATTCCCCGAGAACCGAGGCGAATTGTGATCCACCGGGGCTCCACGGGCCTGGGCTTCAACATCGTGGGTGGCGAGGACGGTGAAGGCATCTTCATCTCCTTTATCCTGGCCGGGGGCCCTGCAGACCTCAGTGGGGAG")
pdz_b2_nuc_wt_full <- tolower("CTGCGGAAGGGGGACCAGATCCTGTCGGTCAACGGTGTGGACCTCCGAAATGCCAGCCATGAGCAGGCTGCCATTGCCCTGAAGAATGCGGGTCAGACGGTCACGATCATCGCTCAGTATAAACCAGAAGAGTACAGCCGATTCGAGGCC")
cript_nuc_full <- tolower("AAAAACTACAAGCAAACATCTGTC")
new_order_aa <- c("*", "A", "M", "V", "L", "I", "C", "H", "R", "K","E", "D", "P", "G", "S", "T", "N", "Q", "F", "Y", "W")
order_cript <- c("K-7", "N-6", "Y-5", "K-4", "Q-3", "T-2", "S-1", "V0")
index_x_to_y <- function(x,y){order(match(x,y))}
#cript colouring
cript_cols_matrix <- viridis(length(order_cript), option = "B")
names(cript_cols_matrix) <- order_cript
#basic info for cript cis double mutant libraries
cript_n_pos_pairs <- c("1_2", "1_3", "1_4", "2_3", "2_4", "3_4")
cript_c_pos_pairs <- c("5_6", "5_7", "5_8", "6_7", "6_8", "7_8")
cript_cross_pos_pairs <- c("1_5","1_6","1_7","1_8","2_5", "2_6", "2_7", "2_8", "3_5", "3_6", "3_7", "3_8", "4_5", "4_6", "4_7", "4_8")
####################################################################################################################################
###############general calculations#################################################################################################
####################################################################################################################################
#get top percent of x variants (sequence)
top_percent_vars <- function(x, perc){
  x$aa_seq[which(x$fitness > quantile(x$fitness, prob=1-perc/100))]
}
#get bottom percent of x variants (sequence)
bottom_percent_vars <- function(x, perc){
  x$aa_seq[-which(x$fitness > quantile(x$fitness, prob=1-perc/100))]
}
#index number of variants with >10 counts in any input column
index_count_cols_input <- function(x){unique(which(x$input1>=10), which(x$input2>=10), which(x$input3>=10))}
####################################################################################################################################
###############model evaluation#####################################################################################################
####################################################################################################################################
#general function for evaluating model
#take out the wt and add use only held-out prediction
add_heldout_pred <- function(pred_file){
  pred_file_no_wt <- pred_file[-which(pred_file$WT=="True"),]
  predicted_phen <- NULL
  for (i in 1:nrow(pred_file_no_wt)){
    predicted_phen[i] <- (pred_file_no_wt[i,paste("fold_",pred_file_no_wt$Fold[i], sep="")])
  }
  pred_file_no_wt$predicted_phen <- predicted_phen
  return(pred_file_no_wt)
}

combine_read_delim_heldout_pred <- function(x){
  add_heldout_pred(read.delim(x))
}
get_first_order_coeffs <- function(ddg_file_loc){
  all_dgs <- read.delim(ddg_file_loc)
  sm_dgs <- all_dgs[which(nchar(all_dgs$id)==3),]$mean_kcal.mol
  return(sm_dgs)
}
get_linear_model <- function(x,y, wildtype){
  fit_s2c <- lapply(x$aa_seq, FUN="s2c")
  seq_data <- data.table(do.call(rbind, fit_s2c))
  my_data <- data.table(seq_data, fitness=x$fitness)
  #setting wt to 0.. also look into contrast argument in model.matrix
  my_data <- data.frame(my_data)
  wt_seq <- wildtype
  for (i in 1:nchar(wt_seq)){
    my_data[which(my_data[,i]==s2c(wt_seq)[i]), i] <- 0
  }
  my_data <- data.table(my_data)
  if(y=="first"){
    my_features <- model.matrix(fitness~., my_data)
  } 
  else{
    my_features <- model.matrix(fitness~.*., my_data)
  }
  my_modeldata <- data.table(my_features, fitness=my_data[,fitness])
  my_linearmodel <- lm(fitness~0+., my_modeldata)
  return(my_linearmodel)
}
####################################################################################################################################
###############heatmaps#############################################################################################################
####################################################################################################################################
new_order_aa <- c("*", "A", "M", "V", "L", "I", "C", "H", "R", "K","E", "D", "P", "G", "S", "T", "N", "Q", "F", "Y", "W")
order_cript <- c("K-7", "N-6", "Y-5", "K-4", "Q-3", "T-2", "S-1", "V0")
index_x_to_y <- function(x,y){order(match(x,y))}
#cript colouring
cript_cols_matrix <- viridis(length(order_cript), option = "B")
names(cript_cols_matrix) <- order_cript
####################################################################################################################################
###############plots################################################################################################################
####################################################################################################################################
make_gg_density_nor_fit_same_axis <- function(x){
  non_stop <- data.frame(fitness=x$fitness[-union(which(x$STOP==T), which(x$Nham==0))], 
                         type=rep("non_stop", length(x$fitness[-union(which(x$STOP==T), which(x$Nham==0))])),
                         block=x$data_id[-union(which(x$STOP==T), which(x$Nham==0))])
  non_stop$block <- gsub("cript_n", "(dynamic) CRIPT N", non_stop$block)
  non_stop$block <- gsub("cript_c", "(stable) CRIPT C", non_stop$block)
  all_stop <- data.frame(fitness=x$fitness[which(x$STOP==T)], type=rep("all_stop", length(x$fitness[which(x$STOP==T)])), block="STOP")
  all_syn <- data.frame(fitness=x$fitness[union(which(x$Nham_aa==0), which(x$WT==F))], type=rep("all_syn", length(x$fitness[union(which(x$Nham_aa==0), which(x$WT==F))])), block=rep("all_syn", length(x$fitness[union(which(x$Nham_aa==0), which(x$WT==F))])))
  wt <- data.frame(fitness=x$fitness[which(x$WT==T)], type=rep("wt", length(x$fitness[which(x$WT==T)])), block=rep("wt", length(x$fitness[which(x$WT==T)])))
  one_percent_n <- as.numeric(quantile(x$fitness[which(x$data_id=="cript_n")], prob=1-1/100))
  one_percent_c <- as.numeric(quantile(x$fitness[which(x$data_id=="cript_c")], prob=1-1/100))
  if(nrow(all_syn)>2){
    df_x <- rbind(non_stop, all_stop, all_syn)
  }else{
    df_x <- rbind(non_stop, all_stop)
  }
  p <- ggplot(df_x, aes(x=fitness, fill=block)) +
    geom_density(alpha=0.4) + 
    theme(element_text(size=12, family="serif"))+
    theme_bw() +
    scale_fill_manual(name=" ",values=custom_scale_discrete_3[c(2,1,3)], label=c("N", "C", "STOP"))+
    xlab("Normalized Fitness")+ylab("Density")+
    scale_x_continuous(limits = c(min(df_x$fitness),max(df_x$fitness)))
  p+ geom_vline(aes(xintercept=mean(wt$fitness),
                    color="WT"), linetype="dashed")+
    geom_vline(aes(xintercept=one_percent_n, color="Top_1_perc_N"), linetype="dashed")+
    geom_vline(aes(xintercept=one_percent_c, color="Top_1_perc_C"), linetype="dashed")+
    scale_color_manual(name=" ",values=c(WT="darkgrey", Top_1_perc_N="midnightblue", Top_1_perc_C="orangered"), breaks=c("WT", "Top_1_perc_N", "Top_1_perc_C"), labels=c("Wildtype", "Top 1% N", "Top 1% C"))
}

#get count replicate correlation matrix/plot of dimsum outputted/mochi outputted dataset
#replicate counts
get_rep_corr_matrix <- function(x){
  r_squared_mat <- (cor(x[,grep("count_e", colnames(x))], use = "pairwise.complete.obs"))
  alt_names_col_names_counts <- colnames(r_squared_mat)
  input_indices <- grep("s0", alt_names_col_names_counts)
  output_indices <- grep("s1", alt_names_col_names_counts)
  alt_names_col_names_counts[input_indices] <- gsub("count_e", "input ", alt_names_col_names_counts[input_indices])
  alt_names_col_names_counts[output_indices] <- gsub("count_e", "output ", alt_names_col_names_counts[output_indices])
  alt_names_col_names_counts <- gsub("_s[01]", "", alt_names_col_names_counts)
  if(identical(rownames(r_squared_mat),colnames(r_squared_mat))==T){
    colnames(r_squared_mat) <- alt_names_col_names_counts
    rownames(r_squared_mat) <- alt_names_col_names_counts
  }else{print("rows and cols don't match")}
  myColor <- colorRamp2(c(0,1),c("white", "orangered"))
  corr_heatmap <- Heatmap(r_squared_mat, col=myColor, heatmap_legend_param = list(title=expression("r")), name = deparse(substitute(x)),
                          width = ncol(r_squared_mat)*unit(7.5, "mm"),
                          height = nrow(r_squared_mat)*unit(7.5, "mm"),
                          cluster_rows=F, cluster_columns = F,rect_gp=gpar(type="none"), row_names_side = "left",
                          cell_fun=function(j, i, x, y, w, h, fill){
                            if(i>=j){
                              grid.rect(x, y, w, h, gp=gpar(fill=fill, col=fill))
                              grid.text(sprintf("%.2f", r_squared_mat[i,j]), x, y, gp=gpar(fontsize=8))
                            }
                          })
  return(corr_heatmap)
}
#for R^2 properly formatted in legend title; heatmap_legend_param = list(title=expression("R"^2))
#fitness uncorrelated
library(circlize)
get_fitness_corr_matrix <- function(x){
  r_squared_mat <- (cor(x[,grep("fitness[0-9]_uncorr", colnames(x))], use = "pairwise.complete.obs"))
  myColor <- colorRamp2(c(0,1),c("white", "orangered"))
  corr_heatmap <- Heatmap(r_squared_mat, col=myColor, heatmap_legend_param = list(title=expression("r")), name = deparse(substitute(x)),
                          cluster_rows=F, cluster_columns = F,rect_gp=gpar(type="none"), row_names_side = "left",
                          width = ncol(r_squared_mat)*unit(7.5, "mm"),
                          height = nrow(r_squared_mat)*unit(7.5, "mm"),
                          cell_fun=function(j, i, x, y, w, h, fill){
                            if(i>=j){
                              grid.rect(x, y, w, h, gp=gpar(fill=fill, col=fill))
                              grid.text(sprintf("%.2f", r_squared_mat[i,j]), x, y, gp=gpar(fontsize=8))
                            }
                          })
  return(corr_heatmap)
}
#aa hamming distance vs number of variants
make_aaham_plot <- function(x){
  df <- data.frame(table(x$Nham_aa))
  ggplot(df)+
    geom_bar(aes(x=Var1, y=Freq), stat="identity", fill="white", colour="black")+
    theme_classic()+
    scale_x_discrete(name="AA Hamming distance")+
    scale_y_continuous(name=(paste("No. of variants | total =", nrow(x), sep="")))
}
#make the ggpairs function able to plot geom hex rather than regular scatterplots
ggally_hexbin <- function (data, mapping, ...)  {
  p <- ggplot(data = data, mapping = mapping) + geom_hex(...) + scale_fill_viridis()
  p
}
#density scatterplot with xy line
make_gg_scatter_with_xy_line <- function(xaxis, yaxis, xname, yname){
  r2_value <- paste(round((cor(xaxis, yaxis)^2), digits=3), " | N =", length(xaxis), sep=" ")
  lm_data_plot <- ggplot()+
    geom_hex(aes(y=yaxis, x=xaxis), bins=60) +
    scale_fill_gradientn(colours = viridis(256, option = "B"))+
    geom_smooth(method='lm', formula=yaxis~xaxis) +
    xlab(as.character(xname)) + ylab(as.character(yname)) +
    theme(plot.title = element_text(size=12, family="serif")) + 
    theme_bw()+
    theme(legend.key.width = unit(0.3, "cm"))+
    ggtitle(bquote(R^2 ~ .(r2_value)))
  #trying to plot the residuals and then label the stuff with big residuals
  lm_data_plot + geom_abline(intercept=0, color="darkgrey", linewidth=0.9)
  #plot(lm_data)
}
#density scatterplot with no correlation
make_gg_scatterplot_no_corr <- function(xaxis, yaxis, xname, yname){
  ggplot()+
    geom_hex(aes(y=yaxis, x=xaxis), bins=60) +
    scale_fill_gradientn(colours = viridis(256, option = "B"))+
    xlab(as.character(xname)) + ylab(as.character(yname)) +
    xlab(as.character(xname)) + ylab(as.character(yname)) +
    theme(plot.title = element_text(size=12, family="serif")) + 
    theme_bw()+
    theme(legend.key.width = unit(0.3, "cm"))
}
make_gg_scatterplot <- function(xaxis, yaxis, xname, yname){
  r2_value <- paste(round((cor(xaxis, yaxis)^2), digits=3), " | N =", length(xaxis), sep=" ")
  ggplot()+
    geom_hex(aes(y=yaxis, x=xaxis), bins=60) +
    scale_fill_gradientn(colours = viridis(256, option = "B"))+
    xlab(as.character(xname)) + ylab(as.character(yname)) +
    theme(plot.title = element_text(size=12, family="serif")) + 
    theme_bw()+
    theme(legend.key.width = unit(0.3, "cm"))+
    ggtitle(bquote(R^2 ~ .(r2_value)))
}
#scatterplot with error bars and labels
make_scatterplot_errorbars <- function(df, xaxis, yaxis, lbl, error_x, error_y, xname, yname){
  scplot <- ggplot(df, aes(x=xaxis, y=yaxis, label=lbl))+
    geom_point(size=2)+
    geom_errorbar(aes(ymin=yaxis-error_y, ymax=yaxis+error_y))+
    geom_errorbar(aes(xmin=xaxis-error_x, xmax=xaxis+error_x))+
    xlab(as.character(xname))+
    ylab(as.character(yname))+
    theme(plot.title = element_text(size=12, family="serif")) + 
    theme_bw()+
    theme(legend.key.width = unit(0.3, "cm"))+
    geom_text(hjust=0, vjust=0, size=3)
  scplot + geom_abline(intercept=0, color="darkgrey", size=0.9)
}
#density plot for variant fitness
make_gg_density_nor_fit_general <- function(x){
  non_stop <- data.frame(fitness=x$nor_fitness[-union(which(x$STOP==T), which(x$Nham==0))], type=rep("non_stop", length(x$nor_fitness[-union(which(x$STOP==T), which(x$Nham==0))])))
  all_stop <- data.frame(fitness=x$nor_fitness[which(x$STOP==T)], type=rep("all_stop", length(x$nor_fitness[which(x$STOP==T)])))
  all_syn <- data.frame(fitness=x$nor_fitness[union(which(x$Nham_aa==0), which(x$WT==F))], type=rep("all_syn", length(x$nor_fitness[union(which(x$Nham_aa==0), which(x$WT==F))])))
  wt <- data.frame(fitness=x$nor_fitness[which(x$WT==T)], type=rep("wt", length(x$nor_fitness[which(x$WT==T)])))
  if(nrow(all_syn)>2){
    df_x <- rbind(non_stop, all_stop, all_syn)
  }else{
    df_x <- rbind(non_stop, all_stop)
  }
  p <- ggplot(df_x, aes(x=fitness, fill=type)) +
    geom_density(alpha=0.4) + 
    theme(element_text(size=12, family="serif"))+
    theme_bw() +
    scale_fill_manual(name=" ",values=custom_scale_discrete_1)+
    xlab("Fitness")+ylab("Density")+
    scale_x_continuous(limits = c(min(df_x$fitness),max(df_x$fitness)))
  #theme(text=element_text(size=9))
  p+ geom_vline(aes(xintercept=wt$fitness,
                    color="WT"), linetype="dashed")+
    scale_color_manual(name=" ",values=c(WT="black"), labels=c("Wildtype"))
}
#single mutant fitness heatmap
#wanna see a single mutant heat map also
make_pheat_fitness_sm <- function(x,y){
  single_mut_only <- x[which(x$Nham_aa==1),]
  #compute double mutant positions/ids
  single_mut_diffs <- lapply(single_mut_only$aa_seq, FUN=list.string.diff, x$aa_seq[which(x$WT==T)])
  single_mut_df_pos <- NULL
  for (i in 1:length(single_mut_diffs)){
    single_mut_df_pos[i] <- paste(single_mut_diffs[[i]][1,]$poly.seq.b, as.numeric(single_mut_diffs[[i]][1,]$position)+y, sep="")
  }
  single_mut_df_mut <- do.call(rbind,single_mut_diffs)$poly.seq.a
  single_mut_df <- data.frame(single_mut_only$fitness)
  colnames(single_mut_df) <- "fitness"
  single_mut_df$pos <- single_mut_df_pos
  single_mut_df$mut <- single_mut_df_mut
  #make wt part of the df so you can outline the wt residues on the heatmap
  wt_df_addition <- data.frame(fitness=x$fitness[which(x$WT==T)][1], pos=unique(single_mut_df$pos), mut=unlist(lapply(unique(single_mut_df$pos), str_sub,1,1)))
  single_mut_df <- rbind(single_mut_df, wt_df_addition)
  is_wt <- rep(0, nrow(single_mut_df))
  for(i in 1:nrow(single_mut_df)){is_wt[i] <- if(str_sub(single_mut_df$pos[i],1,1)==single_mut_df$mut[i]){1}else{0}}
  single_mut_df$is_wt <- is_wt
  single_mut_df$mut <- factor(single_mut_df$mut, levels=new_order_aa)
  single_mut_df$pos <- factor(single_mut_df$pos, levels=paste(s2c(x$aa_seq[which(x$WT==T)]), (1+y):(nchar(x$aa_seq[1])+y), sep=""))
  pheat_sm <- ggplot(single_mut_df,aes(x=pos,y=mut))+
    geom_tile(aes(fill=fitness), )+
    scale_fill_gradient2(low="midnightblue", mid="white", high="sienna1", midpoint=0, name="fitness")+
    geom_tile(data=single_mut_df[which(single_mut_df$is_wt==1),],aes(x=pos, y=mut), color="black", fill="transparent", linewidth=0.1)+
    theme(axis.text.x = element_text(size = 8,
                                     vjust = 0.5,
                                     hjust = 0.5,
                                     angle = 90),
          axis.text.y = element_text(size=10),
          legend.key.width = unit(0.3, "cm"),
          panel.border = element_rect(fill=NA,colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  pheat_sm
}
#make single mutant matrix out of fitness df
make_matrix_fitness_sm <- function(x,y){
  single_mut_only <- x[which(x$Nham_aa==1),]
  #compute double mutant positions/ids
  single_mut_diffs <- lapply(single_mut_only$aa_seq, FUN=list.string.diff, x$aa_seq[which(x$WT==T)])
  single_mut_df_pos <- NULL
  for (i in 1:length(single_mut_diffs)){
    single_mut_df_pos[i] <- paste(single_mut_diffs[[i]][1,]$poly.seq.b, as.numeric(single_mut_diffs[[i]][1,]$position)+y, sep="")
  }
  single_mut_df_mut <- do.call(rbind,single_mut_diffs)$poly.seq.a
  single_mut_df <- data.frame(single_mut_only$fitness)
  colnames(single_mut_df) <- "fitness"
  single_mut_df$pos <- single_mut_df_pos
  single_mut_df$mut <- single_mut_df_mut
  single_mut_df$mut <- factor(single_mut_df$mut, levels=new_order_aa)
  single_mut_df$pos <- factor(single_mut_df$pos, levels=paste(s2c(x$aa_seq[which(x$WT==T)]), (1+y):(nchar(x$aa_seq[1])+y), sep=""))
  return(single_mut_df)
}
make_gg_density_fit_general <- function(x){
  non_stop <- data.frame(fitness=x$fitness[-union(which(x$STOP==T), which(x$Nham==0))], type=rep("non_stop", length(x$fitness[-union(which(x$STOP==T), which(x$Nham==0))])))
  all_stop <- data.frame(fitness=x$fitness[which(x$STOP==T)], type=rep("all_stop", length(x$fitness[which(x$STOP==T)])))
  all_syn <- data.frame(fitness=x$fitness[union(which(x$Nham_aa==0), which(x$WT==F))], type=rep("all_syn", length(x$fitness[union(which(x$Nham_aa==0), which(x$WT==F))])))
  wt <- data.frame(fitness=x$fitness[which(x$WT==T)], type=rep("wt", length(x$fitness[which(x$WT==T)])))
  if(nrow(all_syn)>2){
    df_x <- rbind(non_stop, all_stop, all_syn)
  }else{
    df_x <- rbind(non_stop, all_stop)
  }
  p <- ggplot(df_x, aes(x=fitness, fill=type)) +
    geom_density(alpha=0.4) + 
    theme(element_text(size=12, family="serif"))+
    theme_bw() +
    scale_fill_manual(name=" ",values=custom_scale_discrete_1)+
    xlab("Fitness")+ylab("Density")+
    scale_x_continuous(limits = c(min(df_x$fitness),max(df_x$fitness)))
  #theme(text=element_text(size=9))
  p+ geom_vline(aes(xintercept=wt$fitness,
                    color="WT"), linetype="dashed")+
    scale_color_manual(name=" ",values=c(WT="black"), labels=c("Wildtype"))
}
####################################################################################################################################
###############matching################################################################################################################
####################################################################################################################################
match_datasets <- function(x, y){
  x_matched_indices <- match(intersect(x$aa_seq, y$aa_seq), x$aa_seq)
  y_matched_indices <- match(intersect(x$aa_seq, y$aa_seq), y$aa_seq)
  if(identical(x$aa_seq[x_matched_indices],y$aa_seq[y_matched_indices])==T){
    df_x_y <- data.frame(x_fitness=x$fitness[x_matched_indices], y_fitness=y$fitness[y_matched_indices])
    return(df_x_y)
  }else{print("x and y do not match, plz try again")}
}
match_datasets_dgs <- function(x, y){
  x_matched_indices <- match(intersect(x$id, y$id), x$id)
  y_matched_indices <- match(intersect(x$id, y$id), y$id)
  if(identical(x$id[x_matched_indices],y$id[y_matched_indices])==T){
    df_x_y <- data.frame(lab=x$id[x_matched_indices],
                         x_mean_kcal.mol=x$mean_kcal.mol[x_matched_indices], 
                         y_mean_kcal.mol=y$mean_kcal.mol[y_matched_indices],
                         x_95ci=x$ci95_kcal.mol[x_matched_indices],
                         y_95ci=y$ci95_kcal.mol[y_matched_indices])
    return(df_x_y)
  }else{print("x and y do not match, plz try again")}
}
####################################################################################################################################
###############other################################################################################################################
####################################################################################################################################

make_wt_col_na_free <- function(x){x$WT[which(is.na(x$WT))] <- ""; return(x)}
