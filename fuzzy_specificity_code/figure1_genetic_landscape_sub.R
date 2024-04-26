################################################################################################
#Source code for Figure1 (+ related supplementary material) from 
#"A complete map of specificity encoding for a partially fuzzy protein interaction"
#Zarin, T and Lehner, B. 2024
################################################################################################
#please download required data (Zenodo link) + see instructions/required packages at 
#https://github.com/lehner-lab/fuzzy_specificity/
################################################################################################
#remember to set directory 
#e.g. setwd("dir_containing_fuzzy_specificity_code_dir_AND_fuzzy_specificity_data_dir")
################################################################################################
#load basic functions
source("./fuzzy_specificity_code/source_functions/basic_functions_fuzzy_specificity_TZ.R")
#PDB analysis/structure
#analyze pdb files
#based on analyze_pdb_files.R -- my v1 script
library("bio3d")
home_dir <- getwd()
setwd("./fuzzy_specificity_data/pdb_pdz_ligands/")
pdb_files <- list.files(pattern="*.ent.gz", recursive = T)
read_pdb_files <- lapply(pdb_files, read.pdb, verbose=F)
setwd(home_dir)
chain_lengths <- unlist(lapply(read_pdb_files, FUN=function(x){
  length(unique(x$atom$chain))
}))
barplot(table(chain_lengths), col="white", ylab="Number of PDB entries", xlab="Number of chains", main="PDZ family PDB entries")

#filtering for things that have 2 chains, looking at B factor as you go along the residue
pdbs_chain_length_2 <- which(chain_lengths==2)
read_pdb_files_chain_length_2 <- read_pdb_files[pdbs_chain_length_2]
pdb_bfactor_values <- lapply(read_pdb_files_chain_length_2, FUN=function(x){
  x$atom$b[intersect(which(x$atom$chain=="B"),which(x$calpha))]
})
names(pdb_bfactor_values) <- unlist(lapply(basename.pdb(pdb_files[pdbs_chain_length_2], ext=".ent.gz"), FUN=function(x){strsplit(x, "pdb")[[1]][2]}))
#things that actually have bfactors
pdb_bfactor_values_wbfactor <- pdb_bfactor_values[-which(lapply(pdb_bfactor_values, sum)==0)]
ordered_pdb_bfactor_values_wbfactor <- pdb_bfactor_values_wbfactor[order(unlist(lapply(pdb_bfactor_values_wbfactor, length)), decreasing=T)]
#scale it to min b factor
scaled_ordered_pdb_bfactor <- lapply(ordered_pdb_bfactor_values_wbfactor, FUN=function(x){as.numeric(x)/min(as.numeric(x))})

#chunk up the data into quantiles based on length
bfact_sub1 <- which(lapply(scaled_ordered_pdb_bfactor, length) %in% quantile(unlist(lapply(scaled_ordered_pdb_bfactor, length)))[1]:quantile(unlist(lapply(scaled_ordered_pdb_bfactor, length)))[2])
bfact_sub1 <- bfact_sub1[-length(bfact_sub1)]
bfact_sub2 <- which(lapply(scaled_ordered_pdb_bfactor, length) %in% quantile(unlist(lapply(scaled_ordered_pdb_bfactor, length)))[2]:quantile(unlist(lapply(scaled_ordered_pdb_bfactor, length)))[3])
bfact_sub2 <- bfact_sub2[-length(bfact_sub2)]
bfact_sub3 <- which(lapply(scaled_ordered_pdb_bfactor, length) %in% quantile(unlist(lapply(scaled_ordered_pdb_bfactor, length)))[3]:quantile(unlist(lapply(scaled_ordered_pdb_bfactor, length)))[4])
bfact_sub3 <- bfact_sub3[-length(bfact_sub3)]
bfact_sub4 <- which(lapply(scaled_ordered_pdb_bfactor, length) %in% quantile(unlist(lapply(scaled_ordered_pdb_bfactor, length)))[4]:quantile(unlist(lapply(scaled_ordered_pdb_bfactor, length)))[5])
bfact_sub4 <- bfact_sub4[-length(bfact_sub4)]
scaled_ordered_pdb_bfactor_subset <- scaled_ordered_pdb_bfactor[bfact_sub1]
res_numbers_scaled_ordered_pdb_bfactor_subset <- lapply(scaled_ordered_pdb_bfactor_subset, FUN=function(x){
  1:length(x)
})
#optional:convert everything to the numbering system for pdz ligands, e.g. 0 is the c-terminal residue, -1, -2, -3... etc.
res_numbers_scaled_ordered_pdb_bfactor_subset_c0 <- lapply(res_numbers_scaled_ordered_pdb_bfactor_subset, FUN=function(x){
  x-x[length(x)]
})
melted_df_bfactor_res <- cbind(reshape2::melt(scaled_ordered_pdb_bfactor_subset), reshape2::melt(res_numbers_scaled_ordered_pdb_bfactor_subset_c0))
melted_df_bfactor_res <- melted_df_bfactor_res[,-2]
colnames(melted_df_bfactor_res) <- c("B_factor", "residue_NtoC", "pdb_id")
melted_df_bfactor_5heb <- melted_df_bfactor_res[which(melted_df_bfactor_res$pdb_id=="5heb"),]
#plot
melted_df_bfactor_res$is_5heb <- rep("no", nrow(melted_df_bfactor_res))
melted_df_bfactor_res$is_5heb[which(melted_df_bfactor_res$pdb_id=="5heb")] <- "yes"
lp <- ggplot(melted_df_bfactor_res, aes(x=residue_NtoC, y=B_factor, shape=pdb_id, color=is_5heb)) +
  geom_line()+
  geom_point(aes(shape=pdb_id))+
  scale_shape_manual(values=0:14)+
  scale_colour_manual(values=c(yes="#DA142A", no="black"))+
  guides(colour="none")+
  xlab("Position of residue (N to C-terminus)")+
  ylab("Normalized B-factor")+
  theme(legend.position="bottom", 
        legend.title = element_blank(), axis.title = element_text(size=12, family="sans"), 
        axis.text = element_text(size=12, family="sans"), legend.text = element_text(size=12, family="sans", ))
bfactor_vs_pos_pdz_ligand <- lp + theme_bw()
#also want to output the bfactor stats for 5heb to use later
write.table(melted_df_bfactor_5heb, "./fuzzy_specificity_data/figures/figure1_genetic_landscape/bfactor_5heb_ligand.txt", sep="\t", quote=F, row.names=F)
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
###
n_c_density <- make_gg_density_nor_fit_same_axis(n_c_combined)
###
#QC -- data quality/summary plots
n_repcorr <- (draw(get_rep_corr_matrix(n_c_combined[which(n_c_combined$data_id=="cript_n"),]), column_title="CRIPT N"))
n_fitcorr <- (draw(get_fitness_corr_matrix(n_c_combined[which(n_c_combined$data_id=="cript_n"),])))
c_repcorr <- (draw(get_rep_corr_matrix(n_c_combined[which(n_c_combined$data_id=="cript_c"),]), column_title="CRIPT C"))
c_fitcorr <- (draw(get_fitness_corr_matrix(n_c_combined[which(n_c_combined$data_id=="cript_c"),])))
cript_c_alive_top_percent <- which(n_c_combined$data_id=="cript_c" & n_c_combined$fitness>quantile(n_c_combined[which(n_c_combined$data_id=="cript_c"),]$fitness, prob=1-1/100))
c_fitcorr_alive <- (draw(get_fitness_corr_matrix(n_c_combined[cript_c_alive_top_percent,]), column_title="CRIPT C top 1% vars"))
#QC -- data quality/summary with density matrix + correlations
n_fitcorr_mat <- ggpairs(n_c_combined[which(n_c_combined$data_id=="cript_n"),grep("fitness[0-9]",colnames(n_c_combined))],
                         lower = list(continuous = "hexbin", combo = "facethist", discrete = "facetbar", na =  "na"))+theme_bw()
c_fitcorr_mat <- ggpairs(n_c_combined[which(n_c_combined$data_id=="cript_c"),grep("fitness[0-9]",colnames(n_c_combined))],
                         lower = list(continuous = "hexbin", combo = "facethist", discrete = "facetbar", na =  "na"))+theme_bw()
c_fitcorr_alive_mat <- ggpairs(n_c_combined[cript_c_alive_top_percent,grep("fitness[0-9]",colnames(n_c_combined))],
                               lower = list(continuous = "hexbin", combo = "facethist", discrete = "facetbar", na =  "na"))+theme_bw()
#plot the AA hamming distance versus fitness for both
n_library <- make_aaham_plot(n_c_combined[which(n_c_combined$data_id=="cript_n"),])
c_library <- make_aaham_plot(n_c_combined[which(n_c_combined$data_id=="cript_c"),])
n_c_combined$Nham_aa <- as.factor(n_c_combined$Nham_aa)
n_c_combined$data_id <- factor(n_c_combined$data_id, levels=c("cript_n", "cript_c"))
n_c_aaham <- ggplot(n_c_combined[which(n_c_combined$Nham_aa %in% c(1,2,3,4)),], aes(x=Nham_aa, y=fitness)) + 
  geom_violin(aes(fill= data_id), draw_quantiles = 0.5) + 
  theme_bw() +
  scale_fill_manual(values=custom_scale_discrete_4, labels=c("N", "C"), breaks=c("cript_n", "cript_c"))+
  xlab("AA Hamming distance")+
  ylab("Fitness")+
  theme(legend.title = element_blank())
#top 1% PWMs based on fitness
ggseq_oneperc_n <- ggseqlogo(top_percent_vars(n_c_combined[which(n_c_combined$data_id=="cript_n"),],1), col_scheme=zappo_col_scheme)+
  scale_x_continuous(breaks=1:4, labels=-7:-4)+
  scale_y_continuous(limits=c(0,3.5))+
  ggtitle(paste("N top 1% (N=", length(top_percent_vars(n_c_combined[which(n_c_combined$data_id=="cript_n"),],1)), ")", sep=""))+
  theme(plot.title = element_text(size = 10, margin=margin(0,0,-10,0)))
ggseq_oneperc_c <- ggseqlogo(top_percent_vars(n_c_combined[which(n_c_combined$data_id=="cript_c"),],1), col_scheme=zappo_col_scheme)+
  scale_x_continuous(breaks=1:4, labels=-3:0)+
  scale_y_continuous(limits=c(0,3.5))+
  ggtitle(paste("C top 1% (N=", length(top_percent_vars(n_c_combined[which(n_c_combined$data_id=="cript_c"),],1)), ")", sep=""))+
  theme(plot.title = element_text(size = 10, margin = margin(0,0,-10,0)))
  
###
pwms_n_c <- ggarrange(ggseq_oneperc_n, ggseq_oneperc_c, nrow=1, ncol=2, legend="none")
###
#plot the number of variants in each hamming distance category that are in top 1% as a fraction of the number of variants in that hamming dist. category
cript_n_subset <- n_c_combined[which(n_c_combined$data_id=="cript_n"),]
cript_c_subset <- n_c_combined[which(n_c_combined$data_id=="cript_c"),]
#take out stop codons
cript_n_subset <- cript_n_subset[-which(cript_n_subset$STOP==T),]
cript_c_subset <- cript_c_subset[-which(cript_c_subset$STOP==T),]
norm_hamm_n <- (table(subset(cript_n_subset$Nham_aa, cript_n_subset$fitness > quantile(cript_n_subset$fitness, prob = 1 - 1/100)))/table(c(cript_n_subset$Nham_aa)))*100
norm_hamm_n <- norm_hamm_n[c("1", "2", "3", "4")]
norm_hamm_c <- (table(subset(cript_c_subset$Nham_aa, cript_c_subset$fitness > quantile(cript_c_subset$fitness, prob = 1 - 1/100)))/table(c(cript_c_subset$Nham_aa)))*100
norm_hamm_c <- norm_hamm_c[c("1", "2", "3", "4")]
fraction_hamm_nc <- rbind(norm_hamm_n, norm_hamm_c)
fraction_hamm_nc_melted <- reshape2::melt(fraction_hamm_nc)
###
aa_ham_top_one_percent <- ggplot(data=fraction_hamm_nc_melted, aes(x=Var2, y=value, fill=Var1))+
  geom_bar(stat="identity", position=position_dodge())+
    scale_fill_manual(name=" ", values=custom_scale_discrete_4, breaks=c("norm_hamm_n", "norm_hamm_c"), labels=c("N", "C"))+
    theme_bw()+
    xlab("AA Hamming distance")+
    ylab("% var top 1% fitness")
###
#plot the number of variants in each hamming distance
no_vars_hamm_n <- table(n_c_combined[which(n_c_combined$data_id=="cript_n"),]$Nham_aa)
no_vars_hamm_n <- no_vars_hamm_n[c("1", "2", "3", "4")]
no_vars_hamm_c <- table(n_c_combined[which(n_c_combined$data_id=="cript_c"),]$Nham_aa)
no_vars_hamm_c <- no_vars_hamm_c[c("1", "2", "3", "4")]
no_vars_bound <- rbind(no_vars_hamm_n, no_vars_hamm_c)
no_vars_bound_melted <- reshape2::melt(no_vars_bound)
###
aa_ham_no_variants <- ggplot(data=no_vars_bound_melted, aes(x=Var2, y=value, fill=Var1))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(name=" ", values=custom_scale_discrete_4, breaks=c("no_vars_hamm_n", "no_vars_hamm_c"), labels=c("N", "C"))+
  theme_bw()+
  xlab("AA Hamming distance")+
  ylab("Number of variants")
###
#looking at the actual top sequences delineated by hamming distance
cript_n_top_one_perc_subset <- data.frame(aa_seq=subset(cript_n_subset$aa_seq, cript_n_subset$fitness > quantile(cript_n_subset$fitness, prob = 1 - 1/100)), nham_aa = subset(cript_n_subset$Nham_aa, cript_n_subset$fitness > quantile(cript_n_subset$fitness, prob = 1 - 1/100)))
cript_n_top_one_perc_subset_split <- split.data.frame(cript_n_top_one_perc_subset, f=cript_n_top_one_perc_subset$nham_aa)
cript_n_top_one_perc_subset_split <- cript_n_top_one_perc_subset_split[-which(names(cript_n_top_one_perc_subset_split)=="0")]
cript_n_top_one_perc_subset_split_seqs <- lapply(cript_n_top_one_perc_subset_split, FUN=function(x){x$aa_seq})
names(cript_n_top_one_perc_subset_split_seqs) <- paste("n = ", as.vector(unlist(lapply(cript_n_top_one_perc_subset_split_seqs, length))), sep=" ")

cript_c_top_one_perc_subset <- data.frame(aa_seq=subset(cript_c_subset$aa_seq, cript_c_subset$fitness > quantile(cript_c_subset$fitness, prob = 1 - 1/100)), nham_aa = subset(cript_c_subset$Nham_aa, cript_c_subset$fitness > quantile(cript_c_subset$fitness, prob = 1 - 1/100)))
cript_c_top_one_perc_subset_split <- split.data.frame(cript_c_top_one_perc_subset, f=cript_c_top_one_perc_subset$nham_aa)
cript_c_top_one_perc_subset_split <- cript_c_top_one_perc_subset_split[-which(names(cript_c_top_one_perc_subset_split)=="0")]
cript_c_top_one_perc_subset_split_seqs <- lapply(cript_c_top_one_perc_subset_split, FUN=function(x){x$aa_seq})
names(cript_c_top_one_perc_subset_split_seqs) <- paste("n = ", as.vector(unlist(lapply(cript_c_top_one_perc_subset_split_seqs, length))), sep=" ")
#put them together
cript_all_logos_aaham <- c(cript_n_top_one_perc_subset_split_seqs[1], cript_c_top_one_perc_subset_split_seqs[1],
                           cript_n_top_one_perc_subset_split_seqs[2], cript_c_top_one_perc_subset_split_seqs[2],
                           cript_n_top_one_perc_subset_split_seqs[3], cript_c_top_one_perc_subset_split_seqs[3],
                           cript_n_top_one_perc_subset_split_seqs[4], cript_c_top_one_perc_subset_split_seqs[4])
#make pwms
seq_logos_cript_oneperc <- ggseqlogo(cript_all_logos_aaham, nrow=1, col_scheme=zappo_col_scheme)+
  theme(legend.position="top")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#show the number of seqs in each seqlogo
#more seq pwm analysis
num_pwm_c <- lapply(cript_c_top_one_perc_subset_split, nrow)
names(num_pwm_c) <- c("1", "2", "3", "4")
num_pwm_c <- do.call(cbind, num_pwm_c)
num_pwm_n <- lapply(cript_n_top_one_perc_subset_split, nrow)
names(num_pwm_n) <- c("1", "2", "3", "4")
num_pwm_n <- do.call(cbind, num_pwm_n)
num_pwm_bound <- rbind(num_pwm_n, num_pwm_c)
rownames(num_pwm_bound) <- c("num_pwm_n", "num_pwm_c")
num_pwm_bound_melted <- reshape2::melt(num_pwm_bound)
###
num_variants_matching_pwm <- ggplot(data=num_pwm_bound_melted, aes(x=Var2, y=value, fill=Var1))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(name=" ", values=custom_scale_discrete_4, breaks=c("num_pwm_n", "num_pwm_c"), labels=c("N", "C"))+
  theme_bw()+
  xlab("AA Hamming distance")+
  ylab("Number of variants")+
  scale_y_continuous(breaks=seq(0,500,200))

###
############amino acid properties analysis
###from the figure2_combo_n_c file 
n_no_stop <- cript_n_subset 
c_no_stop <- cript_c_subset
#get fitness
n_fit <- n_no_stop$fitness
c_fit <- c_no_stop$fitness
#aa files
aa_properties_file <- "./fuzzy_specificity_data/aa_properties/amino_acid_properties_annotated_supplementary.txt"
selected_identifiers_file <- "./fuzzy_specificity_data/aa_properties/selected.amino_acid_properties.txt"
#Load amino acid properties
aa_df <- read.table(file=aa_properties_file, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
#Subset to selected amino acid properties if supplied
selected_identifiers <- readLines(selected_identifiers_file)
if(!is.null(selected_identifiers)){
  aa_df <- aa_df[rownames(aa_df) %in% selected_identifiers | grepl("BBSUPPL", rownames(aa_df)),]
}
#Evidences list
aa_evidences <- as.list(aa_df[,1])
names(aa_evidences) <- rownames(aa_df)
#Reformat and add AA identities
aa_df <- aa_df[,-1]
#Remove properties with NAs
aa_evidences <- aa_evidences[apply(is.na(aa_df), 1, sum)==0]
aa_df <- aa_df[apply(is.na(aa_df), 1, sum)==0,]
#Normalise
aa_mat <- scale(t(aa_df))
#function to get score for each aa property and sum it over the given aa_string
get_score <- function(aa_string){
  split_aa_string <- str_split(aa_string, "")
  list_of_values <- sapply(split_aa_string[[1]], FUN=function(x){
    which(rownames(aa_mat)==x)
  })
  return(apply(aa_mat[c(list_of_values),], 2, sum))
}
n_scores <- sapply(n_no_stop$aa_seq, FUN=get_score, USE.NAMES=F)
c_scores <- sapply(c_no_stop$aa_seq, FUN=get_score, USE.NAMES = F)
#look at correlations even though it doesn't quite make sense with discrete distributions..
get_corrs <- function(indep, dep, corr_type){
  corr_obj <- cor.test(x=indep, y=dep, method=corr_type)
  return(corr_obj$estimate)
}
#do the correlations for n
corrs_n_pearson <- NULL;
for (i in 1:nrow(n_scores)){
  corrs_n_pearson[i] <- get_corrs(n_scores[i,], n_fit, "pearson")
}
corrs_n_spearman <- NULL;
for (i in 1:nrow(n_scores)){
  corrs_n_spearman[i] <- get_corrs(n_scores[i,], n_fit, "spearman")
}
#do the correlations for c
corrs_c_pearson <- NULL;
for (i in 1:nrow(c_scores)){
  corrs_c_pearson[i] <- get_corrs(c_scores[i,], c_fit, "pearson")
}
corrs_c_spearman <- NULL;
for (i in 1:nrow(c_scores)){
  corrs_c_spearman[i] <- get_corrs(c_scores[i,], c_fit, "spearman")
}
corr_df <- cbind(aa_evidences,corrs_n_pearson, corrs_n_spearman, corrs_c_pearson, corrs_c_spearman)
names_of_feats <- paste(rownames(corr_df), aa_evidences, sep=": ")
#better plot as vioplot for eg the discrete variables
corr_threshold <- 0.4
feats_with_high_cor_c <- union(which(abs(corrs_c_pearson)>corr_threshold), which(abs(corrs_c_spearman)>corr_threshold))
feats_with_high_cor_n <- union(which(abs(corrs_n_spearman)>corr_threshold), which(abs(corrs_n_spearman)>corr_threshold))

feats_n <- names_of_feats[feats_with_high_cor_n]
feats_c <- names_of_feats[feats_with_high_cor_c]

#present the top correlations
names_top_corrs_n_c <- c(feats_n, feats_c)
values_top_corrs <- corr_df[c(feats_with_high_cor_n,feats_with_high_cor_c),]
df_to_plot <- rbind(as.vector(data.frame(values_top_corrs)$corrs_n_spearman), as.vector(data.frame(values_top_corrs)$corrs_c_spearman))
rownames(df_to_plot) <- c("n_spearman", "c_spearman")
colnames(df_to_plot) <- names_top_corrs_n_c
df_to_plot <- data.frame(df_to_plot)
df_v2_to_plot <- apply(df_to_plot, 1, as.numeric)
rownames(df_v2_to_plot) <- names_top_corrs_n_c
melted_df_v2_to_plot <- reshape2::melt(df_v2_to_plot)
melted_df_v2_to_plot$Var1 <- factor(melted_df_v2_to_plot$Var1, levels=unique(melted_df_v2_to_plot$Var1[order(melted_df_v2_to_plot$value, decreasing=T)]))
top_corrs_plot <- ggplot(melted_df_v2_to_plot, aes(x=Var1, y=value, fill=Var2))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(name=" ", values=custom_scale_discrete_4, breaks=c("n_spearman", "c_spearman"), labels=c("N", "C"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("AA Property")+
  ylab("Spearman's rho")

###
#net charge important in N but not C -- show in detail
#function to calculate net charge
calc_net_charge <- function(x){
  pos_charge <- str_count(x, "[RK]")
  neg_charge <- str_count(x, "[DE]")
  net_charge <- pos_charge - neg_charge
  return(net_charge)
}
net_charge_scores_simple_n <- unlist(lapply(n_no_stop$aa_seq, calc_net_charge))
net_charge_scores_simple_c <- unlist(lapply(c_no_stop$aa_seq, calc_net_charge))
net_charge_n_df <- data.frame(score=net_charge_scores_simple_n, fitness=n_fit, data_id=rep("cript_n", length(net_charge_scores_simple_n)))
net_charge_c_df <- data.frame(score=net_charge_scores_simple_c, fitness=c_fit, data_id=rep("cript_c", length(net_charge_scores_simple_c)))
net_charge_nc_bound <- rbind(net_charge_n_df, net_charge_c_df)
net_charge_nc_bound$score <- as.factor(net_charge_nc_bound$score)
net_charge_nc_bound$data_id <- as.factor(net_charge_nc_bound$data_id)
net_charge_vioplots <- ggplot(data=net_charge_nc_bound, aes(x=score, y=fitness))+
  geom_violin(aes(fill=data_id), draw_quantiles = 0.5)+
  theme_bw()+
  scale_fill_manual(values=custom_scale_discrete_4, labels=c("N", "C"), breaks=c("cript_n", "cript_c"))+
  xlab("Net charge")+
  ylab("Fitness")+
  theme(legend.title = element_blank())



