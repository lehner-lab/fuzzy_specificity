################################################################################################
#Source code for Figure7 (+ related supplementary material) from 
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
#summarize in the same visualization the specificity changes between cript cis interactions and trans interactions with pdz
cript_pdz_trans_intxns <- read.delim("./fuzzy_specificity_data/figures/figure7_summary_intxns/num_specificity_muts_per_pair_pdz_cript.txt")
cript_cript_cis_intxns <- read.delim("./fuzzy_specificity_data/figures/figure7_summary_intxns/num_specificity_muts_per_pair_cript_cript.txt")
#read in distances for trans intxns
cript_pdz_distances <- read.delim("./fuzzy_specificity_data/figures/figure7_summary_intxns/df_dist_fdr_spec_change.txt")
#translate coordinates to pdb id coordinates
cript_pdz_trans_intxns$position_a_pdb <- cript_pdz_trans_intxns$position_a
cript_pdz_trans_intxns$position_b_pdb <- cript_pdz_trans_intxns$position_b+9
cript_cript_cis_intxns$position_a_pdb <- cript_cript_cis_intxns$position_a+9
cript_cript_cis_intxns$position_b_pdb <- cript_cript_cis_intxns$position_b+9
#add distance column in case you want to filter later
if(identical(cript_pdz_trans_intxns$pair_name, cript_pdz_distances$pair_name)==T){
  cript_pdz_trans_intxns$distance_pdb <- cript_pdz_distances$dist_to_lig
}else{
  print("dfs do not match")
}
cript_cript_cis_intxns$distance_pdb <- 10 
#write the template file
#make all the pseudobonds from the alpha carbon [ca]
cript_pdz_trans_intxns$pseudobond_command <- paste(paste("/A:", cript_pdz_trans_intxns$position_a, "@ca", sep=""),
                                                   paste("/B:", cript_pdz_trans_intxns$position_b_pdb, "@ca", sep=""), sep=" ")

cript_cript_cis_intxns$pseudobond_command <- paste(paste("/B:", cript_cript_cis_intxns$position_a_pdb, "@ca", sep=""),
                                                   paste("/B:", cript_cript_cis_intxns$position_b_pdb, "@ca", sep=""), sep=" ")
#merge the two dfs
cript_cis_trans_intxns <- rbind(cript_pdz_trans_intxns, cript_cript_cis_intxns)
#separate based on number of hits
#20+
cript_cis_trans_intxns_20 <- cript_cis_trans_intxns[which(cript_cis_trans_intxns$hit_in_samp>=20),]
cript_cis_trans_intxns_10 <- cript_cis_trans_intxns[which(cript_cis_trans_intxns$hit_in_samp>=10 & cript_cis_trans_intxns$hit_in_samp < 20),]
cript_cis_trans_intxns_5 <- cript_cis_trans_intxns[which(cript_cis_trans_intxns$hit_in_samp>=5 & cript_cis_trans_intxns$hit_in_samp < 10),]
#cript_cis_trans_intxns_1 <- cript_cis_trans_intxns[which(cript_cis_trans_intxns$hit_in_samp>=1 & cript_cis_trans_intxns$hit_in_samp < 5),]
cript_cis_trans_intxns_1 <- cript_cis_trans_intxns[which(cript_cis_trans_intxns$hit_in_samp>=2 & cript_cis_trans_intxns$hit_in_samp < 5),]
spacer_var <- "; space to specify colors etc"
#write them into a file
writeLines(c(spacer_var,cript_cis_trans_intxns_20$pseudobond_command, spacer_var, cript_cis_trans_intxns_10$pseudobond_command, spacer_var, cript_cis_trans_intxns_5$pseudobond_command, spacer_var, cript_cis_trans_intxns_1$pseudobond_command), 
           con = "./fuzzy_specificity_data/figures/figure7_summary_intxns/cis_trans_pseudobonds_more_than_two.pb", sep="\n") 

#separate based on number of hits + allosteric vs. not
#20+
cript_cis_trans_intxns_20_direct <- cript_cis_trans_intxns[which(cript_cis_trans_intxns$hit_in_samp>=20 & cript_cis_trans_intxns$distance_pdb<5),]
cript_cis_trans_intxns_20_allosteric <- cript_cis_trans_intxns[which(cript_cis_trans_intxns$hit_in_samp>=20 & cript_cis_trans_intxns$distance_pdb>=5),]
cript_cis_trans_intxns_10_direct <- cript_cis_trans_intxns[which(cript_cis_trans_intxns$hit_in_samp>=10 & cript_cis_trans_intxns$hit_in_samp < 20 & cript_cis_trans_intxns$distance_pdb<5),]
cript_cis_trans_intxns_10_allosteric <- cript_cis_trans_intxns[which(cript_cis_trans_intxns$hit_in_samp>=10 & cript_cis_trans_intxns$hit_in_samp < 20 & cript_cis_trans_intxns$distance_pdb>=5),]
cript_cis_trans_intxns_5_direct <- cript_cis_trans_intxns[which(cript_cis_trans_intxns$hit_in_samp>=5 & cript_cis_trans_intxns$hit_in_samp < 10 & cript_cis_trans_intxns$distance_pdb<5),]
cript_cis_trans_intxns_5_allosteric <- cript_cis_trans_intxns[which(cript_cis_trans_intxns$hit_in_samp>=5 & cript_cis_trans_intxns$hit_in_samp < 10 & cript_cis_trans_intxns$distance_pdb>=5),]
spacer_var <- "; space to specify colors etc"
#write them into a file
writeLines(c(spacer_var,cript_cis_trans_intxns_20_direct$pseudobond_command, 
             spacer_var, cript_cis_trans_intxns_20_allosteric$pseudobond_command,
             spacer_var, cript_cis_trans_intxns_10_direct$pseudobond_command, 
             spacer_var, cript_cis_trans_intxns_10_allosteric$pseudobond_command,
             spacer_var, cript_cis_trans_intxns_5_direct$pseudobond_command, 
             spacer_var, cript_cis_trans_intxns_5_allosteric$pseudobond_command), 
           con = "./fuzzy_specificity_data/figures/figure7_summary_intxns/cis_trans_pseudobonds_more_than_five_plus_dist.pb", sep="\n") 

