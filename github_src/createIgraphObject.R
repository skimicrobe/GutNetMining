#!/usr/bin/r

################################################################################
# createIgraphObject.R                                                         #
# Name: Suyeon Kim                                                             #
# Date: Feb-20-2023                                                            #
# Last updated: Feb-20-2023                                                    #
# Goal: read co-occurrence network in igraph object                            #
################################################################################
library(stringr)
library(igraph)

# This function will split taxonomy information and create new columns for each
# taxon level. 
preProcessing_coOccur_mat<- function(cocur_mat){
  # Add taxonomy info onto the 'cocur' dataframe to split in the next step. 
  cocur_mat$sp1_taxon<- cocur_mat$sp1
  cocur_mat$sp2_taxon<- cocur_mat$sp2
  
  # Split 'sp1' and 'sp2' columns into number of taxons in taxonomy lineage.
  split_sp1_tx<- data.frame(str_split_fixed(cocur_mat$sp1_taxon, "\\|", n=7))
  colnames(split_sp1_tx)<- c("sp1_k","sp1_p","sp1_c","sp1_o","sp1_f", "sp1_g",
                             "sp1_s")
  split_sp2_tx<- data.frame(str_split_fixed(cocur_mat$sp2_taxon, "\\|", n=7))
  colnames(split_sp2_tx)<-  c("sp2_k","sp2_p","sp2_c","sp2_o","sp2_f", "sp2_g",
                              "sp2_s")
  
  # Combine previously created splitted taxon info onto the 'original dataframe'
  analysis_tbl<- cbind(split_sp1_tx, split_sp2_tx, cocur_mat)
  
  return(analysis_tbl)
}

# This function creates an igraph graph from co-occurrence dataframes containing 
# the edge list and edge/verte attributes. 
# 'edge_col' parameter was added.
creating_igraph_Graphs<- function(cooccur_mat, tx1, tx2, edge_col){
  tx_level_tbl<- data.frame(cooccur_mat[,c(tx1, tx2, edge_col)])
  colnames(tx_level_tbl)[3]<- "weight"
  gp_taxon<- graph_from_data_frame(tx_level_tbl, directed=FALSE, vertices=NULL)
  
  return(gp_taxon)
}

creatingUNWEIGHTED_igraph_Graphs<- function(cooccur_mat, tx1, tx2){
  tx_level_tbl<- data.frame(cooccur_mat[,c(tx1, tx2)])
  #colnames(tx_level_tbl)[3]<- "weight"
  gp_taxon<- graph_from_data_frame(tx_level_tbl, directed=FALSE, vertices=NULL)
  
  return(gp_taxon)
}
