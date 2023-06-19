#!/usr/bin/r

################################################################################
# processingHUMAnN3bugfile.R                                                   #
# Name: Suyeon Kim                                                             #
# Date: MAR-10-2023                                                            #
# Last updated: Jun-13-2023                                                    #
# Goal: Read processed taxonomic/bug file and combine them in dataframe format #
#       (row: taxonomy, column: sample IDs)                                    #                            
################################################################################

mergingBugfiles<- function(files){
  # While we read an each MGX abundance input file, we want to append output of 
  # Taxon node dataframe using 'append=TRUE' parameter in write.table function. 
  for(i in 1:length(files)){
    # Let's read Metegenomic abundance input data in the list of files.  
    print(paste("This is file #", i))
    mgx_dat<- read.table(files[i], sep="\t",header=F, check.names=F)
    
    # grep "taxonomy" and "Abundance" column from mgx_dat file. 
    interest_col <- mgx_dat[,c(1,3)]
    colName <- gsub(".*/", "", files[[i]])
    colName <- gsub("\\_metaphlan.*", "",colName)
    
    colnames(interest_col) <- c("Taxonomy", colName)
    
    if (i == 1){
      merged_dat <- interest_col
    }
    else{
      merged_dat <- merge(merged_dat, interest_col, by = "Taxonomy", all = T)
    }
  }
  
  # Check the merged bug list abundance data                                     
  rownames(merged_dat)<- merged_dat$Taxonomy
  bug_dat<-merged_dat[,-1]
  
  # Filtering the previously merged taxonomic profiles                           
  # Remove 'Archaea' profile 
  filt_bug<- bug_dat[!grepl("k__Archaea", rownames(bug_dat)),]
  filt_bug<-filt_bug[!grepl("k__Eukaryota", rownames(filt_bug)),]
  # Replace NA value with '0' value 
  filt_bug[is.na(filt_bug)]<-0

  return(filt_bug)
}