#!/usr/bin/r

################################################################################
# Script: getKoBySpecies.R                                                     #
# Name: Suyeon Kim                                                             #
# Created: May-16-2023                                                         #
# Last edited: Jun-16-2023                                                     #
################################################################################
library(igraph)
library(MicrobiomeProfiler)

# Note that parameter 'bacList' for this function is needed to be the 'character'
getKoListForSpeciesOfInterest<- function(koDf, groupname, bacList){
  # Fisrt, get ko-IDs and species information for each group 
  groupKO<- koDf[koDf$GROUP == groupname,]
  # Since each species have multiple KO-IDs, let's aggregate KO-IDs by the list 
  # of species in dataframe 'groupKO'. 
  summaryTbl<- aggregate(groupKO$Koid, 
                         by=list(groupKO$Species), function(x) unique(x))
  colnames(summaryTbl)<- c("species","ko")
  # Then, let's get the information, containing our species of interest 
  # and their list of KO-IDs. 
  koList<- summaryTbl[match(bacList, summaryTbl$species),]
  
  return(koList)
} 

