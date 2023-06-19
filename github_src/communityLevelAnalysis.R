#!/usr/bin/r

################################################################################
# Script: communityLevelAnalysis.R                                             #
# Name: Suyeon Kim                                                             #
# Created: May-16-2023                                                         #
# Last edited: Jun-16-2023                                                     #
################################################################################
library(igraph)
library(MicrobiomeProfiler)

args = commandArgs(trailingOnly=TRUE)
print(args)

sourceDir <- args[1]
InFile1<- args[2]
InFile2<- args[3]
InFile3<- args[4]
InFile4<- args[5]
InFile5<- args[6]
InFile6<- args[7]
InFile7<- args[8]
InFile8<- args[9]

################################################################################
# List of functions                                                            #
################################################################################
## In this function, we perform 'Leiden community detection algorithm' to detect
## species communities in the network 
performLeidenAlgorithm<- function(graph){
  set.seed(1)
  r <- quantile(strength(graph))[2] / (gorder(graph) - 1)
  resModularity<- cluster_leiden(
    graph,
    objective_function = c("CPM"),
    weights = NULL,
    resolution_parameter = r,
    beta = 0.01,
    initial_membership = NULL,
    n_iterations = 2,
    vertex_weights = NULL
  )
  return(resModularity)
}

listSpeciesByCommunity<- function(inputGraph, resComm){
  # Get node label from input graph 
  nodeNames<- V(inputGraph)$name
  # Get the total number of clusters 
  totalCluster<- resComm$nb_clusters 
  
  # Get list of species for each detected community, 
  listSpecies<- list()
  for( i in 1:totalCluster){
    numOfcommunity<- nodeNames[resComm$membership == i]
    listSpecies[[i]]<- numOfcommunity
  }
  
  return(listSpecies)
}

################################################################################
# Read the Input file into R.                                                  #
################################################################################
# This should result in a dataframe called 'ctrl 1' that contains 522 species
# co-occurrences in control group from IBD 1 dataset. 
ctrl1<- read.csv(InFile1, row.names=1)
# This should result in a dataframe called 'ctrl 2' that contains 2743 species
# co-occurrences in control group from IBD 2 dataset. 
ctrl2<- read.csv(InFile2, row.names=1)
# This should result in a dataframe called 'cd1' that contains 721 species
# co-occurrences in control group from IBD 1 dataset. 
cd1<- read.csv(InFile3, row.names=1)
# This should result in a dataframe called 'cd 2' that contains 3722 species
# co-occurrences in control group from IBD 2 dataset. 
cd2<- read.csv(InFile4, row.names=1)
# This should result in a dataframe called 'uc 1' that contains 329 species
# co-occurrences in control group from IBD 1 dataset. 
uc1<- read.csv(InFile5, row.names=1)
# This should result in a dataframe called 'uc 1' that contains 3987 species
# co-occurrences in control group from IBD 1 dataset. 
uc2<- read.csv(InFile6, row.names=1)

# This should result in a dataframe called 'ko1' that contains ko-IDs, genus and
# species name, sample-IDs, group information. 
ko1<- read.csv(InFile7, sep=",", header=T)
# This should result in a dataframe called 'ko1' that contains ko-IDs, genus and
# species name, sample-IDs, group information. 
ko2<- read.csv(InFile8, sep=",", header=T)

################################################################################
# Main Program                                                                 #
################################################################################
source(paste0(sourceDir, "/", "createIgraphObject.R"))
source(paste0(sourceDir, "/", "getKoBySpecies.R"))

## Create igrpah graphs objects from dataframes         
ctrl1graph<- creatingUNWEIGHTED_igraph_Graphs(ctrl1, tx1="sp1", tx2="sp2")
ctrl2graph<- creatingUNWEIGHTED_igraph_Graphs(ctrl2, tx1="sp1", tx2="sp2")

cd1graph<- creatingUNWEIGHTED_igraph_Graphs(cd1, tx1="sp1", tx2="sp2")
cd2graph<- creatingUNWEIGHTED_igraph_Graphs(cd2, tx1="sp1", tx2="sp2")

uc1graph<- creatingUNWEIGHTED_igraph_Graphs(uc1, tx1="sp1", tx2="sp2")
uc2graph<- creatingUNWEIGHTED_igraph_Graphs(uc2, tx1="sp1", tx2="sp2")

## Perform the 'community detection algorithm' for each group
ctrl1Comm<- performLeidenAlgorithm(ctrl1graph)
ctrl2Comm<- performLeidenAlgorithm(ctrl2graph)

cd1Comm<- performLeidenAlgorithm(cd1graph)
cd2Comm<- performLeidenAlgorithm(cd2graph)

uc1Comm<- performLeidenAlgorithm(uc1graph)
uc2Comm<- performLeidenAlgorithm(uc2graph)

## Get list of species for each community in each phenotype
speciesCommCtrl1<- listSpeciesByCommunity(ctrl1graph, ctrl1Comm)
speciesCommCtrl2<- listSpeciesByCommunity(ctrl2graph, ctrl2Comm)

speciesCommCD1<- listSpeciesByCommunity(cd1graph, cd1Comm)
speciesCommCD2<- listSpeciesByCommunity(cd2graph, cd2Comm)

speciesCommUC1<- listSpeciesByCommunity(uc1graph, uc1Comm)
speciesCommUC2<- listSpeciesByCommunity(uc2graph, uc2Comm)

## Get list of KO-IDs based on our species of interest from the identified 
## community and perform the pathway enrichement analysis.   
ctrl1Enrich <- lapply(
  lapply(speciesCommCtrl1, function(x) {
    getKoListForSpeciesOfInterest(ko1, groupname="Control", x)
  }), function(y) {
    data.frame(enrichKO(unique(unlist(y$ko)), universe = ko1$Koid))
  })

ctrl2Enrich <- lapply(
  lapply(speciesCommCtrl2, function(x) {
    getKoListForSpeciesOfInterest(ko2, groupname="Control", x)
  }), function(y) {
    data.frame(enrichKO(unique(unlist(y$ko)), universe = ko2$Koid))
  })

cd1Enrich <- lapply(
  lapply(speciesCommCD1, function(x) {
    getKoListForSpeciesOfInterest(ko1, groupname="CD", x)
  }), function(y) {
    data.frame(enrichKO(unique(unlist(y$ko)), universe = ko1$Koid))
  })

cd2Enrich <- lapply(
  lapply(speciesCommCD2, function(x) {
    getKoListForSpeciesOfInterest(ko2, groupname="CD", x)
  }), function(y) {
    data.frame(enrichKO(unique(unlist(y$ko)), universe = ko2$Koid))
  })

uc1Enrich <- lapply(
  lapply(speciesCommUC1, function(x) {
    getKoListForSpeciesOfInterest(ko1, groupname="UC", x)
  }), function(y) {
    data.frame(enrichKO(unique(unlist(y$ko)), universe = ko1$Koid))
  })

uc2Enrich <- lapply(
  lapply(speciesCommUC2, function(x) {
    getKoListForSpeciesOfInterest(ko2, groupname="UC", x)
  }), function(y) {
    data.frame(enrichKO(unique(unlist(y$ko)), universe = ko2$Koid))
  })

## Comparative Analysis: 1) Compare the corresponding lists of enriched pathways 
##                          between disease phenotype and healthy phenotype 
##                       2) Identify the pathways unique to 'CD' group that are 
##                          compared across two datasets. 

cd1minusctrl1 <- setdiff(unique(unlist(lapply(cd1Enrich,function(x) 
                          if(nrow(x)>0) x[[2]]))), 
                         unique(unlist(lapply(ctrl1Enrich,function(x) 
                          if(nrow(x)>0) x[[2]]))))

cd2minusctrl2 <- setdiff(unique(unlist(lapply(cd2Enrich,function(x) 
                          if(nrow(x)>0)  x[[2]]))), 
                         unique(unlist(lapply(ctrl2Enrich,function(x)
                          if(nrow(x)>0) x[[2]]))))

## print the results of CD group 
cat(paste0("Enriched Pathways in CD group: ", 
           intersect(cd1minusctrl1, cd2minusctrl2), "\n"))

## Comparative Analysis: 1) Compare the corresponding lists of enriched pathways 
##                          between disease phenotype and healthy phenotype 
##                       2) Identify the pathways unique to 'UC' group that are 
##                          compared across two datasets. 


uc1minusctrl1 <- setdiff(unique(unlist(lapply(uc1Enrich,function(x) 
                          if(nrow(x)>0) x[[2]]))), 
                         unique(unlist(lapply(ctrl1Enrich,function(x) 
                          if(nrow(x)>0) x[[2]]))))

uc2minusctrl2 <- setdiff(unique(unlist(lapply(uc2Enrich,function(x) 
                          if(nrow(x)>0)  x[[2]]))), 
                         unique(unlist(lapply(ctrl2Enrich,function(x)
                          if(nrow(x)>0) x[[2]]))))

## print the results of UC group 
cat(paste0("Enriched Pathways in UC group: ", 
           intersect(uc1minusctrl1, uc2minusctrl2), "\n"))


