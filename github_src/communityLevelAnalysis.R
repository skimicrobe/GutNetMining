#!/usr/bin/r

################################################################################
# Script: communityLevelAnalysis.R                                             #
# Name: Suyeon Kim                                                             #
# Created: May-16-2023                                                         #
# Last edited: Jan-30-2024                                                     #
################################################################################
library(igraph)
library(MicrobiomeProfiler)
library(pathview)

args = commandArgs(trailingOnly=TRUE)
print(args)

sourceDir<- args[1]
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
  totalCluster<- length(unique(resComm$membership)) 
  
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
print("Opening the Input files...")
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
print("Performing the pathway enrichment analysis...")
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

allEnrich <- list("CD1"= cd1Enrich, "CD2"=cd2Enrich, "UC1"=uc1Enrich, "UC2" = uc2Enrich)


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
cat(paste0("Enriched Pathways in CD groups: ", 
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
cat(paste0("Enriched Pathways in UC groups: ", 
           intersect(uc1minusctrl1, uc2minusctrl2), "\n"))

################################################################################
# For table4 and table5                                                        #                                   
################################################################################
# In this function, create a list of co-occurring speies from the co-occurrence
# network 
listOfcooccurSpecies<- function(cooccurNetw){
  species1<- cooccurNetw[,1]
  class(species1)
  species2<- cooccurNetw[,2]
  class(species2)
  
  listOfSpecies<- unique(union(species1, species2))
  
  return(listOfSpecies)
}

# In this function, get the number of co-occurrence species and their name for 
# each cluster belonging to the enriched pathway. 
findEnrichmentTable <- function(pathwayToLocate, enrichList, KoLongDf, clusterList, group) {
  res <- data.frame()
  for(cluster in 1:length(enrichList)) {
    print(cluster)
    ibdMapList = enrichList[[cluster]]
    
    if(nrow(ibdMapList)>0){
      clusterMembers = clusterList[[cluster]]
      print(clusterMembers)
      
      if(length(which(ibdMapList$Description == pathwayToLocate)) > 0) {
        print(ibdMapList$Description)
        geneList<- unlist(strsplit(ibdMapList[which(ibdMapList$Description == pathwayToLocate),"geneID"],"\\/"))
        print(geneList)
        matchedGene<- KoLongDf[which(KoLongDf$Koid %in% geneList & 
                                       KoLongDf$GROUP == group & 
                                       KoLongDf$cooccur == TRUE),]
        print(dim(matchedGene))
        nameOfsp<- intersect(unique(matchedGene$Species), as.character(clusterMembers))
        
        row <- cbind(group, pathwayToLocate, cluster, "specieslist" = list(nameOfsp), "count" = length(nameOfsp))
        res <- rbind(res,row)
      }
    }
  }
  return(res)
}


findKOInfo <- function(pathway, soi, enrichTable, KoLongDf, group) {
  geneList <- unlist(strsplit(enrichTable[enrichTable$Description == pathway, "geneID"],"\\/"))
  matchedGene<- KoLongDf[which(KoLongDf$Koid %in% geneList & 
                                 KoLongDf$cooccur == TRUE),]
  speciesGenes <- matchedGene[matchedGene$Species %in% soi,]
  
  speciesKO <- aggregate(speciesGenes$Koid,by=list(speciesGenes$Species),function(x) list(unique(x))) 
  colnames(speciesKO) <- c("species","KOlist")
  speciesKO$count <- unlist(lapply(speciesKO$KOlist,length))
  return(speciesKO)
}

print("Creating the result for Table 4...")
listSpeciesCD1<- listOfcooccurSpecies(cd1)
listSpeciesCD2<- listOfcooccurSpecies(cd2)

listSpeciesUC1<- listOfcooccurSpecies(uc1)
listSpeciesUC2<- listOfcooccurSpecies(uc2)

cd1Group<- ko1[ko1$GROUP == "CD",]
cd1Group$cooccur<- cd1Group$Species %in% listSpeciesCD1

cd2Group<- ko2[ko2$GROUP == "CD",]
cd2Group$cooccur<- cd2Group$Species %in% listSpeciesCD2

uc1Group<- ko1[ko1$GROUP == "UC",]
uc1Group$cooccur<- uc1Group$Species %in% listSpeciesUC1

uc2Group<- ko2[ko2$GROUP == "UC",]
uc2Group$cooccur<- uc2Group$Species %in% listSpeciesUC2

allGroups <- rbind(cbind("groupName"="CD1",cd1Group), cbind("groupName"="CD2",cd2Group), 
                   cbind("groupName"="UC1",uc1Group), cbind("groupName"="UC2",uc2Group))

## Find common enriched pathways across two IBD datasets                        
listPathwaysCD<- intersect(cd1minusctrl1, cd2minusctrl2)
listPathwaysUC<- intersect(uc1minusctrl1, uc2minusctrl2)

## Find the number of co-occurring species on these common enriched pathways                                                                                   
cdClusterRes <- list()
ucClusterRes <- list()

for(pathway in listPathwaysCD){
  cdClusterRes[[pathway]][["CD1"]] <- findEnrichmentTable(pathway, cd1Enrich, 
                                                          cd1Group, 
                                                          cd1Comm, group="CD")
  cdClusterRes[[pathway]][["CD2"]] <- findEnrichmentTable(pathway, cd2Enrich, 
                                                          cd2Group, 
                                                          cd2Comm, group="CD")
}

for(pathway in listPathwaysUC){
  ucClusterRes[[pathway]][["UC1"]] <- findEnrichmentTable(pathway, uc1Enrich, 
                                                          uc1Group, 
                                                          uc1Comm, group="UC")
  ucClusterRes[[pathway]][["UC2"]] <- findEnrichmentTable(pathway, uc2Enrich, 
                                                          uc2Group, 
                                                          uc2Comm, group="UC")
}
keys <- unique(c(names(cdClusterRes), names(ucClusterRes)))
allClusterRes <- setNames(mapply(c, cdClusterRes[keys], ucClusterRes[keys]), keys)

printTable <- function(outputClusterRes) {
  # saving output table
  table4 <- data.frame()
  for(i in 1:length(outputClusterRes)){
    pathwayTable <- outputClusterRes[[i]]
    for(j in 1:length(pathwayTable)){
      pathwayGroupName <- names(pathwayTable)[j]
      ptable <- pathwayTable[[j]]
      out <- cbind("groupName"=pathwayGroupName,
                   "pathway" = as.character(ptable$pathwayToLocate), 
                   "cluster"= as.character(ptable$cluster), 
                   "count" = as.character(ptable$count))
      table4 <- rbind(table4, out)
    }
  }
  return(table4)
}
# for both cd and uc
finalTable4 <- rbind(printTable(cdClusterRes), printTable(ucClusterRes))
            
print("Saving the Table 4 output...")         
write.csv(finalTable4, file="outputTable4.csv",row.names = F)

## further dig into each pathway
print("Creating the output for table 5...")
table5 <- data.frame()
pathwayOfInterest <- unique(finalTable4$pathway)
groupNames <- unique(finalTable4$groupName)
for(pathway in pathwayOfInterest){
  #pathway, soi, enrichTable, KoLongDf, group
  print(pathway)
  for(group in groupNames) {
    print(group)
    if(!is.null(allClusterRes[[pathway]][[group]])){
      soi <- as.character(unlist(allClusterRes[[pathway]][[group]]$specieslist))
      cluster <- as.numeric(unlist(allClusterRes[[pathway]][[group]]$cluster))
      enrichTable <- allEnrich[[group]][[cluster]]
      KoLongDf <- allGroups[allGroups$groupName==group,]
      res <- findKOInfo(pathway, soi, enrichTable, KoLongDf, group)
      table5 <- rbind(table5,cbind(group,pathway,res))
    }
  }
 
}
print("Saving the Table 5 output...") 
write.csv(table5, file="outputTable5.csv", row.names=F)

################################################################################
# For figure 4                                                                 #                                   
################################################################################
print("Creating the Figure 4...")

# Prepare the list of pathways of interest
pathwayOfInterest<- unique(table5$pathway)
print(pathwayOfInterest)

cdTable<- table5[grep("CD1|CD2", table5$group), c("group","pathway","KOlist")]
ucTable<- table5[grep("UC1|UC2", table5$group), c("group","pathway","KOlist")]

# For pathview function, 'pathway.id' parameter requires to have 'pathway IDs'  
# from KEGG Pathway database. Once you obtained 'pathwayOfInterest' varaible, 
# go visit (https://www.genome.jp/kegg/pathway.html). Then, search the name of 
# pathway in 'Entery keywords' bar. 

## 'Sulfur relay system' | KEGG pathway.ids = map04122 | CD and UC group 
cdSulfurRelaySys<- unlist(cdTable[cdTable$pathway==pathwayOfInterest[3],"KOlist"])
pathview(cdSulfurRelaySys, pathway.id="04122", species="ko", 
         out.suffix="cd_SulfurRelaySys")

ucSulfurRelaySys<- unlist(ucTable[ucTable$pathway==pathwayOfInterest[3],"KOlist"])
pathview(ucSulfurRelaySys, pathway.id="04122", species="ko", 
         out.suffix="uc_SulfurRelaySys")

################################################################################
# For Supplement figures                                                       #                                   
################################################################################
## 'ABC transporter' | KEGG pathway.ids = map02010 | CD group 
cdABCTransport<- unlist(cdTable[cdTable$pathway==pathwayOfInterest[1],"KOlist"])
pathview(cdABCTransport, pathway.id="02010", species="ko", 
         out.suffix="cd_ABCTransport")

## 'Two-component system' | KEGG pathway.ids = map02020 | CD group 
cdTwoCompSys<- unlist(cdTable[cdTable$pathway==pathwayOfInterest[2],"KOlist"])
pathview(cdTwoCompSys, pathway.id="02020", species="ko", 
         out.suffix="cd_TwoCompSys")