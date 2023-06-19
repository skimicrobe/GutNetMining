#!/usr/bin/r

################################################################################
# Script: keyElementLevelAnalysis.R                                            #
# Name: Suyeon Kim                                                             #
# Created: May-15-2023                                                         #
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
# In this function, measure centrality in the graph using various centrality 
# measurements
measureCentrality<- function(inputGraph){

  # Assortivity centrality 
  assort<- assortativity.degree(inputGraph, directed=F)
  # Closeness centrality
  clos<- closeness(inputGraph, mode = "all", weights=NA)
  # Betweenness centrality
  betw<- betweenness(inputGraph, directed = F, weights=NA)
  # Degree centrality 
  alldeg<- degree(inputGraph, v=V(inputGraph), mode="all")
  # Eigenvector centrality 
  eigenV<- eigen_centrality(inputGraph, directed = F, weights = NA)
  # Hub centrality 
  hub<- hub.score(inputGraph)$vector
  
  centralityMeasure<- data.frame(cbind(assort,clos, betw, alldeg, 
                                       "eig"= eigenV$vector, hub))
  
  return(centralityMeasure)
}

# Get the top 'N' number of scores by each measurement 
getTopNscores<- function(centralityTbl, mcol, topN){
  
  sortedtbl<- head(centralityTbl[order(-centralityTbl[,mcol]),], n=topN)
  
  return(sortedtbl)
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

## Perform centrality analysis 
ctrl1NodeCentral<- measureCentrality(ctrl1graph)
ctrl2NodeCentral<- measureCentrality(ctrl2graph)

cd1NodeCentral<- measureCentrality(cd1graph)
cd2NodeCentral<- measureCentrality(cd2graph)

uc1NodeCentral<- measureCentrality(uc1graph)
uc2NodeCentral<- measureCentrality(uc2graph)

## Key element-level analysis: Top 'N' hub. Here we choose N for 5. 
ctrl1Hub<- getTopNscores(ctrl1NodeCentral, mcol="hub", topN=5)
ctrl2Hub<- getTopNscores(ctrl2NodeCentral, mcol="hub", topN=5)

cd1Hub<- getTopNscores(cd1NodeCentral, mcol="hub", topN=5)
cd2Hub<- getTopNscores(cd2NodeCentral, mcol="hub", topN=5)

uc1Hub<- getTopNscores(uc1NodeCentral, mcol="hub", topN=5)
uc2Hub<- getTopNscores(uc2NodeCentral, mcol="hub", topN=5)

## Key element-level analysis: Top 'N' betweenness. Here we choose N for 5. 
ctrl1Btw<- getTopNscores(ctrl1NodeCentral, mcol="betw", topN=5)
ctrl2Btw<- getTopNscores(ctrl2NodeCentral, mcol="betw", topN=5)

cd1Btw<- getTopNscores(cd1NodeCentral, mcol="betw", topN=5)
cd2Btw<- getTopNscores(cd2NodeCentral, mcol="betw", topN=5)

uc1Btw<- getTopNscores(uc1NodeCentral, mcol="betw", topN=5)
uc2Btw<- getTopNscores(uc2NodeCentral, mcol="betw", topN=5)

## Get list of KO-IDs based on our species of interest from top 5 hub. 
ctrl1KoHub<- getKoListForSpeciesOfInterest(ko1, groupname = "Control", 
                                           rownames(ctrl1Hub))
ctrl2KoHub<- getKoListForSpeciesOfInterest(ko2, groupname = "Control", 
                                           rownames(ctrl2Hub))

cd1KoHub<- getKoListForSpeciesOfInterest(ko1, groupname = "CD", 
                                         rownames(cd1Hub))
cd2KoHub<- getKoListForSpeciesOfInterest(ko2, groupname = "CD", 
                                         rownames(cd2Hub))

uc1KoHub<- getKoListForSpeciesOfInterest(ko1, groupname = "UC", 
                                         rownames(uc1Hub))
uc2KoHub<- getKoListForSpeciesOfInterest(ko2, groupname = "UC", 
                                         rownames(uc2Hub))

## Get list of KO-IDs based on our species of interest from top 5 betweenness. 
ctrl1KoBtw<- getKoListForSpeciesOfInterest(ko1, groupname = "Control", 
                                           rownames(ctrl1Btw))
ctrl2KoBtw<- getKoListForSpeciesOfInterest(ko2, groupname = "Control", 
                                           rownames(ctrl2Btw))

cd1KoBtw<- getKoListForSpeciesOfInterest(ko1, groupname = "CD", 
                                         rownames(cd1Btw))
cd2KoBtw<- getKoListForSpeciesOfInterest(ko2, groupname = "CD", 
                                         rownames(cd2Btw))

uc1KoBtw<- getKoListForSpeciesOfInterest(ko1, groupname = "UC", 
                                         rownames(uc1Btw))
uc2KoBtw<- getKoListForSpeciesOfInterest(ko2, groupname = "UC", 
                                         rownames(uc2Btw))


## Perform the enrichment analysis the list of KO-IDS from top 5 hub 
ctrl1EnrichHub<- data.frame(enrichKO(unique(unlist(ctrl1KoHub$ko)), 
                                     universe = ko1$Koid))
ctrl2EnrichHub<- data.frame(enrichKO(unique(unlist(ctrl2KoHub$ko)),
                                     universe = ko2$Koid))

cd1EnrichHub<- data.frame(enrichKO(unique(unlist(cd1KoHub$ko)),
                                   universe = ko1$Koid))
cd2EnrichHub<- data.frame(enrichKO(unique(unlist(cd2KoHub$ko)),
                                   universe = ko2$Koid))

uc1EnrichHub<- data.frame(enrichKO(unique(unlist(uc1KoHub$ko)),
                                   universe = ko1$Koid))
uc2EnrichHub<- data.frame(enrichKO(unique(unlist(uc2KoHub$ko)),
                                   universe = ko2$Koid))

## Perform the enrichment analysis the list of KO-IDS from top 5 betweenness  
ctrl1EnrichBtw<- data.frame(enrichKO(unique(unlist(ctrl1KoBtw$ko)), 
                                     universe = ko1$Koid))
ctrl2EnrichBtw<- data.frame(enrichKO(unique(unlist(ctrl2KoBtw$ko)),
                                     universe = ko2$Koid))

cd1EnrichBtw<- data.frame(enrichKO(unique(unlist(cd1KoBtw$ko)),
                                   universe = ko1$Koid))
cd2EnrichBtw<- data.frame(enrichKO(unique(unlist(cd2KoBtw$ko)),
                                   universe = ko2$Koid))

uc1EnrichBtw<- data.frame(enrichKO(unique(unlist(uc1KoBtw$ko)),
                                   universe = ko1$Koid))
uc2EnrichBtw<- data.frame(enrichKO(unique(unlist(uc2KoBtw$ko)),
                                   universe = ko2$Koid))


## Comparative Analysis: 1) Compare the corresponding lists of enriched pathways 
##                          between disease phenotype and healthy phenotype 
##                       2) Identify the pathways unique to 'CD' group that are 
##                          compared across two datasets. 

resCDhub<- sort(intersect(setdiff(cd1EnrichHub$Description, 
                                  ctrl1EnrichHub$Description),
                          setdiff(cd2EnrichHub$Description, 
                                  ctrl2EnrichHub$Description)))

resUChub<- sort(intersect(setdiff(uc1EnrichHub$Description, 
                                  ctrl1EnrichHub$Description),
                          setdiff(uc2EnrichHub$Description, 
                                  ctrl2EnrichHub$Description)))

## print the results of hub scores for CD and UC group 
if(length(resCDhub) == 0){
  cat(paste0("Enriched Pathways in CD group (Top 5 Hub nodes): ", "None", "\n"))
} else{
  cat(paste0("Enriched Pathways in CD group (Top 5 Hub nodes) : ", resCDhub, 
             "\n"))
}

if(length(resUChub) == 0){
  cat(paste0("Enriched Pathways in UC group (Top 5 Hub nodes): ", "None", "\n"))
} else{
  cat(paste0("Enriched Pathways in UC group (Top 5 Hub nodes) : ", resUChub, 
             "\n"))
}

## Comparative Analysis: 1) Compare the corresponding lists of enriched pathways 
##                          between disease phenotype and healthy phenotype 
##                       2) Identify the pathways unique to 'UC' group that are 
##                          compared across two datasets. 

resCDBetweenness<- sort(intersect(setdiff(cd1EnrichBtw$Description, 
                                          ctrl1EnrichBtw$Description),
                                  setdiff(cd2EnrichBtw$Description,
                                          ctrl2EnrichBtw$Description)))

resUCBetweenness<- sort(intersect(setdiff(uc1EnrichBtw$Description, 
                                          ctrl1EnrichBtw$Description),
                                  setdiff(uc2EnrichBtw$Description,
                                          ctrl2EnrichBtw$Description)))

## print the results of betweenness scores for CD and UC group 
if(length(resCDBetweenness) == 0){
  cat(paste0("Enriched Pathways in CD group (Top 5 Betweenness) : ", 
             "None", "\n"))
} else{
  cat(paste0("Enriched Pathways in CD group (Top 5 Betweenness) : ", 
             resCDBetweenness, "\n"))
}

if(length(resUCBetweenness) == 0){
  cat(paste0("Enriched Pathways in UC group (Top 5 Betweenness) : ", 
             "None", "\n"))
} else{
  cat(paste0("Enriched Pathways in UC group (Top 5 Betweenness) : ", 
             resUCBetweenness, "\n"))
}