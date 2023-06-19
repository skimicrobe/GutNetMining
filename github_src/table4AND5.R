#!/usr/bin/r

################################################################################
# tabe4AND5.R                                                                  #
# Name: Suyeon Kim                                                             #
# Date: May-31-2023                                                            #
# Last updated: June-19-2023                                                   #
# Note: This script is to produce the results for Table 4 and 5.               #
################################################################################

args = commandArgs(trailingOnly=TRUE)
print(args)

sourceDir <- args[1]
InFile1<- args[2]
InFile2<- args[3]

################################################################################
# List of functions                                                            #
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
        matchedGene<- ko.long.ibd[which(KoLongDf$Koid %in% geneList & 
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
                                    KoLongDf$GROUP == group & 
                                    KoLongDf$cooccur == TRUE),]
  speciesGenes <- matchedGene[matchedGene$Species %in% soi,]
  
  speciesKO <- aggregate(speciesGenes$Koid,by=list(speciesGenes$Species),function(x) list(unique(x))) 
  colnames(speciesKO) <- c("species","KOlist")
  speciesKO$count <- unlist(lapply(speciesKO$KOlist,length))
  return(speciesKO)
}

################################################################################
# Read the Input file into R.                                                  #
################################################################################
# This should result in a dataframe called 'ko1' that long format of KOs table 
# with species and sampleIDs. 
ko1<- read.csv(InFile1, sep=",", header=T)
# This should result in a dataframe called 'ko2' that long format of KOs table 
# with species and sampleIDs.
ko2<- read.csv(InFile2, sep=",", header=T)

################################################################################
# Main Program                                                                 #
################################################################################
source(paste0(sourceDir, "/", "communityLevelAnalysis.R"))

## Make a list of co-occurring species from the co-occurrence networks
listSpeciesCD1<- listOfcooccurSpecies(cd1)
listSpeciesCD2<- listOfcooccurSpecies(cd2)

listSpeciesUC1<- listOfcooccurSpecies(uc1)
listSpeciesUC2<- listOfcooccurSpecies(uc2)

## Create a new column onto the KEGG Orthologue table (long format) to add      
## boolean if species which are associated with KO are identified in the        
## co-occurrence network                                                        
cd1Group<- ko1[ko1$GROUP == "CD",]
cd1Group$cooccur<- cd1Group$Species %in% listSpeciesCD1

cd2Group<- ko2[ko2$GROUP == "CD",]
cd2Group$cooccur<- cd2Group$Species %in% listSpeciesCD2

uc1Group<- ko1[ko1$GROUP == "UC",]
uc1Group$cooccur<- ko1$Species %in% listSpeciesUC1

uc2Group<- ko2[ko2$GROUP == "UC",]
uc2Group$cooccur<- ko2$Species %in% listSpeciesUC2

## Find common enriched pathways across two IBD datasets                        
listPathwaysCD<- intersect(cd1minusctrl1, cd2minusctrl2)
listPathwaysUC<- intersect(uc1minusctrl1, uc2minusctrl2)

## Find the number of co-occurring species on these common enriched pathways                                                                                   
cdClusterRes <- list()
ucClusterRes <- list()

for(pathway in listPathwaysCD){
  cdClusterRes[[pathway]][["CD1"]] <- findEnrichmentTable(pathway, cd1Enrich, 
                                                          listSpeciesCD1, 
                                                          cd1Comm, group="CD")
  cdClusterRes[[pathway]][["CD2"]] <- findEnrichmentTable(pathway, cd2Enrich, 
                                                          listSpeciesCD2, 
                                                          cd2Comm, group="CD")
}

for(pathway in listPathwaysUC){
  ucClusterRes[[pathway]][["UC1"]] <- findEnrichmentTable(pathway, uc1Enrich, 
                                                            listSpeciesUC1, 
                                                            uc1Comm, group="UC")
  ucClusterRes[[pathway]][["UC2"]] <- findEnrichmentTable(pathway, uc2Enrich, 
                                                            listSpeciesUC2, 
                                                            uc2Comm, group="UC")
}

print(cdClusterRes)
print(ucClusterRes)


## This is result for Table 5. List of co-occurring species and count of KOs 
## contributing to the Sulfur Relay System pathway.
cd1KoRes <- findKOInfo("Sulfur relay system", 
                         unlist(cdClusterRes$`Sulfur relay system`$CD1$specieslist), 
                         cd1Enrich[[unlist(cdClusterRes$`Sulfur relay system`$CD1$cluster)]],
                         listSpeciesCD1, 
                         group = "CD")

cd2KoRes <- findKOInfo("Sulfur relay system", 
                         unlist(cdClusterRes$`Sulfur relay system`$CD2$specieslist), 
                         cd2Enrich[[unlist(cdClusterRes$`Sulfur relay system`$CD2$cluster)]],
                         listSpeciesCD2, 
                         group = "CD")


uc1KoRes <- findKOInfo("Sulfur relay system", 
                         unlist(ucClusterRes$`Sulfur relay system`$UC1$specieslist), 
                         uc1Enrich[[unlist(ucClusterRes$`Sulfur relay system`$UC1$cluster)]],
                         listSpeciesUC1, 
                         group = "UC")

uc2KoRes <- findKOInfo("Sulfur relay system", 
                         unlist(ucClusterRes$`Sulfur relay system`$UC2$specieslist), 
                         uc2Enrich[[unlist(ucClusterRes$`Sulfur relay system`$UC2$cluster)]],
                         listSpeciesUC2, 
                         group = "UC")
