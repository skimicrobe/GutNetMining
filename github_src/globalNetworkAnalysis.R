#!/usr/bin/r

################################################################################
# Script: globalNetworkAnalysis.R                                              #
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

## Get list of KO-IDs based on our species of interest from global network  
ctrl1KO<- getKoListForSpeciesOfInterest(ko1, groupname = "Control", 
                                           unique(c(ctrl1$sp1, ctrl1$sp2)))
ctrl2KO<- getKoListForSpeciesOfInterest(ko2, groupname = "Control", 
                                        unique(c(ctrl2$sp1, ctrl2$sp2)))

cd1KO<- getKoListForSpeciesOfInterest(ko1, groupname = "CD", 
                                         unique(c(cd1$sp1, cd1$sp2)))
cd2KO<- getKoListForSpeciesOfInterest(ko2, groupname = "CD", 
                                         unique(c(cd2$sp1, cd2$sp2)))

uc1KO<- getKoListForSpeciesOfInterest(ko1, groupname = "UC", 
                                         unique(c(uc1$sp1, uc1$sp2)))
uc2KO<- getKoListForSpeciesOfInterest(ko2, groupname = "UC", 
                                         unique(c(uc2$sp1, uc2$sp2)))

## Perform the enrichment analysis the list of KO-IDS from global network 
ctrl1EnrichNW<- data.frame(enrichKO(unique(unlist(ctrl1KO$ko)), 
                                     universe = ko1$Koid))
ctrl2EnrichNW<- data.frame(enrichKO(unique(unlist(ctrl2KO$ko)),
                                     universe = ko2$Koid))

cd1EnrichNW<- data.frame(enrichKO(unique(unlist(cd1KO$ko)), 
                                    universe = ko1$Koid))
cd2EnrichNW<- data.frame(enrichKO(unique(unlist(cd2KO$ko)),
                                    universe = ko2$Koid))

uc1EnrichNW<- data.frame(enrichKO(unique(unlist(uc1KO$ko)), 
                                  universe = ko1$Koid))
uc2EnrichNW<- data.frame(enrichKO(unique(unlist(uc2KO$ko)),
                                  universe = ko2$Koid))

## Find common enriched pathways across two IBD datasets 
resCDnetwork<- sort(intersect(setdiff(cd1EnrichNW$Description, 
                                  ctrl1EnrichNW$Description),
                          setdiff(cd2EnrichNW$Description, 
                                  ctrl2EnrichNW$Description)))

resUCnetwork<- sort(intersect(setdiff(uc1EnrichNW$Description, 
                                  ctrl1EnrichNW$Description),
                          setdiff(uc2EnrichNW$Description, 
                                  ctrl2EnrichNW$Description)))

## print the results of hub scores for CD and UC group 
if(length(resCDnetwork) == 0){
  cat(paste0("Enriched Pathways in CD group by global network: ", "None", "\n"))
} else{
  cat(paste0("Enriched Pathways in CD group by global network : ", resCDnetwork, 
             "\n"))
}

if(length(resUCnetwork) == 0){
  cat(paste0("Enriched Pathways in UC group by global network: ", "None", "\n"))
} else{
  cat(paste0("Enriched Pathways in UC group by global network : ", resUCnetwork, 
             "\n"))
}
