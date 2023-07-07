################################################################################
# cooccurrenceNetwork.R                                                        #
# Name: Suyeon Kim                                                             #
# Date: July-07-2023                                                           #
# Last updated: July-07-2023                                                   #
################################################################################
library(cooccur)

args = commandArgs(trailingOnly=TRUE)
print(args)

sourceDir <- args[1]
InFile1<- args[2]
InFile2<- args[3]

################################################################################
# Read the Input file into R.                                                  #
################################################################################
ibdMgx<- read.csv(InFile1, sep=",",header=T, stringsAsFactors = FALSE, 
                    check.names=FALSE, row.names=1)

ibdMeta<- read.csv(InFile2, sep=",",header=T, stringsAsFactors = FALSE, 
                     check.names=FALSE)

################################################################################
# Main program                                                                 #
################################################################################
## Subset metagenome abundance matrix based on the taxon of user interest. 
## In this study, we chose the 'species' level (tx_level = 6). 

source(paste0(sourceDir, "/", "mainPreProcessing.R"))
speciesAbundanceMat<- extractTaxonAndFilter(ibdMgx, tx_level=6, ibdMeta,
                                               frequentN = 1)

## Group samples based on its phenotypic information.                                         
controlSpecies<- getSamplebasedOnDiagnosis(speciesAbundanceMat, 
                                            ibdMeta, diagnosisType="Control")
cdSpecies<- getSamplebasedOnDiagnosis(speciesAbundanceMat, 
                                            ibdMeta, diagnosisType = "CD")
ucSpecies<- getSamplebasedOnDiagnosis(speciesAbundanceMat,
                                            ibdMeta, diagnosisType = "UC")

## Convert relative abundance matrix to presence/basence matrix                       
## (If proportion >0 => 1, else proportion = 0 => 0)                                            
controlSpecies[controlSpecies > 0]<- 1
cdSpecies[cdSpecies > 0]<- 1
ucSpecies[ucSpecies > 0]<- 1

## Create the species co-occurrence network                                               
## Note: 1) We suggest to save the result so it can be easily loaded and used in 
##          further our network analysis.
##       2) Please remove the sample if sample has no abundance value before 
##          conducting 'cooccur' function. 

## Co-occurrence network in control group 
controlCooccur<- cooccur(controlSpecies, thresh=T, spp_names=T)
save(controlCooccur, file="control.cooccur.bin")

# Co-occurrence network in CD group 
cdCooccur<- cooccur(cdSpecies, thresh = T, spp_names = T)
save(cdCooccur, file="cd.cooccur.bin")

# Co-occurrence network in UC group 
ucCooccur<- cooccur(ucSpecies, thresh=T, spp_names = T)
save(ucCooccur, file="uc.cooccur.bin")