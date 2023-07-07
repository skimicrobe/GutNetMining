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
InFile3<- args[4]
InFile4<- args[5]

################################################################################
# Read the Input file into R.                                                  #
################################################################################
ibd1.mgx<- read.csv(InFile1, sep=",",header=T, stringsAsFactors = FALSE, 
                    check.names=FALSE, row.names=1)

ibd1.meta<- read.csv(InFile2, sep=",",header=T, stringsAsFactors = FALSE, 
                     check.names=FALSE)

ibd2.mgx<- read.csv(InFile3, sep=",",header=T, stringsAsFactors = FALSE, 
                    check.names=FALSE, row.names=1 )

ibd2.meta<- read.csv(InFile4, sep=",",header=T, stringsAsFactors = FALSE, 
                     check.names=FALSE)

################################################################################
# Main program                                                                 #
################################################################################
## Subset metagenome abundance matrix based on the taxon of user interest. 
## In this study, we chose the 'species' level (tx_level = 6). 

source(paste0(sourceDir, "/", "mainPreProcessing.R"))
species.DefaultAb.ibd1 <- extractTaxonAndFilter(ibd1.mgx, tx_level=6, ibd1.meta,
                                                frequentN = 1)
species.DefaultAb.ibd2<- extractTaxonAndFilter(ibd2.mgx, tx_level=6, ibd2.meta,
                                               frequentN = 1)

## Group samples based on its phenotypic information.                                         
control.Species.Def.ibd1<- getSamplebasedOnDiagnosis(species.DefaultAb.ibd1, 
                                            ibd1.meta, diagnosisType="Control")
cd.Species.Def.ibd1<- getSamplebasedOnDiagnosis(species.DefaultAb.ibd1, 
                                            ibd1.meta, diagnosisType = "CD")
uc.Species.Def.ibd1<- getSamplebasedOnDiagnosis(species.DefaultAb.ibd1,
                                            ibd1.meta, diagnosisType = "UC")

control.Species.Def.ibd2<- getSamplebasedOnDiagnosis(species.DefaultAb.ibd2,
                                          ibd2.meta, diagnosisType = "Control")
cd.Species.Def.ibd2<- getSamplebasedOnDiagnosis(species.DefaultAb.ibd2,
                                          ibd2.meta, diagnosisType = "CD")
uc.Species.Def.ibd2<- getSamplebasedOnDiagnosis(species.DefaultAb.ibd2,
                                          ibd2.meta, diagnosisType = "UC")

## Convert relative abundance matrix to presence/basence matrix                       
## (If proportion >0 => 1, else proportion = 0 => 0)                                            
control.Species.Def.ibd1[control.Species.Def.ibd1 > 0]<- 1
cd.Species.Def.ibd1[cd.Species.Def.ibd1 > 0]<- 1
uc.Species.Def.ibd1[uc.Species.Def.ibd1 > 0]<- 1

control.Species.Def.ibd2[control.Species.Def.ibd2 > 0]<- 1
cd.Species.Def.ibd2[cd.Species.Def.ibd2 > 0]<- 1
uc.Species.Def.ibd2[uc.Species.Def.ibd2 > 0]<- 1

## Create the species co-occurrence network                                               
## Note: 1) We suggest to save the result so it can be easily loaded and used in 
##          further our network analysis.
##       2) Please remove the sample if sample has no abundance value before 
##          conducting 'cooccur' function. 

## IBD dataset 1 control group 
control1<- cooccur(control.Species.Def.ibd1, thresh=T, spp_names=T)
save(control1, file="ibdDataset1.control.cooccur.bin")

## IBD dataset 2 control group 
control2<- cooccur(control.Species.Def.ibd2, thresh = T, spp_names=T)
save(control2, file="ibdDataset2.control.cooccur.bin")

## IBD dataset 1 CD group 
filtCD1<- cd.Species.Def.ibd1[,-10] # This is because the sample has no abundance value
cd1<- cooccur(filtCD1, thresh = T, spp_names = T)
save(cd1, file="ibdDataset1.cd.cooccur.bin")

# IBD datset 2 CD group 
cd2<- cooccur(cd.Species.Def.ibd2, thresh = T, spp_names = T)
save(cd2, file="ibdDataset2.cd.cooccur.bin")

# IBD dataset 1 UC group 
filtUC1<- uc.Species.Def.ibd1[,-16] # This is because the sample has no abundance value
uc1<- cooccur(filtUC1, thresh = T, spp_names = T)
save(uc1, file="ibdDataset1.uc.cooccur.bin")

# IBD dataset 2 UC group 
uc2<- cooccur(uc.Species.Def.ibd2, thresh=T, spp_names = T)
save(uc2, file="ibdDataset2.uc.cooccur.bin")