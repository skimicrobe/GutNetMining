# GutNetMining User Manaul 

This repository contains all the code and processed datsets necessary to reproduce the results
Contact the author at suyeonkim [at] unomaha.edu. This version has been tested for OSX. 

----------------------------------------------------------------------
## Features 
1. Pathway-based profile based on species co-occurring community 
2. Finding a robust pathway which are enriched for captured co-occurring community 
3. Easy to implement
----------------------------------------------------------------------
## Requirements
Software R 

----------------------------------------------------------------------
## Workflows
### Main workflow
![Overview](https://github.com/skimicrobe/GutNetMining/blob/main/Overviewpipeline.png)

### How to construct the species co-occurrence network 
We mainly utilize the R package called 'cooccur' to create the co-occurrence network. More details to this package: _['cooccur'](https://www.jstatsoft.org/article/view/v069c02)_. 
_Note: We will construct the co-occurrence network using 'cooccur' function. Input data must be presence/absence matrix. Abundance values greater than zero are considered "present"._

Run the 'cooccurrenceNetwork.R' script in github_src directory with the following parameters:
```
./cooccurrenceNetwork.R MGX_TAXONPROFILE, METADATA
```
A sample run using the example data. 
```
Rscript ./github_src/cooccurrenceNetwork.R ./github_src/ ./sampledata/rawdata/IBD_mgxdata.csv ./sampledata/rawdata/IBD_metadata.csv 
```
_Note: This sample dataset is publicly available ([PRJNA400072](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP129027&o=acc_s%3Aa)). Our provided data is pre-processed and filtered by our sample filtering process._

### Network analysis 
Run the 'NETWORK_ANALYSIS.R' script in github_src directory with the following parameters: 
```
./NETWORK_ANALYSIS.R SRC_DIR CONTROL_NETW_1 CONTROL_NETW_2 CD_NETW_1 CD_NETW_2 UC_NETW_1 UC_NETW_2 KO_1 KO_2
```
* NETWORK_ANALYSIS.R: Perform the following scripts( _glboalNetworkAnalysis.R_, _communityLevelAnalysis.R_, and _keyElementLevelAnalysis.R_ )

* SRC_DIR: directory containing other necessary scripts for performing these scripts.
 
* CONTROL_NETW_1: file with the species co-occurrence network at the level of significance less than 0.05 in Control (dataset 1)

* CONTROL_NETW_2: file with the species co-occurrence network at the level of significance less than 0.05 in Control (dataset 2)

* CD_NETW_1: file with the species co-occurrence network at the level of significance less than 0.05 in CD (dataset 1)

* CD_NETW_2: file with the species co-occurrence network at the level of significance less than 0.05 in CD (dataset 2)

* UC_NETW_1: file with the species co-occurrence network at the level of significance less than 0.05 in UC (dataset 1)

* UC_NETW_2: file with the species co-occurrence network at the level of significance less than 0.05 in UC (dataset 2)

* KO_1: file with the long format KEGG Orthologue containing KO-IDs, Genus, Species, Sample-IDs, GroupName (control, cd, and uc) from dataset 1

* KO_2: file with the long format KEGG Orthologue containing KO-IDs, Genus, Species, Sample-IDs, GroupName (control, cd, and uc) from dataset 2

A sample run using the example data. The level of significance level for the species co-occurrence network can be changed by the users.

#### _Global Network analysis_

```
Rscript ./github_src/globalNetworkAnaysis.R ./github_src/ ./sampledata/cooccurNetwork/ctrl1.cooccur.network.csv ./sampledata/cooccurNetwork/ctrl2.cooccur.network.csv ./sampledata/cooccurNetwork/cd1.cooccur.network.csv ./sampledata/cooccurNetwork/cd2.cooccur.network.csv ./sampledata/cooccurNetwork/uc1.cooccur.network.csv ./sampledata/cooccurNetwork/uc2.cooccur.network.csv ./sampledata/ko/IBD1_KOdata.csv ./sampledata/ko/IBD2_KOdata.csv
```
#### _Community-level analysis_ 
Run the 'communityLevelAnalysis.R'
```
Rscript ./github_src/communityLevelAnalysis.R ./github_src/ ./sampledata/cooccurNetwork/ctrl1.cooccur.network.csv ./sampledata/cooccurNetwork/ctrl2.cooccur.network.csv ./sampledata/cooccurNetwork/cd1.cooccur.network.csv ./sampledata/cooccurNetwork/cd2.cooccur.network.csv ./sampledata/cooccurNetwork/uc1.cooccur.network.csv ./sampledata/cooccurNetwork/uc2.cooccur.network.csv ./sampledata/ko/IBD1_KOdata.csv ../sampledata/ko/IBD2_KOdata.csv
```
#### _Key element-level analysis_ 
Run the 'keyElementLevelAnalysis.R'
```
Rscript ./github_src/keyElementLevelAnalysis.R ./github_src/ ./sampledata/cooccurNetwork/ctrl1.cooccur.network.csv ./sampledata/cooccurNetwork/ctrl2.cooccur.network.csv ./sampledata/cooccurNetwork/cd1.cooccur.network.csv ./sampledata/cooccurNetwork/cd2.cooccur.network.csv ./sampledata/cooccurNetwork/uc1.cooccur.network.csv ./sampledata/cooccurNetwork/uc2.cooccur.network.csv ./sampledata/ko/IBD1_KOdata.csv ../sampledata/ko/IBD2_KOdata.csv
```
### Comparative analysis
```
Rscript ./github_src/table4AND5.R ./github_src/ ../sampledata/ko/IBD1_KOdata.csv ../sampledata/ko/IBD2_KOdata.csv
```
