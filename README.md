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
![overviewPipeline.png](overviewPipeline_IP_OP.png)

### Network analysis 

#### _Global Network analysis_
Run the 'globalNetworkAnaysis.R' script in github_src directory with the following parameters: 

#### _Community-level analysis_ 
Run the 'communityLevelAnalysis.R'
#### _Key element-level analysis_ 
Run the 'keyElementLevelAnalysis.R'
```
Rscript ./github_src/keyElementLevelAnalysis.R output/ctrl1.cooccur.network.csv output/ctrl2.cooccur.network.csv output/cd1.cooccur.network.csv output/cd2.cooccur.network.csv output/uc1.cooccur.network.csv output/uc2.cooccur.network.csv output/s1_preprocess/IBD1_KOdata.csv output/s1_preprocess/IBD2_KOdata.csv
```
### Comparative analysis
