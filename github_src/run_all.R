#!/usr/bin/r

################################################################################
# Script: communityLevelAnalysis.R                                             #
# Name: Suyeon Kim                                                             #
# Created: Feb-01-2024                                                         #
# Last edited: Feb-01-2024                                                     #
################################################################################

# Check the number of command-line arguments
if (length(commandArgs(trailingOnly = TRUE)) < 8) {
  stop("Usage: Rscript run_all.R source directory path input_file1 input_file2 
       input_file3 input_file4 input_file5 input_file6 input_file7 input_file8")
}

# Extract command-line arguments
sourceDir<- commandArgs(trailingOnly = TRUE)[1]
InFile1 <- commandArgs(trailingOnly = TRUE)[2]
InFile2 <- commandArgs(trailingOnly = TRUE)[3]
InFile3 <- commandArgs(trailingOnly = TRUE)[4]
InFile4 <- commandArgs(trailingOnly = TRUE)[5]
InFile5 <- commandArgs(trailingOnly = TRUE)[6]
InFile6 <- commandArgs(trailingOnly = TRUE)[7]
InFile7 <- commandArgs(trailingOnly = TRUE)[8]
InFile8 <- commandArgs(trailingOnly = TRUE)[9]

# Source or run the individual scripts with the specified input files
source(paste0(sourceDir, "globalNetworkAnalysis.R"), local = TRUE)
source(paste0(sourceDir, "communityLevelAnalysis.R"), local = TRUE)
source(paste0(sourceDir, "keyElementLevelAnalysis.R"), local = TRUE)
