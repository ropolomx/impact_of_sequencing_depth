# Rarefaction analysis of 2-4-8 study data
# Processing of AMR and Kraken data

# Load necessary packages

library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(metagenomeSeq)
library(vegan)

# Source script of utility functions

source('rarefaction_utility_functions.R')

# Read AMR and Kraken data

# Get the paths of the files that we want to process

# Update path accordingly

amrResultsFiles <- Sys.glob(file.path("~",
                                      "amr",
                                      "2-4-8_results",
                                      "bwa_aln",
                                      "*",
                                      "*",
                                      "*rarefied*.tab*"))


# Split the path names and extract the sample name only

amrResultsNames <- str_split(amrResultsFiles, pattern = "\\/")

# Extraction of sample name is being done by extracting the 8th element out of each list element.

# TODO: Explore if this can be done with a print statement

amrResultsNames <- amrResultsNames %>% 
  map(function(x){
    sample <- x[8]
    sample
  })

amrResultsNames <- unlist(amrResultsNames)

# Let's now read all the Coverage Sampler tabular files
# We are using readr (read_tsv)
# We are also using the list of sample names extracted in the previous function
# to set the names of the list elements
# This will make life so much easier!

amrResults <- amrResultsFiles %>%
  map(read_tsv) %>%
  set_names(nm=amrResultsNames)

# Need to add a column with sample name to the Coverage Sampler output
# Also need to add a column with sample depth

#TODO: Attempt to use map2 or pmap

amrArgs <- list(amrResults, amrResultsNames)

amrResults <- amrArgs %>%
  map2(function(x) {
   x$Sample <- rep(amrResultsNames, nrow(x))
   x
  })


# Join all the datasets into one dataframe that will be analyzed

amrResults <- do.call(amrResults, "rbind")

# TODO: If output is generated from Rarefaction Analyzer,
# see how it can be applied to the analysis of Kraken data, 
# although it seems that program can only be applied to process SAM files

# Convert AMR results to "tidy" format
# Kraken table is already in tidy format.

amrResultsTidy <- amrResults %>% 
  gather(Level, LevelName, c(1,6:8))

# Change the names of the AMR sequencing depths columns

amrResultsTidy$Depth <- str_replace(amrResultsTidy$Sample, "\\d+_","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "\\d$","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "full", "D1")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "half", "D0.5")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "quar*", "D0.25")

# Change the names of the Kraken sequencing depths columns

krakenResults$Depth <- str_replace(krakenResults$Sample, "_filtered2_report$","")
krakenResults$Depth <- str_replace(krakenResults$Depth, "^\\d+_","")
krakenResults$Depth <- str_replace(krakenResults$Depth, "\\d$","")
krakenResults$Depth <- str_replace(krakenResults$Depth, "full", "D1")
krakenResults$Depth <- str_replace(krakenResults$Depth, "half", "D0.5")
krakenResults$Depth <- str_replace(krakenResults$Depth, "quar*", "D0.25")

# Convert the LevelName to Factor

amrResultsTidy$LevelName <- as.factor(amrResultsTidy$LevelName)

# Keep results with over 80 % coverage ratio

amrResultsTidy <- amrResultsTidy %>% 
  filter(Coverage_Ratio >= 0.80)

# Vectors containing amr Levels and taxon levels to analyze
#amrLevels <- c("Class", "Mechanism", "Group", "Gene Id")

amrResultsTidy$Level <- factor(amrResultsTidy$Level, 
                               levels = c('Class', 'Mechanism', 'Group', 'Name'))

amrLevels <- levels(amrResultsTidy$Level)

taxonLevels <- c("D", "K", "C", "O", "F", "G", "S", "-")


# Split tidy data frame of AMR results according to the Level vector

amrResultsList <- amrResultsTidy %>% 
  split(.$Level) %>% #splitting dataframe using base R but with purrr's piping
  set_names(nm=amrLevels)

# Summarize results: add counts by Sample and AMR Level

amrResultsSummary <- lapply(amrResultsList, function(x){
  summarizeAMRbySample(x)
})

# Convert AMR results to wide format

amrResultsWide <- lapply(amrResultsSummary, function(x){
  widenAMR(x)
})

# Convert AMR results from wide format to matrix

amrResultsMat <- lapply(amrResultsWide, function(x){
  matrixAMR(x)
})


# Transpose matrices
# Awesome!!!

amrResultsMat <- lapply(amrResultsMat, function(x){
  
  t(x)
  
})

# Example: one specific rarefaction curve

#rarec(amrResultsMat[['Name']], 
#      step=100, 
#      sample=min(rowSums(amrResultsMat[['Name']])))


# Rarefaction curves for all members of the list
# This step requires optimization (parallel mapping or compiling, maybe?)

# purrr version

#amrRarefy <- map(amrResultsMat, function(x){
#  raremax <- min(rowSums(x))
#  rarecurve_ROP(x, step=5, sample=raremax)
#})

# mclapply version

amrRarefy <- mclapply(amrResultsMat, function(x){
  raremax <- min(rowSums(x))
  rarecurve(x, step=5, sample=raremax)
},mc.cores=3)

# Use microbenchmark to compare between the two approaches

#mbm <- microbenchmark(
#  mapping = map(amrResultsMat, function(x){
#    raremax <- min(rowSums(x))
#    rarecurve(x, step=5, sample=raremax)}),
#  multicore = mclapply(amrResultsMat, function(x){
#    raremax <- min(rowSums(x)) 
#    rarecurve(x, step=5, sample=raremax)}, mc.cores=12),
#  times=2
#)

# Results
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# mapping 432.4910 432.4910 439.0898 439.0898 445.6887 445.6887     2
# multicore 140.0049 140.0049 140.8782 140.8782 141.7515 141.7515     2

# mclapply wins!

# Rename list of rarefied data using the sample names
# Need to think of a better way to extract the sample names

samples <- rownames(amrResultsMat[["Class"]])

# Review syntax here
# Might need to add function to utility function list

amrRarefy <- mclapply(amrRarefy, function(x) set_names(x,samples), mc.cores=3)

# Generate list of dataframes
# Hacky, but it works. Need to make it cleaner and faster.

amrRarefy <- amrRarefy %>% 
  map(as_vector)

amrRarefyDF <- map(amrRarefy, function(x) as_tibble(x, attr(x, "names")))

# Split rownames in order to generate columns with useful information

amrRarefyDF <- mclapply(amrRarefyDF, function(x) {
  x$Sample <- row.names(x)
  x$Subsample <- str_extract(x$Sample, "N\\d+")
  x$Subsample <- as.numeric(str_replace(x$Subsample, "N",""))
  x$Depth <- str_replace(x$Sample, "\\.N.*$", "")
  x$Depth <- str_replace(x$Depth, "^\\d+_", "")
  x$SampleNumber <- str_replace(x$Sample, "_.*", "")
  x$Sample <- str_replace(x$Sample, "\\.N.*$", "")
  x
}, 
mc.cores=3
)

# Generate one single dataframe and create column for AMR Level

amrRarefyDF <- do.call("rbind", amrRarefyDF)
amrRarefyDF$AMRLevel <- row.names(amrRarefyDF)
amrRarefyDF$AMRLevel <- str_extract(amrRarefyDF$AMRLevel, "^\\w+")

