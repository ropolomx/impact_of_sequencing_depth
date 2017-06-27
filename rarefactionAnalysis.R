# Rarefaction analysis of 2-4-8 study data

# Load necessary packages

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(metagenomeSeq)
library(vegan)

# Source script with functions
source('rarefaction_functions.R')

# Clean AMR data and convert it to "tidy" format

amrResults <- read.csv('noelle/AMR/amr_new_dataframe_ROP.csv')

# Tidy results and change the names of the sequencing depths
amrResultsTidy <- amrResults %>% gather(Level, LevelName, c(1,6:8))
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Sample, "\\d+_","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "\\d$","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "full", "D1")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "half", "D0.5")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "quar*", "D0.25")

# Convert the LevelName to Factor
amrResultsTidy$LevelName <- as.factor(amrResultsTidy$LevelName)

# Keep results with over 80 % coverage ratio
amrResultsTidy <- amrResultsTidy %>% filter(Coverage_Ratio >= 0.80)

# Vector containing amrLevels to analyze
amrLevels <- c("Class", "Mechanism", "Group", "Gene Id")

# Split tidy data frame of AMR results according to the Level vector
amrResultsList <- split(amrResultsTidy,amrResultsTidy$Level)

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

# Lapply version

#amrRarefy <- lapply(amrResultsMat2, function(x){
#raremax <- min(rowSums(x))
#rarecurve(x, step=5, sample=raremax)
#})

# purrr version

amrRarefy <- map(amrResultsMat2, function(x){
  raremax <- min(rowSums(x))
  rarecurve(x, step=5, sample=raremax)
})

# mclapply version

amrRarefy <- mclapply(amrResultsMat, function(x){
  raremax <- min(rowSums(x))
  rarecurve(x, step=5, sample=raremax)
}
,mc.cores=10)




# Unlisting rarefied data and isolating DFs
amrRarefyDF2 <- lapply(amrRarefy, unlist)

amrRarefyDF2 <- lapply(amrRarefyDF2, function(x){
  data.frame(otus=x,subsample=attr(x, "names"))
})

# Isolate Class DF
amrRarefyClass <- amrRarefyDF2[['Class']]

# Avoid the work below by making sure we have a named list

amrRarefyClassDF$Test <- ifelse(amrRarefyClassDF$subsample == "N1", 
                                amrRarefyClassDF$Test == sapply(sampleNames, function(x){x}), NA)


amr_results_class_wide <- amr_results_class %>%
    spread(Class, Hits, fill = 0)

amr_results_class_wide <- amr_results_class %>%
  spread(Sample, Hits, fill = 0)

table(annotations$class)

amr_results_class_wide <- amr_results_class %>% 
  spread(Class, Hits, fill = 0)

amr_results_class_mat <- amr_results_class_wide[,2:ncol(amr_results_class_wide)]
amr_class_norm <- data.frame(MRcounts(amr_class_experiment, norm=T))
amr_class_norm <- t(amr_class_norm)
row.names(amr_class_norm) <- str_replace(row.names(amr_class_norm), "X", "")
