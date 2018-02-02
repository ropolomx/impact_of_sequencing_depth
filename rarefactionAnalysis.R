# Rarefaction analysis of 2-4-8 study data
# Analysis of alpha-diversity and species richness

# Load packages -----------------------------------------------------------

library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(foreach)
library(ggplot2)
library(metagenomeSeq)
library(vegan)
library(scales)
library(PMCMR) # For statistical tests

# Source utility functions ------------------------------------------------

source('rarefaction_utility_functions.R')


# Colour palettes ---------------------------------------------------------

# Same colour palette as the matplotlib_venn default color palette
# Consider adding them to the utility functions script

vennPalette <- c("#F8766D", "#00BA38", "#619CFF")
vennPalette <- rev(vennPalette)

# Colourblind friendly palette

cbPalette <- c("#FFCD48", # Mango from Crayola palette
               dichromat::colorschemes$Categorical.12[8], # Blue from dichromat package
                 dichromat::colorschemes$Categorical.12[12]) # Red from dichromat package

# Load and filter AMR and Kraken data -------------------------------------

# Results generated with Coverage Sampler and filtered with Python-Pandas
# Filtering involved keeping results with gene fraction >= 75% and 
# removing all those results with genes that require SNP confirmation.

amrResultsFiltered <- read_csv(file.path(
  '~',
  'amr',
  '2-4-8_results',
  '2_4_8_study_RZ',
  'amrResults_Aug2017_75_gene_frac/cov_sampler_parsed',
  'amrFiltered_75_genefrac.csv')
)

amrReadstoHitRatio <- read_tsv('~/aafc/amr/amrplusplus_rarefaction_analysis/2_4_8_study_RZ/Results_Aug2017/reads_and_hits.tsv')

# Read Kraken concatenated and filtered file (no Eukaryotes, and no PhiX)

krakenResultsFiltered <- read.delim(file.path(
  '~',
  'amr',
  '2-4-8_results',
  '2_4_8_study_RZ',
  'krakenResults_Aug2017',
  'allKraken_FHQ',
  'kraken_filtered',
  'krakenConcat.tsv'),
  stringsAsFactors = FALSE
)

# Remove D0.25_seqtk data from filtered data frame

amrResultsFiltered <- amrResultsFiltered %>% 
  filter(Sample_type != "D0.25_seqtk")

# Replace "Name" with "Gene" in headers of the filtered data frame

names(amrResultsFiltered) <- str_replace(names(amrResultsFiltered), "Name", "Gene")

# Convert AMR filtered results to "tidy" format in order to have a column with 
# all AMR levels

amrResultsTidy <- amrResultsFiltered %>% 
  gather(Level, LevelName, c(1,6:8))

# Change the names of the AMR sequencing depths columns
# Make this more functional

amrResultsTidy$Depth <- str_replace(amrResultsTidy$Sample, "\\d+_","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "\\d$","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "full", "D1")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "half", "D0.5")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "quar*", "D0.25")

# Load filtered AMR data --------------------------------------------------

# Start here if AMR results have been filtered already

# Tidy the dataset

amrResultsTidy <- amrResultsFiltered %>% 
  gather(Category, CategoryName, c(1:4))

amrResultsTidy$Category <- factor(amrResultsTidy$Category, 
                                  levels = c('Class', 
                                             'Mechanism', 
                                             'Group', 
                                             'Gene'))

amrCategories <- levels(amrResultsTidy$Category)

krakenTaxa <- levels(krakenResultsFiltered$TaxRank)

# Split results by categories ---------------------------------------------

# Will adopt the split, apply, combine strategy
# Using the names of the AMR categories and taxonomic ranks as names of the 
# elements of the lists

amrResultsList <- amrResultsTidy %>% 
  split(.$Category) %>% #splitting dataframe using base R but with purrr's magrittr piping
  set_names(nm=amrCategories)

krakenResultsList <- krakenResultsFiltered %>% 
  split(.$TaxRank) %>% 
  set_names(nm=krakenTaxa)

# Summarize results(sum) --------------------------------------------------

amrResultsSummary <- lapply(amrResultsList, function(x){
  summarizeAMRbyCategory(x)
})

krakenResultsSummary <- lapply(krakenResultsList, function(x){
  summarizeKrakenbyTaxID(x)
})

# Convert results to wide format ------------------------------------------

# Potentially can use the same function and a list/vector of AMR and Kraken 
# results to avoid duplication

amrResultsWide <- lapply(amrResultsSummary, function(x){
  widenAMR(x)
})

krakenResultsWide <- lapply(krakenResultsSummary, function(x){
  widenKraken(x)
})

# Convert wide to matrix --------------------------------------------------

amrResultsMat <- lapply(amrResultsWide, function(x){
  matrixAMR(x)
})

krakenResultsMat <- lapply(krakenResultsWide, function(x){
  matrixKraken(x)
})

# Transpose matrices
# Awesome!!!

amrResultsMat <- lapply(amrResultsMat, function(x){
  
  t(x)
  
})

krakenResultsMat <- lapply(krakenResultsMat, function(x){
  
  t(x)
  
})


# CSS Normalization -------------------------------------------------------

# Transform dataframes into wide format
# Fill NAs with a value of zero

# Splitting by Depth

amrResultsbyDepth <- amrResultsFiltered %>%
  split(.$Sample_type)

amrAnalytical <- lapply(amrResultsbyDepth, function(x){
  amrAnalyticWide <- x %>% 
    select(Gene, Hits, Sample) %>%
    spread(Sample, Hits, fill = 0, convert = TRUE)
  return(amrAnalyticWide)
})
  
  
#   lapply(amrResultsbyDepth, function(x){
#   amrResultsWide <- x %>% 
#     group_by(Sample) %>%
#     select(Sample, Hits, CategoryName) %>%
#     mutate(id=1:n()) %>%
#     spread(key=Sample, value=Hits, fill = 0) %>% 
#     select(-id)
#   return(amrResultsWide)
# })
#  

krakenResultsbyDepth <- krakenResultsFiltered %>%
  split(.$Sample_Type)

krakenAnalytical <- lapply(krakenResultsbyDepth, function(x){
  krakenAnalyticWide <- x %>% 
    select(TaxID, CladeReads, Sample) %>%
    spread(Sample, CladeReads, fill = 0, convert = TRUE)
  return(krakenAnalyticWide)
})
  
amrAnalytical <- lapply(amrAnalytical, function(x){
  amrAnaMat <- matrixAMRanalytical(x)
  return(amrAnaMat)
})

krakenAnalytical <- lapply(krakenAnalytical, function(x){
  krakenAnaMat <- matrixKraken(x)
  return(krakenAnaMat)
})

krakenTaxInfo <- krakenResultsFiltered  %>% select(2:ncol(.))

krakenTaxInfo$TaxID <- as.character(krakenTaxInfo$TaxID)

amrExp <- lapply(amrAnalytical, function(x){
  amrMR <- newMRexperiment(x[rowSums(x) > 0, ])
  return(amrMR)
})
                 
krakenExp <- lapply(krakenAnalytical, function(x){
  krakenMR <- newMRexperiment(x[rowSums(x) > 0, ])
  return(krakenMR)
})
                 

amrNorm <- lapply(amrExp, function(x){
  amrCSS <- cumNorm(x)
})

amrNorm <- lapply(amrNorm, function(x){
  amrCSSdf <- data.frame(MRcounts(x, norm = T))
})


krakenNorm <- lapply(krakenExp, function(x){
  krakenCSS <- cumNorm(x)
})

krakenNorm <- lapply(krakenNorm, function(x){
  krakenCSSdf <- data.frame(MRcounts(x, norm = T))
})

amrAnnotations <- read_tsv('amr_genes.tabular_parsed.tab')

# amrQuarter <- amrResultsFiltered %>%
#   filter(Sample_type == "D0.25") %>%
#   select(Gene) %>%
#   distinct(Gene)
# 
# amrHalf <- amrResultsTidy %>%
#   filter(Sample_type == "D0.5" & Category == "Gene") %>%
#   select(CategoryName) %>%
#   distinct(CategoryName)
# 
# amrFull <- amrResultsTidy %>%
#   filter(Sample_type == "D1" & Category == "Gene") %>%
#   select(CategoryName) %>%
#   distinct(CategoryName)
# 

names(amrAnnotations) <- c("Gene", "Class", "Mechanism", "Group")

amrNormAnnot <- lapply(amrNorm, function(x){
  x$Gene <- row.names(x)
  amrAnnotated <- left_join(amrAnnotations, x, by="Gene") %>% 
    na.omit()
  return(amrAnnotated)
})

krakenNormAnnot <- lapply(krakenNorm, function(x){
  x$TaxID <- row.names(x)
  krakenAnnotated <- left_join(krakenTaxInfo, x, by="TaxID") %>% 
    na.omit() %>%
    select(-matches("Sample", "SampleType"))
  return(krakenAnnotated)
})


# Tidy normalized, annotated datasets -------------------------------------

amrNormTidy <- lapply(amrNormAnnot, function(x){
  x %>% 
    gather(key = samples, value = normCounts, 5:ncol(x)) %>%
    gather(key = category, value = categoryNames, 1:4)
})
  
krakenNormTidy <- lapply(krakenNormAnnot, function(x){
  x %>% 
    gather(key = samples, value = normCounts, 8:ncol(x))
})
  

# Split, apply, combine: aggregate AMR and Kraken -------------------------

# Join the data from all depths, then split by AMR category and taxonomy

amrAllDepths <- do.call("rbind", amrNormTidy)

amrAllDepths$SampleType <- str_extract(amrAllDepths$samples, "^[A-Z]+")

krakenAllDepths <- do.call("rbind", krakenNormTidy)

krakenAllDepths$SampleType <- str_extract(krakenAllDepths$samples, "^[A-Z]+")


amrNormAgg <- amrAllDepths %>% 
    split(.$category)
  
amrNormAgg <- lapply(amrNormAgg, function(x){
  group_by(x, categoryNames, samples) %>%
    summarise(normCountsSum = sum(normCounts))
})

amrNormDivMat <- lapply(amrNormAgg, function(x){
  amrNormWide <- spread(x, key=categoryNames, value=normCountsSum, fill=0)
  row.names(amrNormWide) <- amrNormWide$samples
  amrNormWide <- amrNormWide %>%
    select(2:ncol(amrNormWide))
  return(amrNormWide)
})

# Normalized diversity indices --------------------------------------------

amrDiversity <- lapply(amrNormDivMat, function(x){
  observed_richness <- specnumber(x, MARGIN=1)
  invsimpson <- diversity(x, index="invsimpson", MARGIN=1)
  simpson <- diversity(x, index="simpson", MARGIN=1)
  shannon <- diversity(x, index="shannon", MARGIN=1)
  evenness <- shannon/log(observed_richness)
  return(list(observed_richness=observed_richness,
              invsimpson = invsimpson,
              simpson = simpson,
              shannon = shannon,
              evenness=evenness))
})

amrEstimated <- lapply(amrNormDivMat, function(x){
  specpool(x, sampleMetadata$SampleType)
})

amrDiversityDF <- lapply(amrDiversity, function(x) data.frame(
  ID=names(x$observed_richness), 
  InvSimpson=as.numeric(x$invsimpson),
  Simpson = as.numeric(x$simpson),
  Shannon=as.numeric(x$shannon), 
  Evenness=as.numeric(x$evenness)
))

amrDiversityDF <- do.call("rbind", amrDiversityDF)

amrDiversityDF <- amrDiversityDF %>% 
  mutate(Level=row.names(amrDiversityDF)) %>%
  mutate(Level=str_extract(Level, "^\\w+")) %>%
  mutate(Depth=str_extract(ID, "^[A-Z]+")) %>%
  mutate(Depth=str_replace(Depth,"F", "D1")) %>%
  mutate(Depth=str_replace(Depth,"H", "D0.5")) %>%
  mutate(Depth=str_replace(Depth, "Q", "D0.25"))

amrDiversityDF$Level <- factor(amrDiversityDF$Level, 
                                  levels = c('Class',
                                             'Mechanism',
                                             'Group',
                                             'Gene'))

amrEstimatedDF <- do.call("rbind", amrEstimated)

amrEstimatedDF <- amrEstimatedDF %>% 
  mutate(amrLevel=row.names(amrEstimatedDF)) %>% 
  separate(amrLevel, 
           into = c("Level", "Depth"),
           sep="\\.")

krakenNormAgg <- krakenAllDepths %>% 
    split(.$TaxRank)
  
krakenNormAgg <- lapply(krakenNormAgg, function(x){
  group_by(x, TaxID, samples) %>%
    summarise(normCountsSum = sum(normCounts))
})

krakenNormDivMat <- lapply(krakenNormAgg, function(x){
  krakenNormWide <- spread(x, key=TaxID, value=normCountsSum, fill=0)
  row.names(krakenNormWide) <- krakenNormWide$samples
  krakenNormWide <- krakenNormWide %>%
    select(2:ncol(krakenNormWide))
  return(krakenNormWide)
})

krakenDiversity <- lapply(krakenNormDivMat, function(x){
  observed_richness <- specnumber(x, MARGIN=1)
  invsimpson <- diversity(x, index="invsimpson", MARGIN=1)
  simpson <- diversity(x, index="simpson", MARGIN=1)
  shannon <- diversity(x, index="shannon", MARGIN=1)
  evenness <- shannon/log(observed_richness)
  return(list(observed_richness=observed_richness,
              invsimpson = invsimpson,
              simpson = simpson,
              shannon = shannon,
              evenness=evenness))
})

krakenEstimated <- lapply(krakenNormDivMat, function(x){
  specpool(x, sampleMetadata$SampleType)
})


krakenDiversityDF <- lapply(krakenDiversity, function(x) data.frame(
  ID=names(x$observed_richness), 
  InvSimpson=as.numeric(x$invsimpson),
  Simpson = as.numeric(x$simpson),
  Shannon=as.numeric(x$shannon), 
  Evenness=as.numeric(x$evenness)
))

krakenDiversityDF <- do.call("rbind", krakenDiversityDF)

krakenDiversityDF <- krakenDiversityDF %>% 
  mutate(taxLevel=row.names(krakenDiversityDF)) %>%
  separate(taxLevel, into=c("Level", "sample_number"), sep="\\.") %>%
  select(-sample_number) %>%
  filter(!Level %in% c("-", "D")) %>%
  mutate(Depth = str_extract(ID, "^[A-Z]+"))
 
krakenDiversityDF$Level <- str_replace(krakenDiversityDF$Level, "P", "Phylum")
krakenDiversityDF$Level <- str_replace(krakenDiversityDF$Level, "C", "Class")
krakenDiversityDF$Level <- str_replace(krakenDiversityDF$Level, "O", "Order")
krakenDiversityDF$Level <- str_replace(krakenDiversityDF$Level, "F", "Family")
krakenDiversityDF$Level <- str_replace(krakenDiversityDF$Level, "G", "Genus")
krakenDiversityDF$Level <- str_replace(krakenDiversityDF$Level, "S", "Species")

krakenDiversityDF$Depth <- str_replace(krakenDiversityDF$Depth, "F", "D1")
krakenDiversityDF$Depth <- str_replace(krakenDiversityDF$Depth, "H", "D0.5")
krakenDiversityDF$Depth <- str_replace(krakenDiversityDF$Depth, "QD", "D0.25")

krakenDiversityDF$Level <- factor(krakenDiversityDF$Level, 
                                         levels = c('Phylum', 
                                                    'Class', 
                                                    'Order', 
                                                    'Family', 
                                                    'Genus', 
                                                    'Species'))
 
krakenEstimatedDF <- do.call("rbind", krakenEstimated)

krakenEstimatedDF <- krakenEstimatedDF %>% 
  mutate(krakenLevel=row.names(krakenEstimatedDF)) %>% 
  separate(krakenLevel,
           into = c("Level", "Depth"),
           remove = TRUE)

write_csv(krakenEstimatedDF, '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/kraken_estimated_richness.csv')

write_csv(amrEstimatedDF, '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/amr_estimated_richness.csv')

# Normalized boxplots -----------------------------------------------------

amrShannonBoxplot <- amrDiversityDF %>%
    amrShannon()

ggsave(filename = 'amrShannon.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/alphaDiversity/',
       plot = amrShannonBoxplot,
       width = 10.50,
       height = 8.50,
       units = "in")

amrInvSimpsonBoxplot <- amrDiversityDF %>% 
  amrAlphaDiv()

ggsave(filename = 'amrInvSimpson.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/alphaDiversity/',
       plot = amrInvSimpsonBoxplot,
       width = 10.50,
       height = 8.50,
       units = "in")

krakenShannonBoxplot <- krakenDiversityDF %>%
    krakenShannon()

ggsave(filename = 'krakenShannon.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/alphaDiversity',
       plot = krakenShannonBoxplot,
       width = 10.50,
       height = 8.50,
       units = "in")

krakenInvSimpsonBoxplot <- krakenDiversityDF %>%
  krakenAlphaDiv()

ggsave(filename = 'krakenInvSimpson.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/alphaDiversity',
       plot = krakenInvSimpsonBoxplot,
       width = 10.50,
       height = 8.50,
       units = "in")


# Construction of rarefaction curves --------------------------------------

# AMR curves built from the Coverage Sampler results

# This step requires optimization (parallel mapping or compiling, maybe?)

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

# purrr version

#amrRarefy <- map(amrResultsMat, function(x){
#  raremax <- min(rowSums(x))
#  rarecurve_ROP(x, step=5, sample=raremax)
#})

# mclapply version (and the winner!)

amrRarCurve <- mclapply(amrResultsMat, function(x){
  raremax <- min(rowSums(x))
  rarecurve(x, step=10000, sample=raremax)
},mc.cores=10)

krakenRarCurve <- mclapply(krakenResultsMat, function(x){
  raremax <- min(rowSums(x))
  rarecurve(x, step=1000, sample=raremax)
},mc.cores=10)

krakenRarCurveNoSample <- mclapply(krakenResultsMat, function(x){
  rarecurve(x, step=1000)
},mc.cores=10)

# Rename list of rarefied data using the sample names
# Need to think of a better way to extract the sample names

amrSamples <- rownames(amrResultsMat[["Class"]])

krakenSamples <- rownames(krakenResultsMat[["D"]])

# Review syntax here
# Might need to add function to utility function list

amrRarCurve <- mclapply(amrRarCurve, function(x){
  set_names(x,amrSamples)}, mc.cores=3)

krakenRarCurve <- mclapply(krakenRarCurve, function(x){
  set_names(x,krakenSamples)}, mc.cores=10)

# Generate list of dataframes
# Hacky, but it works. Need to make it cleaner and faster.

amrRarCurve <- amrRarCurve %>% 
  map(as_vector)

krakenRarCurve <- krakenRarCurve %>% 
  map(as_vector)

amrRarCurveDF <- map(amrRarCurve, function(x) as_tibble(x, attr(x, "names")))

krakenRarefyDF <- map(krakenRarCurve, function(x) as_tibble(x, attr(x, "names")))

# Split rownames in order to generate columns with useful information

amrRarCurveDF <- mclapply(amrRarCurveDF, function(x) {
  x$Sample <- row.names(x)
  x$Subsample <- str_extract(x$Sample, "N\\d+")
  x$Subsample <- as.numeric(str_replace(x$Subsample, "N",""))
  x$SampleID <- str_replace(x$Sample, "\\.N.*$", "")
  x$Depth <- str_replace(x$Sample, "_.*", "")
  x
}, 
mc.cores=3
)

krakenRarefyDF <- mclapply(krakenRarefyDF, function(x) {
  x$Sample <- row.names(x)
  x$Subsample <- str_extract(x$Sample, "\\.N\\d+")
  x$Subsample <- as.numeric(str_replace(x$Subsample, "\\.N",""))
  x$Depth <- str_replace(x$Sample, "_.*", "")
  x$SampleID <- str_extract(x$Sample, "_.*\\.")
  x$SampleID <- str_replace(x$SampleID, "_", "")
  x$SampleID <- str_replace(x$SampleID, "\\.", "")
  x$Sample <- str_replace(x$Sample, "\\.N.*$", "")
  x
}, 
mc.cores=10
)

amrRarCurveDF <- do.call("rbind", amrRarCurveDF)
amrRarCurveDF$AMRLevel <- row.names(amrRarCurveDF)
amrRarCurveDF$AMRLevel <- str_extract(amrRarCurveDF$AMRLevel, "^\\w+")
amrRarCurveDF$Depth <- str_replace(amrRarCurveDF$Depth, "F", "D1")
amrRarCurveDF$Depth <- str_replace(amrRarCurveDF$Depth, "H", "D0.5")
amrRarCurveDF$Depth <- str_replace(amrRarCurveDF$Depth, "QD", "D0.25")

krakenRarefyDF <- do.call("rbind", krakenRarefyDF)
krakenRarefyDF$krakenLevel <- row.names(krakenRarefyDF)
krakenRarefyDF$krakenLevel <- str_extract(krakenRarefyDF$krakenLevel, "^.\\.")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "\\.", "")
krakenRarefyDF$Depth <- str_replace(krakenRarefyDF$Depth, "F", "D1")
krakenRarefyDF$Depth <- str_replace(krakenRarefyDF$Depth, "H", "D0.5")
krakenRarefyDF$Depth <- str_replace(krakenRarefyDF$Depth, "QD", "D0.25")

amrRarCurveDF$AMRLevel <- factor(amrRarCurveDF$AMRLevel, 
                                  levels = c('Class', 'Mechanism', 'Group', 'Gene'))

amrAlphaRarefactionDF$Level <- factor(amrAlphaRarefactionDF$Level, 
                                 levels = c('Class', 'Mechanism', 'Group', 'Gene'))


# Need to apply more functional programming
 
krakenRarefyDF <- krakenRarefyDF %>% filter(!krakenLevel %in% c("-", "D"))

krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "P", "Phyla")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "C", "Classes")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "O", "Orders")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "F", "Families")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "G", "Genera")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "S", "Species")

amrRarCurveDF <- amrRarCurveDF %>%
  mutate(AMRLevel=str_replace(AMRLevel,"Class","Classes")) %>%
  mutate(AMRLevel=str_replace(AMRLevel,"Mechanism","Mechanisms")) %>%
  mutate(AMRLevel=str_replace(AMRLevel,"Group","Groups")) %>%
  mutate(AMRLevel=str_replace(AMRLevel,"Gene","Genes"))

# Split AMR data frame and generate rarefaction curves for each AMR level

amrAllList <- amrRarCurveDF %>% 
  split(.$AMRLevel)

amrAllRarCurves <- amrAllList %>%
  map(function(x){
    amrRarefactionCurve(x)
  })

# Split Kraken data frame and generate rarefaction curves for each taxonomic level

krakenAllRarCurveList <- krakenRarefyDF %>% 
  split(.$krakenLevel)

krakenAllRarCurvesLarger <- krakenAllRarCurveList %>%
  map(function(x){
    krakenRarefactionCurve(x)
  })

foreach(i=krakenAllRarCurvesLarger) %do%
  ggsave(filename=paste('rarefaction',unique(i$data$krakenLevel),'CB','noFacet','.png', sep='', collapse=''),
         path = '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/rarefaction',
         plot = i,
         height=8.50,
         width=10.50,
         units="in",
         device="png")

amrAllAlphaBoxPlots <- amrAlphaRarefactionDF %>%
  amrAlphaDiv()

ggsave('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/alphaDiversity/amrAlphaDiversityCB.png',
       plot = amrAllAlphaBoxPlots,
       width = 10.50,
       height = 8.50,
       units = "in")

# Alpha Diversity (rarefied) and Species Richness calculations -------------

krakenAlphaRarefaction <- mclapply(krakenResultsMat, function(x){
  alpha_rarefaction(x, minlevel=0)
},mc.cores = 10)

amrAlphaRarefaction <- mclapply(amrResultsMat, function(x){
  alpha_rarefaction(x, minlevel=0)
}, mc.cores=10)

# Consider change to data frame/tibble 

krakenAlphaRarefaction <- lapply(krakenAlphaRarefaction2, function(x) data.table(
  ID=names(x$raw_species_abundance), 
  RawSpeciesAbundance=as.numeric(x$raw_species_abundance), 
  RarSpeciesAbundance=as.numeric(x$rarefied_species_abundance), 
  AlphaDiv=as.numeric(x$alphadiv), 
  Shannon=as.numeric(x$shannon), 
  Evenness=as.numeric(x$evenness)
))

amrAlphaRarefaction <- lapply(amrAlphaRarefaction, function(x) data.table(
  ID=names(x$raw_species_abundance), 
  RawSpeciesAbundance=as.numeric(x$raw_species_abundance), 
  RarSpeciesAbundance=as.numeric(x$rarefied_species_abundance), 
  AlphaDiv=as.numeric(x$alphadiv), 
  Shannon=as.numeric(x$shannon), 
  Evenness=as.numeric(x$evenness)
))

amrCategories <- levels(amrResultsTidy$Category)

foreach(i=amrAlphaRarefaction, j=amrCategories) %do%
  rep(j,length(i$ID)) -> i$Level

#amrAlphaRarefaction[[1]]$Level <- rep(amrCategories[[1]], length(amrAlphaRarefaction[[1]]$ID))
#amrAlphaRarefaction[[2]]$Level <- rep(amrCategories[[2]], length(amrAlphaRarefaction[[2]]$ID))
#amrAlphaRarefaction[[3]]$Level <- rep(amrCategories[[3]], length(amrAlphaRarefaction[[3]]$ID))
#amrAlphaRarefaction[[4]]$Level <- rep(amrCategories[[4]], length(amrAlphaRarefaction[[4]]$ID))

foreach(i=krakenAlphaRarefaction, j=krakenTaxa) %do%
  rep(j,length(i$ID)) -> i$Level


foreach(i=krakenAlphaRarefaction2, j=krakenTaxa) %do%
  rep(j,length(i$ID)) -> i$Level

# Seems like unnecessary repetition to convert these to data frames when they
# could have been created as dataframes

amrAlphaRarefactionDF <- lapply(amrAlphaRarefaction, function(x){
  x <- as.data.frame(x)
  x
})

krakenAlphaRarefactionDF <- lapply(krakenAlphaRarefaction, function(x){
  x <- as.data.frame(x)
  x
})

krakenAlphaRarefaction2DF <- lapply(krakenAlphaRarefaction2, function(x){
  x <- as.data.frame(x)
  x
})

amrAlphaRarefactionDF <- do.call("rbind", amrAlphaRarefactionDF)
  
krakenAlphaRarefactionDF <- do.call("rbind", krakenAlphaRarefactionDF)

krakenAlphaRarefaction2DF <- do.call("rbind", krakenAlphaRarefaction2DF)

amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$ID, "_.*", "")

krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$ID, "_.*", "")


krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$ID, "_.*", "")
# Generate one single dataframe and create column for AMR Level

krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$Depth, "F", "D1")
krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$Depth, "H", "D0.5")
krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$Depth, "QD", "D0.25")

amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$Depth, "F", "D1")
amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$Depth, "H", "D0.5")
amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$Depth, "QD", "D0.25")

krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$Depth, "F", "D1")
krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$Depth, "H", "D0.5")
krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$Depth, "QD", "D0.25")


# Alpha Diversity (rarefied) and Species Richness boxplots -------------------

krakenAlphaRarefactionDF <- krakenAlphaRarefactionDF %>% filter(!Level %in% c("-", "D"))
krakenAlphaRarefaction2DF <- krakenAlphaRarefaction2DF %>% filter(!Level %in% c("-", "D"))

krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "P", "Phyla")
krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "C", "Classes")
krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "O", "Orders")
krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "F", "Families")
krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "G", "Genera")
krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "S", "Species")

krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "P", "Phyla")
krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "C", "Classes")
krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "O", "Orders")
krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "F", "Families")
krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "G", "Genera")
krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "S", "Species")

krakenAlphaRarefactionDF$Level <- factor(krakenAlphaRarefactionDF$Level, 
                                  levels = c('Phyla', 
                                             'Classes', 
                                             'Orders', 
                                             'Families', 
                                             'Genera', 
                                             'Species'))

amrAllSpRawBoxPlots <- amrAlphaRarefactionDF %>%
    amrRawSpeciesRich()

ggsave(filename = 'amrSpeciesRichnessCB.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/alphaDiversity',
       plot = amrAllSpRawBoxPlots,
       width = 10.50,
       height = 8.50,
       units = "in")

krakenAllAlphaBoxPlots <- krakenAlphaRarefaction2DF %>%
    krakenAlphaDiv()

ggsave(filename = 'krakenAlphaDivCB.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/alphaDiversity',
       plot = krakenAllAlphaBoxPlots,
       width = 10.50,
       height = 8.50,
       units = "in")

krakenAllShannonBoxPlots <- krakenAlphaRarefaction2DF %>%
  krakenShannon()

ggsave(filename = 'krakenAlphaDivCB.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/alphaDiversity',
       plot = krakenAllAlphaBoxPlots,
       width = 10.50,
       height = 8.50,
       units = "in")


krakenRarefiedNewNames <- krakenAlphaRarefaction2DF %>% 
  mutate(Level=str_replace(Level,"Phyla", "Phylum")) %>%
  mutate(Level=str_replace(Level,"Classes", "Class")) %>%
  mutate(Level = str_replace(Level, "Orders", "Order")) %>%
  mutate(Level=str_replace(Level, "Families", "Family")) %>%
  mutate(Level=str_replace(Level, "Genera", "Genus")) %>%
  mutate(Level=str_replace(Level, "Species", "Species"))

krakenRarefiedNewNames$Level <- factor(krakenRarefiedNewNames$Level, 
                                  levels = c('Phylum', 
                                             'Class', 
                                             'Order', 
                                             'Family', 
                                             'Genus', 
                                             'Species'))

krakenAllSpRawBoxPlots <- krakenRarefiedNewNames %>%
  krakenRawSpeciesRich()

ggsave(filename = 'krakenObservedRichness.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/alphaDiversity',
       plot = krakenAllSpRawBoxPlots,
       width = 10.50,
       height = 8.50,
       units="in")

krakenAllSpRawBoxPlots <- krakenAlphaRarefaction2DF %>%
  krakenRawSpeciesRich()

ggsave(filename = 'krakenSpRichnessCB.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/alphaDiversity',
       plot = krakenAllSpRawBoxPlots,
       width = 10.50,
       height = 8.50,
       units="in")

# Correlation plots -------------------------------------------------------

# Should make into a function for the package

amrReadstoHitRatio$Sample_type <- str_replace(amrReadstoHitRatio$Sample_type, "F", "D1")

amrReadstoHitRatio$Sample_type <- str_replace(amrReadstoHitRatio$Sample_type, "H", "D0.5")

amrReadstoHitRatio$Sample_type <- str_replace(amrReadstoHitRatio$Sample_type, "Q", "D0.25")

amrReadsvsHits <- cor(x = amrReadstoHitRatio$`Number of reads`, amrReadstoHitRatio$`AMR hits`, method = "spearman")

amrReadsvsHitsCorTest <- cor.test(x = amrReadstoHitRatio$`Number of reads`, amrReadstoHitRatio$`AMR hits`, method = "spearman")

amrReadsvsHitsCor <- ggplot(amrReadstoHitRatio, aes(`Number of reads`, `AMR hits`)) + 
  geom_point(aes(fill=Sample_type), 
             alpha=0.6, 
             size=10, 
             pch=21, 
             color="grey") 
  # geom_smooth(aes(group=1, weight=0.2), 
  #             method="lm", 
  #             se=FALSE, 
  #             colour="grey", 
  #             alpha=0.5) 

amrReadsvsHitsCor +
  ylab("Number of AMR Hits\n") +
  xlab("\nNumber of reads") +
  theme(axis.text.y=element_text(size=35),
        axis.title.y=element_text(size=44),
        axis.text.x=element_text(size=35),
        axis.title.x=element_text(size=44),
        legend.title=element_text(size=36),
        legend.text=element_text(size=36, vjust=0.5),
        legend.key = element_rect(size = 2),
        legend.key.size = unit(2, "lines"),
        legend.spacing = unit(0.2,"lines"),
        panel.background = element_rect(fill = "grey90", colour = "grey80")) +
  scale_fill_manual(values=rev(cbPalette),
                     name="Sequencing Depth\n")
ggsave('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/amrReadsvsHitsCorUpdated.png', 
       width=14, 
       height=8.50,
       units="in")

# Generating kraken correlation plot

krakenReadsvsHitsCorTest <- cor.test(x = amrReadstoHitRatio$`Number of reads`, amrReadstoHitRatio$Phylum_hits, method = "spearman")

krakenReadsvsHitsCor <- ggplot(amrReadstoHitRatio, aes(`Number of reads`, Phylum_hits)) + 
  geom_point(aes(fill=Sample_type), alpha=0.6, size=10, pch=21, color="grey") 
  # geom_smooth(aes(group=1, weight=0.2), method="lm", se=FALSE, colour="grey", alpha=0.5) 

krakenReadsvsHitsCor +
  ylab("Number of Kraken Hits\n") +
  xlab("\nNumber of reads") +
  theme(axis.text.y=element_text(size=35),
        axis.title.y=element_text(size=44),
        axis.text.x=element_text(size=35),
        axis.title.x=element_text(size=44),
        legend.title=element_text(size=36),
        legend.text=element_text(size=36, vjust=0.5),
        legend.key = element_rect(size = 2),
        legend.key.size = unit(2, "lines"),
        legend.spacing = unit(0.2,"lines"),
        panel.background = element_rect(fill = "grey90", colour = "grey80")) +
  scale_fill_manual(values=rev(cbPalette),
                     name="Sequencing Depth\n")

ggsave(filename='krakenReadsvsHitsCorUpdated.png', 
       path='~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017',
       width=14, 
       height=8.50,
       units="in")

# Kruskal-Wallis tests ----------------------------------------------------

# Resistome

amrAlphaRarefactionLevels <- amrAlphaRarefactionDF %>%
  split(.$Level)

amrRichnessKruskal <- amrAlphaRarefactionLevels %>%
  map(function(x){
    kruskal.test(RawSpeciesAbundance ~ Depth, data = x)
    })

amrRichnessPosthoc <- amrAlphaRarefactionLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(RawSpeciesAbundance ~ Depth, data=x, dist="Chisq")
  })

amrAlphaKruskal <- amrAlphaRarefactionLevels %>%
  map(function(x){
    kruskal.test(AlphaDiv ~ Depth, data = x)
  })

amrAlphaPosthoc <- amrAlphaRarefactionLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(AlphaDiv ~ Depth, data = x, dist="Chisq")
  })

# Microbiome

krakenAlphaRarefaction2DF$Depth <- as.factor(krakenAlphaRarefaction2DF$Depth)

krakenAlphaRarefactionLevels <- krakenAlphaRarefaction2DF %>%
  split(.$Level)

krakenRichnessKruskal <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    kruskal.test(RawSpeccolor=Sample_typeiesAbundance ~ Depth, data = x)
    })

krakenRichnessPosthoc <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(RawSpeciesAbundance ~ Depth, data=x, dist="Chisq")
  })

krakenAlphaKruskal <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    kruskal.test(AlphaDiv ~ Depth, data = x)
  })

krakenAlphaPosthoc <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(AlphaDiv ~ Depth, data = x, dist="Chisq")
  })


# Kruskal-Wallis tests of normalized diversity data -----------------------

# Resistome

amrDiversityLevels <- amrDiversityDF %>%
  split(.$Level)

amrInvSimpsonKruskall <- amrDiversityLevels %>%
  map(function(x){
    kruskal.test(InvSimpson ~ as.factor(Depth), data = x)
    })

amrISKruskallTidy <- lapply(amrInvSimpsonKruskall, function(x){
  tidy(x)
})

amrISKruskallTidy <- do.call("rbind", amrISKruskallTidy)

write.csv(amrISKruskallTidy, '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/amrInvSimpsonKruskalWallis_dos.csv', row.names=TRUE)

amrShannonKruskall <- amrDiversityLevels %>%
  map(function(x){
    kruskal.test(Shannon ~ as.factor(Depth), data = x)
  })

amrShannonKruskallTidy <- lapply(amrShannonKruskall, function(x){
  tidy(x)
})

amrShannonKruskallTidy <- do.call("rbind", amrShannonKruskallTidy)

write.csv(amrShannonKruskallTidy, '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/amrShannonKruskalWallis_dos.csv', row.names=TRUE)

amrShannonPostHoc <- amrDiversityLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(Shannon ~ as.factor(Depth), data=x, dist="Chisq")
  })


amrAlphaKruskal <- amrAlphaRarefactionLevels %>%
  map(function(x){
    kruskal.test(AlphaDiv ~ Depth, data = x)
  })

amrAlphaPosthoc <- amrAlphaRarefactionLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(AlphaDiv ~ Depth, data = x, dist="Chisq")
  })

# Microbiome

krakenAlphaRarefaction2DF$Depth <- as.factor(krakenAlphaRarefaction2DF$Depth)

krakenAlphaRarefactionLevels <- krakenAlphaRarefaction2DF %>%
  split(.$Level)

krakenRichnessKruskal <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    kruskal.test(RawSpeccolor=Sample_typeiesAbundance ~ Depth, data = x)
    })

krakenRichnessPosthoc <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(RawSpeciesAbundance ~ Depth, data=x, dist="Chisq")
  })

krakenAlphaKruskal <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    kruskal.test(AlphaDiv ~ Depth, data = x)
  })

krakenAlphaPosthoc <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(AlphaDiv ~ Depth, data = x, dist="Chisq")
  })

# Normalized

krakenDiversityLevels <- krakenDiversityDF %>%
  split(.$Level)

krakenInvSimpsonKruskall <- krakenDiversityLevels %>%
  map(function(x){
    kruskal.test(InvSimpson ~ as.factor(Depth), data = x)
    })

krakenISKruskallTidy <- lapply(krakenInvSimpsonKruskall, function(x){
  tidy(x)
})

krakenISKruskallTidy <- do.call("rbind", krakenISKruskallTidy)

write.csv(krakenISKruskallTidy, '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/krakenInvSimpsonKruskalWallis_dos.csv', row.names = TRUE)

krakenShannonKruskall <- krakenDiversityLevels %>%
  map(function(x){
    kruskal.test(Shannon ~ as.factor(Depth), data = x)
  })

krakenShannonKrusnkallTidy <- lapply(krakenShannonKruskall, function(x){
  tidy(x)
})

krakenShannonKruskallTidy <- do.call("rbind", krakenShannonKruskallTidy)

write.csv(krakenShannonKruskallTidy, '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/krakenShannonKruskalWallis_dos.csv', row.names = TRUE)

krakenShannonPostHoc <- krakenDiversityLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(Shannon ~ as.factor(Depth), data=x, dist="Chisq")
  })

# AMR Rarefaction Curves (Rarefaction Analyzer) ---------------------------

# This work was done using the percent sampled vs. number of observations 


amrRarConcat <- read_csv('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/rarefiedConcat.csv')

amrRarConcatFilter <- amrRarConcat %>%
  mutate(Depth=str_replace(Depth,"full", "D1")) %>%
  mutate(Depth=str_replace(Depth,"half[1-2]", "D0.5")) %>%
  mutate(Depth=str_replace(Depth, "quarter", "D0.25")) %>%
  filter(Depth != "seqtk") %>%
  mutate(Level=str_replace(Level,"class","Classes")) %>%
  mutate(Level=str_replace(Level,"mechanism","Mechanisms")) %>%
  mutate(Level=str_replace(Level, "group", "Groups")) %>%
  mutate(Level=str_replace(Level, "gene", "Genes"))

amrRarConcatbyLevel <- amrRarConcatFilter %>%
  split(.$Level)

amrRarCurvesPercent <- amrRarConcatbyLevel %>%
  map(function(x){
    amrRarFigure(x)
  })


# AMR scaled rarefaction curves -------------------------------------------

# Will also use data from rarefaction analyzer, but the difference will be that the percent sampled will be scaled to the number of reads.

# Mutate Depth colum so that it matches the Sample_type column of the amrReadsto
# Hit Ratio dataframe

#> unique(amrReadstoHitRatio$Sample_type)
#[1] "F" "H" "Q"

amrRarConcat <- read_csv('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/rarefiedConcat.csv')

amrRarConcatFilter <- amrRarConcat %>%
  mutate(Depth=str_replace(Depth,"full", "F")) %>%
  mutate(Depth=str_replace(Depth,"half[1-2]", "H")) %>%
  mutate(Depth=str_replace(Depth, "quarter", "Q")) %>%
  filter(Depth != "seqtk") %>%
  rename(Sample_type = Depth) %>%
  mutate(Level=str_replace(Level,"class","Classes")) %>%
  mutate(Level=str_replace(Level,"mechanism","Mechanisms")) %>%
  mutate(Level=str_replace(Level, "group", "Groups")) %>%
  mutate(Level=str_replace(Level, "gene", "Genes"))

# Extract sample ID from the AMR reads dataframe

amrReadstoHitRatio <- amrReadstoHitRatio %>%
  mutate(SampleID=str_extract(Samples, "^\\w\\d+"))

# Merge (left join) AMR total reads dataframe with the concatenated AMR rarefacti
# on dataframe of the the sample ID.

amrScaledRar <- left_join(amrReadstoHitRatio, amrRarConcatFilter, by=c("SampleID","Sample_type"))

amrScaledRarCurve <- amrScaledRar %>% 
  mutate(Sample_type=str_replace(Sample_type,"F", "D1")) %>%
  mutate(Sample_type=str_replace(Sample_type,"H", "D0.5")) %>%
  mutate(Sample_type=str_replace(Sample_type, "Q", "D0.25"))

# Scale reads by multiplying the percent sampling by the total number of reads 
# for each sample

amrScaledRarCurve <- amrScaledRarCurve %>%
  mutate(PercentSampling=0.01*PercentSampling) %>%
  mutate(ScaledReads=`Number of reads`*PercentSampling) %>%
  rename(Depth = Sample_type)

amrScaledbyLevel <- amrScaledRarCurve %>%
  split(.$Level)

amrRarCurvesLines <- amrScaledbyLevel %>%
  map(function(x){
    amrScaledNonSmooth(x)
})


ggsave(filename = 'amrGeneRarefaction.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/rarefaction',
       plot = amrRarCurvesLines$Genes,
       width = 10.50,
       height = 8.50,
       units = "in")

ggsave(filename = 'amrClassRarefaction.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/rarefaction',
       plot = amrRarCurvesLines$Class,
       width = 10.50,
       height = 8.50,
       units = "in")

ggsave(filename = 'amrMechRarefaction.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/rarefaction',
       plot = amrRarCurvesLines$Mechanisms,
       width = 10.50,
       height = 8.50,
       units = "in")

ggsave(filename = 'amrGroupRarefaction.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/rarefaction',
       plot = amrRarCurvesLines$Groups,
       width = 10.50,
       height = 8.50,
       units = "in")


# Updated AMR rarefaction curves (AAFC annotation) ------------------------

amrRarConcatAAFC <- read_csv('~/amr/2-4-8_results/2_4_8_study_RZ/rarefaction_AAFC_annotation/rarefiedConcat_AAFC.csv')

amrRarConcatAAFCFilter <- amrRarConcatAAFC %>%
  mutate(Depth=str_replace(Depth,"F", "D1")) %>%
  mutate(Depth=str_replace(Depth,"H", "D0.5")) %>%
  mutate(Depth=str_replace(Depth, "Q", "D0.25")) %>%
  mutate(Level=str_replace(Level,"class","Classes")) %>%
  mutate(Level=str_replace(Level,"mechanism","Mechanisms")) %>%
  mutate(Level=str_replace(Level, "group", "Groups")) %>%
  mutate(Level=str_replace(Level, "gene", "Genes")) %>%
  rename(Samples = SampleID) %>% 
  rename(SampleID = Sample)

# > head(amrReadstoHitRatio)
# A tibble: 6 x 4
# Samples Depth `Number of reads` `AMR hits`
# <chr>   <chr>             <int>      <int>
# 1 A062    F              96153609     256108
# 2 A070    F             137157234     335299
# 3 N003    F              91882251     230886
# 4 N013    F              87784587     222068
# 5 S034    F             119514751     315742
# 6 S040    F              76537181     220583

# Extract sample ID from the AMR reads dataframe

#amrReadstoHitRatio <- amrReadstoHitRatio %>%
#  mutate(SampleID=str_extract(Samples, "^\\w\\d+"))

# Merge (left join) AMR total reads dataframe with the concatenated AMR rarefacti
# on dataframe of the the sample ID.

amrReadstoHitRatio <- amrReadstoHitRatio %>%
  rename(Depth = Sample_type) %>%
  mutate(Depth=str_replace(Depth,"F", "D1")) %>%
  mutate(Depth=str_replace(Depth,"H", "D0.5")) %>%
  mutate(Depth=str_replace(Depth, "Q", "D0.25"))

amrScaledAAFCRar <- left_join(amrReadstoHitRatio, amrRarConcatAAFCFilter, by=c("SampleID","Depth"))


# Scale reads by multiplying the percent sampling by the total number of reads 
# for each sample

amrAAFCScaledRarCurve <- amrScaledAAFCRar %>%
  mutate(PercentSampling=0.01*PercentSampling) %>%
  mutate(ScaledReads=`Number of reads`*PercentSampling)

amrAAFCScaledbyLevel <- amrAAFCScaledRarCurve %>%
  split(.$Level)

amrAAFCRarCurvesLines <- amrAAFCScaledbyLevel %>%
  map(function(x){
    amrScaledNonSmooth(x)
})

ggsave(filename = 'amrGeneRarefaction.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/rarefaction_AAFC_annotation/',
       plot = amrAAFCRarCurvesLines$Genes,
       width = 10.50,
       height = 8.50,
       units = "in")

ggsave(filename = 'amrClassRarefaction.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/rarefaction_AAFC_annotation/',
       plot = amrAAFCRarCurvesLines$Class,
       width = 10.50,
       height = 8.50,
       units = "in")

ggsave(filename = 'amrMechRarefaction.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/rarefaction_AAFC_annotation/',
       plot = amrAAFCRarCurvesLines$Mechanisms,
       width = 10.50,
       height = 8.50,
       units = "in")

ggsave(filename = 'amrGroupRarefaction.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/rarefaction_AAFC_annotation/',
       plot = amrAAFCRarCurvesLines$Groups,
       width = 10.50,
       height = 8.50,
       units = "in")


# Rarefaction curves of Kraken rarefaction (Waffles) ----------------------

# Microbiome rarefaction data generated in NML cluster (Waffles)

krakenRarWaffles <- Sys.glob("~/amr/2-4-8_results/2_4_8_study_RZ/Kraken_rarefaction_waffles/RF_*")

krakenRarWafflesFileNames <- str_extract(krakenRarWaffles, "RF_.*$")

rarCategories <- c("Rates", 
                   "Reads", 
                   "Domains", 
                   "Phyla", 
                   "Classes", 
                   "Orders", 
                   "Families", 
                   "Genera", 
                   "Species")

krakenRarWafflesAll <- map(krakenRarWaffles, function(x){
  df <- read.csv(x, header=FALSE)
  df <- t(df)
  colnames(df) = rarCategories
  df <- df[2:nrow(df),]
  as.data.frame(df)
}) %>%
  set_names(nm=krakenRarWafflesFileNames)

krakenRarWafflesDF <- bind_rows(krakenRarWafflesAll, .id = "Sample")

krakenRarWafflesDF$Rates <- as.numeric(as.character(krakenRarWafflesDF$Rates))

krakenRarWafflesDF$Reads <- as.numeric(as.character(krakenRarWafflesDF$Reads))

krakenRarWafflesDF$Domains <- as.numeric(as.character(krakenRarWafflesDF$Domains))

krakenRarWafflesDF$Phyla <- as.numeric(as.character(krakenRarWafflesDF$Phyla))

krakenRarWafflesDF$Classes <- as.numeric(as.character(krakenRarWafflesDF$Classes))

krakenRarWafflesDF$Orders <- as.numeric(as.character(krakenRarWafflesDF$Orders))

krakenRarWafflesDF$Families <- as.numeric(as.character(krakenRarWafflesDF$Families))

krakenRarWafflesDF$Genera <- as.numeric(as.character(krakenRarWafflesDF$Genera))

krakenRarWafflesDF$Species <- as.numeric(as.character(krakenRarWafflesDF$Species))

#krakenRarWafflesDF$Sample <- as.factor(krakenRarWafflesDF$Sample)

krakenRarWafflesTidy <- krakenRarWafflesDF %>%
  tidyr::gather(key = Level, value = taxonCount, -Sample,-Rates,-Reads)

krakenRarWafflesTidy <- krakenRarWafflesTidy %>%
  mutate(Depth=Sample) %>%
  mutate(Depth=str_replace(Depth,"RF_F.*", "D1")) %>%
  mutate(Depth=str_replace(Depth,"RF_H.*", "D0.5")) %>%
  mutate(Depth=str_replace(Depth,"RF_Q.*", "D0.25"))

krakenRarWafflesbyTaxon <- krakenRarWafflesTidy %>% 
  split(.$Level)

krakenRarWafflesCurves <- krakenRarWafflesbyTaxon %>%
  map(function(x){
    krakenRarCurveNML(x)
  })

krakenRarCurvePhyla <- ggplot(krakenPhylaCombined, aes(Reads, taxonCount, color=Depth)) +
    geom_line(aes(group=Sample), alpha=0.5, size=4) + 
    #geom_point(aes(group=Sample), size=4) +
    theme(strip.text.x=element_text(size=35),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=42),
          axis.title.y=element_text(size=42),
          legend.position="right",
          legend.title=element_text(size=34),
          legend.text=element_text(size=34, vjust=0.5),
          legend.key.size = unit(2, "lines"),
          legend.spacing = unit(0.2,"lines"),
          plot.title=element_text(size=52, hjust=0.5)) +
    xlab("\nNumber of reads") +
    ylab(paste('Number ', 'of ', krakenRarWafflesbyTaxon$Phyla$Level, '\n')) +
    scale_color_manual(values=rev(cbPalette)) +
    ylim(0,40) +
    #scale_y_log10() +
    scale_x_continuous(labels=scientific) 
  #  scale_y_log10(labels=scientific)
  #facet_grid(. ~ Depth)

ggsave(filename = 'amrGroupRarefaction.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/rarefaction_AAFC_annotation/',
       plot = amrAAFCRarCurvesLines$Groups,
       width = 10.50,
       height = 8.50,
       units = "in")


# Relative abundance plots for review -------------------------------------

amrRelAbund <- read_csv('amrGeneRelAbundance.csv')

amrRelAbund <- tidyr::gather(amrRelAbund, key="Depth", value="Assignments", -Gene)

amrRelAbundPerc <- amrRelAbund %>%
  group_by(Depth) %>%
  mutate(Depth_Total=sum(Assignments)) %>%
  mutate(Rel_Abund=(Assignments/Depth_Total)*100)

amrMissing <- amrRelAbundPerc %>%
  filter(Depth == "D0.25" & Assignments == 0.0) %>%
  distinct(Gene)

amrMissingCompare <- amrRelAbundPerc %>%
  filter(Gene %in% amrMissing$Gene) %>%
  arrange(Gene)

amrMissingFilter <- amrMissingCompare %>%
  filter(Depth != "D0.25")

amrAnnotations$Class <- str_replace(amrAnnotations$Class, "betalactams", "Betalactams")

amrAnnotations$Class <- str_replace(amrAnnotations$Class, "Cationic_antimicrobial_peptides", "Cationic antimicrobial peptides")

amrAnnotations$Class <- str_replace(amrAnnotations$Class, "Multi-drug_resistance", "Multi-drug resistance")

amrMissingAnnotated <- left_join(amrMissingFilter, amrAnnotations, by="Gene") %>%
  arrange(Class)

# Hack to reorder Gene variable by Classes

amrMissingAnnotated$Gene = factor(amrMissingAnnotated$Gene, levels=unique(amrMissingAnnotated$Gene[order(amrMissingAnnotated$Class)]), ordered=TRUE)

newCBScale = dichromat::colorschemes$Categorical.12

amrRelAbundPlot <- ggplot(amrMissingAnnotated, 
                              aes(x=Gene, y=Rel_Abund,fill=Class)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=newCBScale) +
  facet_grid(. ~ Depth) +
  ylim(0.00,0.04) +
  xlab("Gene") +
  ylab("Relative Abundance (%)\n") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x=element_text(size=21),
        axis.text.y=element_text(size=24),
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24),
        legend.position="right",
        legend.title=element_text(size=15),
        legend.text=element_text(size=12, vjust=0.5))
  
ggsave(filename = 'amrMissingGenes.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/reviewer_comments/',
       plot = amrRelAbundPlot,
       width = 10.50,
       height = 8.50,
       units = "in")
