library(readr)
library(stringr)
library(purrr)


# Generating AMR rarefaction curves from Rarefaction Analyzer results

# Reading AMR results that were generated with Rarefaction Analyzer (percent-based)

# Update path accordingly

amrRarefiedFiles <- Sys.glob(file.path("~",
                                       "amr",
                                       "2-4-8_results",
                                       "2_4_8_study_RZ",
                                       "amrResults_Aug2017_75_gene_frac",
                                       "*_results",
                                       "*rarefied*.tabular"))

# Split the path names and extract the sample name only

amrResultsNames <- str_split(amrRarefiedFiles, pattern = "\\/")

# Extraction of sample name is being done by extracting the 9th element out of each list 
# element.

amrResultsNames <- amrResultsNames %>% 
  map(function(x){
    filename <- x[9]
  })

# Split filename strings by dot character. Exclude extension.

amrResultsNames <- amrResultsNames %>%
  map(function(x){
    rarSampleName <- str_split(x, pattern="\\.")
    rarSampleName
  })

# Remove one level of hierarchy to list

amrResultsNames <- flatten(amrResultsNames)

# Extract only first element (sample name and AMR level)

amrResultsNames <- amrResultsNames %>%
  map(function(x){
    x <- x[1]
    x
  })

amrResultsNames <- unlist(amrResultsNames)

# Let's now read all the Coverage Sampler tabular files
# We are using the readr package (read_tsv)
# We are also using the list of sample names extracted in the previous function
# to set the names of the list elements
# This will make life so much easier!

amrResults <- amrRarefiedFiles %>%
  map(read_tsv, col_names = c("Sample", "Counts")) %>%
  set_names(nm=amrResultsNames)

# Join all the datasets into one dataframe that will be analyzed

amrResults <- do.call("rbind", amrResults)

# Need to add a column with sample name to the Coverage Sampler output
# Also need to add a column with sample depth

amrResults$SampleName <- row.names(amrResults)
amrResults$SampleID <- str_extract(amrResults$SampleName, "^.*_rarefied")
amrResults$SampleID <- str_replace(amrResults$SampleID, "_rarefied", "")
amrResults$Depth <- str_extract(amrResults$SampleName, "rarefied_.*_")
amrResults$Depth <- str_replace(amrResults$Depth, "rarefied_", "")
amrResults$Depth <- str_replace(amrResults$Depth, "_","")
amrResults$amrLevel <- str_extract(amrResults$SampleName, "(class|mechanism|group|gene)")

# The other (easier) option is to read csv file containing all rarefied datasets which were concatenated with Python-Pandas

amrRarefiedConcat <- read_csv('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/rarefiedConcat.csv')

# The dataframes amrResults or amrRarefiedConcat that were generated above can be used for 
#plotting rarefaction curves with ggplot2
