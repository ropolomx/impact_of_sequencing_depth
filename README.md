# Code for the analysis of the paper "Impact of sequencing depth on the characterization of the microbiome and resistome"

The code in this repository was written and used for analyzing the data in the study by [Zaheer et al. (2018).](https://www.nature.com/articles/s41598-018-24280-8)
These R and Python scripts were used for processing the output of the AMRPlusPlus-Kraken Galaxy workflows.

The version of the code that was used for the analysis submitted to the journal can be found in the [releases](https://github.com/ropolomx/impact_of_sequencing_depth/releases) page under the [__v1.0__ tag](https://github.com/ropolomx/impact_of_sequencing_depth/releases/tag/v1.0). 

A `refactoring` branch is currently being developed to polish the code and to improve reusability and for making it more maintainable for the future.

## Paper citation

```bibtex

@article{zaheer_impact_2018,
	title = {Impact of sequencing depth on the characterization of the microbiome and resistome},
	volume = {8},
	issn = {2045-2322},
	url = {https://www.nature.com/articles/s41598-018-24280-8},
	doi = {10.1038/s41598-018-24280-8},
	abstract = {
		Developments in high-throughput next generation sequencing (NGS) technology
		have rapidly advanced the understanding of overall microbial ecology as well as
		occurrence and diversity of specific genes within diverse environments. In the
		present study, we compared the ability of varying sequencing depths to generate
		meaningful information about the taxonomic structure and prevalence of
		antimicrobial resistance genes (ARGs) in the bovine fecal microbial community.
		Metagenomic sequencing was conducted on eight composite fecal samples
		originating from four beef cattle feedlots. Metagenomic DNA was sequenced to
		various depths, D1, D0.5 and D0.25, with average sample read counts of 117, 59
		and 26 million, respectively. A comparative analysis of the relative abundance
		of reads aligning to different phyla and antimicrobial classes indicated that
		the relative proportions of read assignments remained fairly constant
		regardless of depth. However, the number of reads being assigned to ARGs as
		well as to microbial taxa increased significantly with increasing depth. We
		found a depth of D0.5 was suitable to describe the microbiome and resistome of
		cattle fecal samples. This study helps define a balance between cost and
		required sequencing depth to acquire meaningful results.
		},
	language = {en},
	journal = {Scientific Reports},
	author = {Zaheer, Rahat and Noyes, Noelle and Ortega Polo, Rodrigo and Cook, Shaun R. and Marinier, Eric and Van Domselaar, Gary and Belk, Keith E. and Morley, Paul S. and McAllister, Tim A.},
	month = apr,
	year = {2018},
	pages = {5890}
}
```

## Software requirements

### R packages

_tidyverse_ packages:

`readr`

`tidyr`

`dplyr`

`stringr`

`purrr`

`ggplot2`


