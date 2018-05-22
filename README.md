# Code for the analysis of the paper "Impact of sequencing depth on the characterization of the microbiome and resistome"

The code in this repository was written and used for analyzing the data in [Zaheer et al. (2018). Impact of sequencing depth on the characterization of the microbiome and resistome. _Sci. Rep._ 8:5890](https://www.nature.com/articles/s41598-018-24280-8).
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


