# polymorphology2

## Overview
`polymorphology2` is an R package that provides an assortment of wrapper functions for working with genomic data. It enables the analysis of polymorphisms in relation to genome features such as gene bodies, SBS profiles of mutations, and epigenome enrichments in genome features. 

This package is particularly useful for finding overlaps between features or between features and sites, such as calculating the number of sites across all features. It also includes various functions used by our lab, such as filtering somatic mutations called by strelka2, creating windows around genome features (like genes), and calculating the enrichment of ChIPseq experiments in genome features.

## Installation
Install the package directly from GitHub using `devtools` with the following command:
  
  ```r
devtools::install_github("greymonroe/polymorphology2")
```

## Basic Architecture
The package mainly works with two types of objects:
  
  1. **Features** - These are `data.table` objects with CHROM, START, STOP, and ID columns.
2. **Sites** - These are `data.table` objects with CHROM, POS, and ID columns.

Both sites and features can contain other columns which can be used for various calculations (e.g. calculating the average depth of ChIP results).

## Future Work
In the future, look for:
  
  - Parsing VCF files with specified formats (e.g. variants called with DeepVariant, strelka2, pbsv, HaplotypeCaller) and with additional functions for VCFs with multiple samples.
- Addition of tutorials for better understanding of the package usage and functionalities.

## Dependencies
The package depends on the following R packages:
  
  - data.table
- ggplot2
- seqinr
- vcfR
- stringr

## License
This package is free to use by anyone. You can repurpose it and do whatever you want with it.

## Contact
If you have any issues or questions, feel free to contact the maintainer:
  
  - Grey Monroe - <greymonroe@gmail.com>
