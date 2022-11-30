# WHO_catalogue_paper

## Conversion from WHO catalogue to HGVS
Pull the entire repository and run the python script to translate the WHO catalogue to HGVS notation

## Exitrif analysis
The `./data/exitrif_all_variants.csv` file contains a row for each exitrif sample and a column for all the variants in the the exitrif dataset that lie in the tier 1 and 2 genes as well as the variants listed in the WHO catalogue. 
The value in a cell indicates whether a variant is present in a sample with where `1` indicates a resistant variant, `2` indicates a sensitive variant, `3` indicates a variant not classified in the catalogue that is sysnonymous, and `4` indicates an unclassified non sysnonynous variant.
