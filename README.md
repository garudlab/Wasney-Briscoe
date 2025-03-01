# Uniform bacterial genetic diversity along the guts of mice inoculated with human stool
*Michael Wasney, Leah Briscoe, Richard Wolff, Hans Ghezzi, Carolina Tropini, Nandita Garud*

## Overview
This is the GitHub repository for Wasney and Briscoe et al., 2024, BioRxiv. This repository contains scripts used in all analyses performed in this study, as well as instructions on how to recreate the figures, tables, and statistics in our study. 

## Table of Contents

1. [Analysis](https://github.com/garudlab/Wasney-Briscoe/analysis): contains instructions in `.md` format on how to run the following pipelines/analyses
2. [Scripts](https://github.com/garudlab/Wasney-Briscoe/scripts): scripts used in all pipelines/analyses. Instructions on how to execute the scripts can be found in Analysis.
3. [Metadata](https://github.com/garudlab/Wasney-Briscoe/metadata): metadata files necessary to run analysis and generate figures.
4. [example_data](https://github.com/garudlab/Wasney-Briscoe/example_data): example MIDAS output for _Bacteroides vulgatus_, which can be used to demonstrate pipelines included in this repository.

## How to get started

### 1. Download data

*Coming soon.*

### 2. Clone this repository

This repository contains all code necessary to analyze the raw data and generate figures/tables. This process typicaly only takes a couple of seconds.

### 3. Run MIDAS

Ensure you are running MIDAS v1.2.2. To download MIDAS v1.2.2, consult the official [MIDAS v1 documentation](https://github.com/snayfach/MIDAS).

Consult [MIDAS.md](https://github.com/garudlab/Wasney-Briscoe/analysis/MIDAS.md) for instructions on running MIDAS on the mouse data. Note that running MIDAS can require significant computational resources, and we therefore recommend that this step is performed on a high performance computing cluster.

### 4. Run post-processing on data

Consult [postprocessing.md](https://github.com/garudlab/Wasney-Briscoe/analysis/postprocessing.md)

## Demonstrate

The following MIDAS outputs have been included in the [example_data](https://github.com/garudlab/Wasney-Briscoe/example_data) directory:
- [MIDAS species outputs](https://github.com/garudlab/Wasney-Briscoe/example_data/merged_data/species/): this includes the MIDAS species output corresponding to all species.
- [MIDAS genes outputs](https://github.com/garudlab/Wasney-Briscoe/example_data/merged_data/genes/): this includes the MIDAS genes output corresponding to _B. vulgatus_.
- [MIDAS genes outputs](https://github.com/garudlab/Wasney-Briscoe/example_data/merged_data/snps/): this includes the MIDAS snps output corresponding to _B. vulgatus_.
