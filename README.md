# Uniform bacterial genetic diversity along the gut
*Michael Wasney, Leah Briscoe, Richard Wolff, Hans Ghezzi, Carolina Tropini, Nandita Garud*

## Overview
This is the GitHub repository for Wasney and Briscoe et al., 2024, BioRxiv. This repository contains scripts used in all analyses performed in this study, as well as instructions on how to recreate the figures, tables, and statistics in our study. 

## Table of Contents

1. [Analysis](https://github.com/garudlab/Wasney-Briscoe/tree/main/analysis): contains instructions in `.md` format on how to run the following pipelines/analyses
2. [Scripts](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts): scripts used in all pipelines/analyses. Instructions on how to execute the scripts can be found in Analysis.
3. [Metadata](https://github.com/garudlab/Wasney-Briscoe/tree/main/metadata): metadata files necessary to run analysis and generate figures.
4. [example_data](https://github.com/garudlab/Wasney-Briscoe/tree/main/example_data): example MIDAS output for _Bacteroides vulgatus_, which can be used to demonstrate pipelines included in this repository.
5. [environments](https://github.com/garudlab/Wasney-Briscoe/tree/main/environments): this folder contains .yml files specifying two conda environments used throughout this pipeline

## How to get started

### 1. Download data

*Coming soon.*

### 2. Clone this repository

This repository contains all code necessary to analyze the raw data and generate figures/tables. This process typicaly only takes a couple of seconds.

### 3. Create necessary conda environments

Code in this repository is built either for a python 2.7- or python 3.6-based environment (`python27_env` and `python_env`, respectively). Create these environments in your analysis workspace like so:

```
conda env create -f ~/Wasney-Briscoe/environments/python27_env.yml
conda env create -f ~/Wasney-Briscoe/environments/python_env.yml
```

Alternatively, you can use your own python 2.7- and python 3.6-based environments, provided they're loaded with the necessary packages.

### 4. Run MIDAS

Ensure you are running MIDAS v1.2.2. To download MIDAS v1.2.2, consult the official [MIDAS v1 documentation](https://github.com/snayfach/MIDAS).

Consult [MIDAS.md](https://github.com/garudlab/Wasney-Briscoe/tree/main/analysis/MIDAS.md) for instructions on running MIDAS on the mouse data. Note that running MIDAS can require significant computational resources, and we therefore recommend that this step is performed on a high performance computing cluster.

### 5. Run post-processing on data

Consult [postprocessing.md](https://github.com/garudlab/Wasney-Briscoe/tree/main/analysis/postprocessing.md)

## Demonstration with example data (_B. vulgatus_)

The following MIDAS outputs have been included in the [example_data](https://github.com/garudlab/Wasney-Briscoe/tree/main/example_data) directory:
- [MIDAS species outputs](https://github.com/garudlab/Wasney-Briscoe/tree/main/example_data/merged_data/species/): this includes the MIDAS species output corresponding to all species.
- [MIDAS genes outputs](https://github.com/garudlab/Wasney-Briscoe/tree/main/example_data/merged_data/genes/Bacteroides_vulgatus_57955): this includes the MIDAS genes output corresponding to _B. vulgatus_.
- [MIDAS genes outputs](https://github.com/garudlab/Wasney-Briscoe/tree/main/example_data/merged_data/snps/Bacteroides_vulgatus_57955): this includes the MIDAS snps output corresponding to _B. vulgatus_.

For demonstrating post-processing, $\pi$, and strain phasing pipeline and plotting the strain and SNV frequencies of _B. vulgatus_, consult [demonstration.MD](https://github.com/garudlab/Wasney-Briscoe/blob/main/analysis/demonstration.MD)
