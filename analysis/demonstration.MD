# Demonstration with example data

We have provided examples of MIDAS species, genes, and snps outputs. The species outputs contain information on all species detected in this mouse dataset; contrastingly, we have only provided the MIDAS genes and snps output files for _B. vulgatus_ (Bacteroides_vulgatus_57955) due to the size of these files. As a result, the provided species files can be used to produce all plots related to taxonomic composition and diversity (e.g., figure 2; supplementary figure 1). Meanwhile, the genes and snps outputs can be used to to run the post-processing, $\pi$, and strain phasing pipelines as well as produce strain and SNV plots for _B. vulgatus_ only. Below, we provide step-by-step unstructions on how to demonstrate these pipelines using the provided _B. vulgatus_-data.

## Step 1: Move example data to the correct directory

This repository has been set up with the assumption that MIDAS data produced by the MIDAS pipeline will be output to your home directory (`~`). To use the example data we've provided, manually move it to your home directory like so:

```
mv ~/Wasney-Briscoe/example_data/merged_data/ ~/.
```

Subsequent steps assume that the `merged_data/` directory can be found in your home directory.

## Step 2: Post-processing pipeline

The post-processing pipeline can be executed almost identically as described in [postprocessing.MD](https://github.com/garudlab/Wasney-Briscoe/blob/main/analysis/postprocessing.md), with the following changes:

### calculating core genes

Run `core_gene_utils.py` with the `--species` flag, i.e., 

```
python core_gene_utils.py --species Bacteroides_vulgatus_57955
```

### Running the post-processing script

Run `post_processing_wrapper.sh` with a species list comprising only `Bacteroides_vulgatus_57955`. Such a list has been included in [example_data](https://github.com/garudlab/Wasney-Briscoe/tree/main/example_data), such that the code can be run like so:

```
qsub post_processing_wrapper.sh -species_list ~/Wasney-Briscoe/example_data/Bacteroides_vulgatus_57955_list.txt
```


