# Post-processing

During post-processing, merged MIDAS outputs are wrangled into data formats that can be used for downstream analyses (e.g., strain inference). The post-processing scripts are written in python 2.7. Here, we execute them in a python 2.7-based conda environment `python27_env`.

## Set up

Project paths and global parameters are defined in `config.py` file in the [`scripts/postprocessing/postprocessing_scripts/`](https://github.com/garudlab/Wasney-Briscoe-2024/tree/main/scripts/postprocessing_scripts/) directory.

This file is currently modifed to reflect the directory strucutre of this repository, and will work *sans* modifications if the repository is cloned to your home directory (i.e., `~/`). To Process in another location or use a different set of global parameters, make those changes to the `config.py` directly.

## Calculate core genes, shared genes, and gene prevalences

In this step, three different outputs of interest are produced:
- `core_genes.txt.gz`: lists of core genes for all species detected in the dataset. Core genes are genes present in $\ge$ 90% of strains.
- `shared_genes.txt.gz`: lists of genes that are likely shared across species boundaries for all species detected in the dataset. Shared genes are genes that have a copy number $\gt$ 3, as this provides evidence for cross-species gene sharing. Putatively shared genes are removed from downstream analysis.
- `*species name*s__gene_freqs.txt.gz`: prevalences of each gene for all species detected in the dataset. A unique file is produced for each species.

For definitons of core and shared genes, see the methods of Wasney & Briscoe et al. and [Garud & Good et al., 2019](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000102).

Because this dataset essentially represents the amount of genetic diversity observed in a single individual's microbiome which mostly comprises one to two strains of given species, and designations of "core" and "shared" rely on the prevalence of genes in the broader population of strains we used definitions of core and shared genes calculated in Garud & Good, adding genes to the shared genes list (and removing them from the core genes list) if they showed evidence for cross-species sharing (i.e., copy number $c \gt 3$) in any sample within our dataset.

### Step 1: Create a `core_genes/` directory and add HMP files

In the `merged_data/` directory, create a directory `core_genes/`:

```
cd ~/merged_data/
mkdir core_genes/
```

In the [`scripts/postprocessing/`](https://github.com/garudlab/Wasney-Briscoe-2024/tree/main/scripts/postprocessing) directory, you will find the files `core_genes_HMP.txt.gz` and `shared_genes_HMP.txt.gz`. These are lists of core and shared genes defined in Garud & Good. **Move these to your. Copy these files into the `core_genes/` directory:

```
cp ~/Wasney-Briscoe-2024/scripts/postprocessing/core_genes_HMP.txt.gz ~/merged_data/core_genes/.
cp ~/Wasney-Briscoe-2024/scripts/postprocessing/shared_genes_HMP.txt.gz ~/merged_data/core_genes/.
```

Next, run `core_gene_utils.py`. The command below will alter the Garud & Good definitions of core and shared genes in the manner described above, and put out modified lists `core_genes.txt.gz` and `shared_genes.txt.gz`. From the [`scripts/postprocessing/`](https://github.com/garudlab/Wasney-Briscoe-2024/tree/main/scripts/postprocessing) directory, run:

```
conda activate python27_env
python core_gene_utils.py
```





