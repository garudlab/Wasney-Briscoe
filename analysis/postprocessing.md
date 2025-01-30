# Post-processing

During post-processing, merged MIDAS outputs are wrangled into data formats that can be used for downstream analyses (e.g., strain inference). The post-processing scripts are written in python 2.7. Here, we execute them in a python 2.7-based conda environment `python27_env`.

## Set up

Project paths and global parameters are defined in `config.py` file in the [`Wasney-Briscoe/scripts/postprocessing/postprocessing_scripts/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/postprocessing/postprocessing_scripts) directory.

This file is currently modifed to reflect the directory strucutre of this repository, and will work *sans* modifications if the repository is cloned to your home directory (i.e., `~/`). To Process in another location or use a different set of global parameters, make those changes to the `config.py` directly.

## Step 1: Calculate core genes, shared genes, and gene prevalences

In this step, three different outputs of interest are produced:
- `core_genes.txt.gz`: lists of core genes for all species detected in the dataset. Core genes are genes present in $\ge$ 90% of strains.
- `shared_genes.txt.gz`: lists of genes that are likely shared across species boundaries for all species detected in the dataset. Shared genes are genes that have a copy number $\gt$ 3, as this provides evidence for cross-species gene sharing. Putatively shared genes are removed from downstream analysis.
- `species_id_gene_freqs.txt.gz`: prevalences of each gene for all species detected in the dataset. A unique file is produced for each species.

For definitons of core and shared genes, see the methods of Wasney & Briscoe et al. and [Garud & Good et al., 2019](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000102).

Because this dataset essentially represents the amount of genetic diversity observed in a single individual's microbiome which mostly comprises one to two strains of given species, and designations of "core" and "shared" rely on the prevalence of genes in the broader population of strains we used definitions of core and shared genes calculated in Garud & Good, adding genes to the shared genes list (and removing them from the core genes list) if they showed evidence for cross-species sharing (i.e., copy number $c \gt 3$) in any sample within our dataset.

### Create a `core_genes/` directory and add HMP files

In the `merged_data/` directory, create a directory `core_genes/`:

```
cd ~/merged_data/
mkdir core_genes/
```

In the [`Wasney-Briscoe/scripts/postprocessing/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/postprocessing) directory, you will find the files `core_genes_HMP.txt.gz` and `shared_genes_HMP.txt.gz`. These are lists of core and shared genes defined in Garud & Good. **Move these to your. Copy these files into the `core_genes/` directory:

```
cp ~/Wasney-Briscoe-2024/scripts/postprocessing/core_genes_HMP.txt.gz ~/merged_data/core_genes/.
cp ~/Wasney-Briscoe-2024/scripts/postprocessing/shared_genes_HMP.txt.gz ~/merged_data/core_genes/.
```

### Run `core_gene_utils.py`

Next, run `core_gene_utils.py`. The command below will alter the Garud & Good definitions of core and shared genes in the manner described above, and put out modified lists `core_genes.txt.gz` and `shared_genes.txt.gz`. From the [`Wasney-Briscoe/scripts/postprocessing/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/postprocessing) directory, run:

```
conda activate python27_env
python core_gene_utils.py
```

## Step 2: Process nucleotide diversity

In this step, files that summarize nucleotide diversity within and across samples are produced.

To run this pipeline, navigate to the [`Wasney-Briscoe/scripts/postprocessing/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/postprocessing) directory and run:

```
qsub post_processing_wrapper.sh
```

This pipeline is performed on all species for which MIDAS-processed SNP data is available (in `~/merged_data/snps/`). A list of those species is present as `species_snps.txt` in [`Wasney-Briscoe/scripts/postprocessing/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/postprocessing) directory. To process a subset of these species, assemble a new list file and pass the path to the `post_processing_wrapper.sh` script using the `-s` or `--species_list` flags, e.g., 

```
qsub post_processing_wrapper.sh -s ~/new_species_list.txt
```

Note that while the normal Garud & Good pipeline uses the panel in it's own dataset to calculate average rates of between strain diversity per sample and species, we use the Human Microbiome Project (HMP) panel used in Garud & Good to calculate average rates of between host diversity (using files found in the `Wasney-Briscoe-2024/scripts/postprocessing/HMP_snp_prevalences/` directory(. To use the mouse dataset itself to estimate rates of between strain diversity, remove the `--use_HMP` flag when running the `calculate_within_person_sfs.py` script in the pipeline. However, doing so will render it impossible to detect monocolonized samples for species that have only a single strain that exists across the entire dataset.

`post_processing_wrapper.sh` should produce the following files:
- In `~/merged_data/snps/`:
  - `marker_coverage.txt.bz2`
  - `gene_coverage.txt.bz2`
  - `coverage_distribution.txt.bz2`
  - `full_coverage_distribution.txt.bz2`
  - `annotated_snps.txt.bz2`
  - `within_sample_sfs.txt.bz2`
- In `~/merged_data/snp_prevalences/` (created by the pipeline)
  - `*species_id*.txt.gz`, where *species_id* is the species id for all species passed to the script in `species_snps.txt`. Note that these files will not be used to generate `within_sample_sfs.txt.bz2` unless removing the  `--use_HMP` flag when running the `calculate_within_person_sfs.py` script.

## Step 3: Calculate evolutionary changes between samples

### Identify SNP changes

Next, we identify two types of evolutionary changes between all pairs of samples:
- SNPs going from low frequency (allele frequency $f \le 0.2$) in one sample to high frequency in another ($f \ge 0.8$)
- Genes going from 0 copies (copy number $c \le 0.05$) in one sample to 1 copy ($0.6 \le c \le 1.2$) in another.

To calculate evolutionary changes, navigate to the [`Wasney-Briscoe/scripts/postprocessing/postprocessing_scripts`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/postprocessing/postprocessing_scripts) directory and run:

```
python calculate_intersample_changes.py
```


Alternatively, from the [`Wasney-Briscoe/scripts/postprocessing/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/postprocessing/) directory, you can run submit the `calculate_intersample_changes.py` script as a job:

```
qsub ./calculate_intersample_changes_WRAPPER.sh
```

`calculate_intersample_changes.py` should produce the following files:
- In `~/merged_data/intersample_change/` (created by the pipeline)
  - `species_id.txt.gz`, where *species_id* is the species id for all species passed to the script in `species_snps.txt`.
  
### Summarize SNP changes and opportunities in dataframe format

Downstream steps require SNP changes to be summarized in a dataframe. To do this, run the [`summarize_snp_changes.py`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/postprocessing/summarize_snp_changes.py) script from the [`Wasney-Briscoe/scripts/postprocessing/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/postprocessing/) directory:

```
conda activate python27_env #If python 2.7 isn't already loaded 
python summarize_snp_changes.py
```

`summarize_snp_changes.py` will make the following: 
- A directory called `evolutionary_changes` in your project folder (set to `~/` in the [`config.py`](https://github.com/garudlab/Wasney-Briscoe/blob/main/scripts/postprocessing/postprocessing_scripts/config.py). Within that directory:
  - `snp_changes.txt.bz2`: a dataframe containing all SNP changes between samples (i.e., SNPs going from allele frequency $f \le 0.2$ in one sample to $f \ge 0.8$ in another sample).
  - `opportunities.txt.bz2`: a dataframe quantifying the the number of loci that have high coverage (i.e., coverage $D \ge 20$ reads) in both samples between all pairs of QP samples.
