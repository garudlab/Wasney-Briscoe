# Post-processing

During post-processing, merged MIDAS outputs are wrangled into data formats that can be used for downstream analyses (e.g., strain inference). 

# Calculate core genes, shared genes, and gene prevalences

In this step, three different outputs of interest are produced:
- `core_genes.txt.gz`: lists of core genes for all species detected in the dataset. Core genes are genes present in $\ge$ 90% of strains.
- `shared_genes.txt.gz`: lists of genes that are likely shared across species boundaries for all species detected in the dataset. Shared genes are genes that have a copy number $\gt$ 3, as this provides evidence for cross-species gene sharing.
- `*species name*s__gene_freqs.txt.gz`: prevalences of each gene for all species detected in the dataset. A unique file is produced for each species.

For definitons of core and shared genes, see the methods of Wasney & Briscoe et al. and [Garud & Good et al., 2019](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000102).

Because this dataset essentially represents the amount of genetic diversity observed in a single individual's microbiome which mostly comprises one to two strains of given species, and designations of "core" and "shared" rely on the prevalence of genes in the broader population of strains we used definitions of core and shared genes calculated in Garud & Good, adding genes to the shared genes list (and removing them from the core genes list) if they showed evidence for cross-species sharing (i.e., copy number $c \gt 3$) in any sample within our dataset.

