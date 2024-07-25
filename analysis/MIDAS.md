# Running MIDAS

In the standard MIDAS workflow, `species`, `gene`, and `SNP` content are called by aligned metagenomic sequencing to a database comprising reference genomes representing 5,952 species. In the final step, outputs are merged across all relevant hosts.

Before executing these steps, please download MIDAS (v.1). Instructions to do so, as well as more information about each step can be found at [MIDAS documentation](https://github.com/snayfach/MIDAS) (Nayfach et al., 2016). All scripts referenced here can be found in [`scripts/MIDAS/`](https://github.com/garudlab/Wasney-Briscoe-2024/tree/main/scripts/MIDAS/) directory. All steps were performed on Hoffman2 using a job scheduler.

MIDAS runs on `python 2.7`. Therefore, all MIDAS code is executed in a custom conda environment called `python27_env` that has `python 2.7` loaded. 

## 1. Species step

First, we infer species presence/absence and relative abundance by running the following script:

```
qsub midas_species.sh
```

Before completing the `gene` and `SNP`, we produce "species union" files, which inform MIDAS which samples belong to the same host:

```
qsub construct_species_union_jobarray.sh
```

## 2. Gene step

Second, we infer gene copy number by running the following script:

```
qsub midas_genes.sh
```

## 3. SNP step

Third, we infer allele frequency for all loci in the genome:

```
qsub midas_genes.sh
```

## 4. Merge MIDAS outputs

Finally, we merge MIDAS outputs across relevant hosts:

```
qsub midas_merge.sh
```

To save space, we compress files into `.bz2` format:

```
qsub bzip2_MIDAS_outputs
```

Merged outputs were used in all subsequent analyses performed in this paper. 
