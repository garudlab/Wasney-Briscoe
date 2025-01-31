# Running MIDAS

In the standard MIDAS workflow, `species`, `gene`, and `SNP` content are called by aligned metagenomic sequencing to a database comprising reference genomes representing 5,952 species. In the final step, outputs are merged across all relevant hosts.

Before executing these steps, please download MIDAS (v.1). Instructions to do so, as well as more information about each step can be found in the [MIDAS documentation](https://github.com/snayfach/MIDAS) (Nayfach et al., 2016). All scripts referenced here can be found in [`Wasney-Briscoe/scripts/MIDAS/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/MIDAS/) directory, and were executed on Hoffman2 using a job scheduler. Make sure to modify the scripts to reflect the location of the MIDAS software and the MIDAS database on your machine. 

MIDAS runs on `python 2.7`. Therefore, all MIDAS code is executed in a custom conda environment called `python27_env` that has `python 2.7` loaded. Further, note that these scripts assume that MIDAS software (`MIDAS_mod`) and the MIDAS database of common gut commensal reference genomes (`midas_db_v1.2`) are located in your home directory (`~/`). Make sure the paths in the [`midas_species.sh`](https://github.com/garudlab/Wasney-Briscoe/blob/main/scripts/MIDAS/MIDAS_species.sh), [`midas_genes.sh`](https://github.com/garudlab/Wasney-Briscoe/blob/main/scripts/MIDAS/MIDAS_species.sh), and [`midas_snps.sh `](https://github.com/garudlab/Wasney-Briscoe/blob/main/scripts/MIDAS/MIDAS_species.sh) are updated to reflect where your MIDAS software and database are actually located.

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
qsub midas_snps.sh
```

## 4. Merge MIDAS outputs

Finally, we merge MIDAS outputs across relevant hosts:

```
qsub midas_merge.sh
```

To save space, we compress files into `.bz2` format:

```
qsub bzip2_MIDAS_outputs.sh
```

Merged outputs were used in all subsequent analyses performed in this paper. 
