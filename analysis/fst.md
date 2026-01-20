# Calculating $F_{ST}$ for bacterial species

## step 1. Inferring $F_{ST}$

Code for calculating pi can be found at [`Wasney-Briscoe/scripts/pi/StatsPipeline/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/) are pertinent to the the pi calculation pipeline. they are:
1. [`qsubs/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/qsubs/)
2. [`configs`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/configs/)
3. [`pop_gen_calculator`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/pop_gen_calculator/)

Similar to the $\pi$ calculation pipeline, you will need the following files to run the $F_{ST}$ pipeline:
- A list of blacklist genes [`cross_species_centroids_clean.csv`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/BlacklistGenes/cross_species_centroids_clean.csv), which includes genes that share >95% average nucleotide identity across species.
- a list of species to process. This pipeline is set up to use [`species_snps.txt`](https://github.com/garudlab/Wasney-Briscoe/tree/main/metadata/species_snps.txt).
- a list of genome lengths for each species lists in [`species_snps.txt`](https://github.com/garudlab/Wasney-Briscoe/tree/main/metadata/species_snps.txt). This pipeline is set up to use [`species_lengths.txt`](https://github.com/garudlab/Wasney-Briscoe/tree/main/metadata/species_lengths.txt). This file lists each species genome length as 20 Mb, which far exceeds the genome size of bacteria found in the human gut. Because the pi-calculation algorithm processes diversity at each nucleotide site, setting such a generously high genome length ensures the algorithm will not be terminated before the end of a bacterial species' genome.
- a list of sample accessions. This pipeline is set up to use [`accessions.txt`](https://github.com/garudlab/Wasney-Briscoe/tree/main/metadata/accessions.txt)
- a [`config.yaml`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/configs/config.yaml) specifying the relevant paths.

Run the following job, which will submit a job array from the [`qsubs/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/qsubs/), which will submit a number of jobs corresponding to the number of species-sample combinations:

```
sh ./qsub_calculate_fst_EXECUTE.sh
```

This calculates $F_{ST}$ for the downsampled MIDAS data. By default, this calculates [qsub_calculate_fst_EXECUTE.sh](https://github.com/garudlab/Wasney-Briscoe/blob/main/scripts/pi/StatsPipeline/qsubs/qsub_calculate_fst_EXECUTE.sh) for normal MIDAS outputs. To calculate $F_{ST}$ on downsampled MIDAS data, change `study=normal` to `study=SiteDownsampled`. The paths to the normal MIDAS data are specified in the [config.yaml](https://github.com/garudlab/Wasney-Briscoe/blob/main/analysis/pi.md) file.'

## Step 2: step 4. Create a single $F_{ST}$ dataframe
