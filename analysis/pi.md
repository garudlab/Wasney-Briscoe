# Calculating $\pi$ for bacterial species

## Step 1. Downsampling nucleotide sites with high coverage in mice 7 and 8

This code produces a midas-format snps_ref_freq.txt.bz2 and snps_depth.txt.bz2 where the columns for samples in Mouse 7 and Mouse 8 are updated on a genomic site-by-genomic site basis if the depth of those samples is above the median depth for samples from Mouse 1 through 6. If the median depth for Mouse 1 through 6 is less than 4, then the sample depth in Mouse 7 and 8 are automatically converted to 0 at that site because a site with less than 4 depth will not be considered in $\pi$ calcultion regardless. Note that $\pi$ can be inferred from MIDAS outputs that have not been downsampled (i.e., step 3 onwards), but we have found that downsampling is an import step for removing bias introduced by read coverage differences between mice 7 and 8 relative to the other mice, as mice 7 and 8 were sequenced separately and to higher depths.

Input: `snps_ref_freq.txt.bz2` and `snps_depth.txt.bz2`
Output: snps_ref_freq.txt.bz2 and snps_depth.txt.bz2 with modifications to columns corresponding to Mouse 7 and 8. These files can be found in a new directory: 

At each nucleotide site, the script [optimized_site_downsampling.py](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/optimized_site_downsampling.pi) does the following
1. Compute median depth of reads for all samples from Mice 1 through 6
2. If median depth < 4, pass to next site. Else, continue
3. This creates a vector of 0s and 1s, representing the total pool of reference and alternative nucleotides for a given sample, respectively.
2. Based on median

To run fo all species with SNP data, submit the [optimized_site_downsampling.py](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/optimized_site_downsampling.pi) script for multiple species simultaneously by executing the following code from the [Wasney-Briscoe/scripts/pi/](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/) directory:

```
while IFS= read -r species;
do
  	qsub -cwd -V -N $species -e logs -o logs -l h_data=4G,time=12:00:00 -b y "./job_script_site_downsampling.sh $species" 
done < ~/Wasney-Briscoe/metadata/species_snps.txt
```

## Step 2. Copying files from normal MIDAS output directory to the downsampled directory

Run [`cp_files_for_pi.sh`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/cp_files_for_pi.sh) in [~/Wasney-Briscoe/scripts/pi/](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/):

```
sh ./cp_files_for_pi.sh
```

This script creates a folder called `merged_data_downsampled` in your home directory (`~/`) and copy the MIDAS files `snps_info.txt.bz2` and `genes_copynum.txt.bz2` for all species from the normal MIDAS `merged_data/` output directory to this new directory. These files are subsequently used to infer pi. 

## step 3. Inferring $\pi$

Code for calculating pi can be found at [`Wasney-Briscoe/scripts/pi/StatsPipeline/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/) are pertinent to the the pi calculation pipeline. they are:
1. [`qsubs/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/qsubs/)
2. [`configs`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/configs/)
3. [`pop_gen_calculator`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/pop_gen_calculator/)

Here's what you'll need to run this pipeline, which are included in this repository:
- A list of blacklist genes [`cross_species_centroids_clean.csv`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/BlacklistGenes/cross_species_centroids_clean.csv), which includes genes that share >95% average nucleotide identity across species.
- a list of species to process. This pipeline is set up to use [`species_snps.txt`](https://github.com/garudlab/Wasney-Briscoe/tree/main/metadata/species_snps.txt).
- a list of genome lengths for each species lists in [`species_snps.txt`](https://github.com/garudlab/Wasney-Briscoe/tree/main/metadata/species_snps.txt). This pipeline is set up to use [`species_lengths.txt`](https://github.com/garudlab/Wasney-Briscoe/tree/main/metadata/species_lengths.txt). This file lists each species genome length as 20 Mb, which far exceeds the genome size of bacteria found in the human gut. Because the pi-calculation algorithm processes diversity at each nucleotide site, setting such a generously high genome length ensures the algorithm will not be terminated before the end of a bacterial species' genome.
- a list of sample accessions. This pipeline is set up to use [`accessions.txt`](https://github.com/garudlab/Wasney-Briscoe/tree/main/metadata/accessions.txt)
- a [`config.yaml`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/configs/config.yaml) specifying the relevant paths.

Run the following job, which will submit a job array from the [`qsubs/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/qsubs/), which will submit a number of jobs corresponding to the number of species-sample combinations:

```
sh ./qsub_single_sample_pi_EXECUTE.sh
```

This calculates pi for the downsampled MIDAS data. To calculate pi for the normal MIDAS outputs, in the [`qsub_single_sample_pi_EXECUTE.sh`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/qsubs/qsub_single_sample_pi_EXECUTE.sh) file, change `study=SiteDownsampled` to `study=normal`. The paths to the normal MIDAS data are specified in the [`config.yaml`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/configs/config.yaml) file.'

In addition, this script assumes the MIDAS database of reference genomes (`midas_db_v1.2`) is located in your home directory. If elsewhere, modify the path in [`config.yaml`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/configs/config.yaml).

## step 4. Create a single $\pi$ dataframe

Step 3 outputs an individual file for each species-sample combination. To create a single dataframe containing pi estimates for all species and samples, run:

```
conda activate python_env
python ~/Wasney-Briscoe/scripts/pi/summarize_pi.py
```

This will output a .csv file called `SinglePi_SchloissnigPi_cov4.csv` in your home directory (`~/`).






