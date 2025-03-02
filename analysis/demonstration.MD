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
qsub post_processing_wrapper.sh --species_list ~/Wasney-Briscoe/example_data/Bacteroides_vulgatus_57955_list.txt
```

Ensure that the number of tasks in the `post_processing_wrapper.sh` task array script reflects there being only a single species in the list (i.e., line 6 reads `#$ -t 1-1:1`)

### Calculating intersample changes

Run `calculate_intersample_changes.py` for Bacteroides_vulgatus_57955_list only like so:

```
python calculate_intersample_changes.py --species Bacteroides_vulgatus_57955
```

### Summarize SNP changes and opportunities in dataframe format

You pass the `summarize_snp_changes.py` script the custom species list with only Bacteroides_vulgatus_57955 like so:

```
python summarize_snp_changes.py --species_list  ~/Wasney-Briscoe/example_data/Bacteroides_vulgatus_57955_list.txt
```

### Extract frequencies of evolutionary SNPs from all samples

No change is necessary here, as the `extract_snv_freqs.py` will extract frequencies for all species included in the output of the last step.


## Step 3: calculating $\pi$

Implement the following modifications to the $\pi$ pipeline to calculate nucleotide diversity for Bacteroides_vulgatus_57955 only. Note that the scripts provided to produce Figure 3 work best with many species. Nonetheless, the dataframe output at the end of this pipeline `SinglePi_SchloissnigPi_cov4.csv` is of the format that is used to generate figure 3, albeit with information for only a single species instead of many.

### Downsampling

To downsample the _B. vulgatus_ snp data by site, pass the job submission code the list to the Bacteroides_vulgatus_57955-only species list:

```
while IFS= read -r species;
do
  	qsub -cwd -V -N $species -e logs -o logs -l h_data=4G,time=12:00:00 -b y "./job_script_site_downsampling.sh $species" 
done < ~/Wasney-Briscoe/example_data/Bacteroides_vulgatus_57955_list.txt
```

### Copying files from normal MIDAS output directory to the downsampled directory

Similarly, you can copy downsampled files into a new folder by passing the Bacteroides_vulgatus_57955-only species list:

```
sh ./cp_files_for_pi.sh --species_list ~/Wasney-Briscoe/example_data/Bacteroides_vulgatus_57955_list.txt
```

Ensure that the number of tasks specified reflects the presence of only a single species (i.e., line 8 should read `#$ -t 1-1:1`)

### Calculating pi

Instead of using the `qsub_single_sample_pi_EXECUTE.sh` in the [`qsubs/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/StatsPipeline/qsubs/) directory, submit the version in the [`example_data`](https://github.com/garudlab/Wasney-Briscoe/tree/main/example_data) directory. From the [`example_data`](https://github.com/garudlab/Wasney-Briscoe/tree/main/example_data) directory, you will submit it the same way:

```
sh ./qsub_single_sample_pi_EXECUTE.sh
```

### step 4. Create the $\pi$ dataframe

Run python script to summarize the dataframe as normal: 

```
conda activate python_env
python ~/Wasney-Briscoe/scripts/pi/summarize_pi.py
```















