# Strain Frequency Inference

Strain frequencies are inferred from SNV frequency trajectories across samples using a method developed in Roodgar et al., 2018 and extended in Wolff et al., 2023. 

## Step 1. Reformat SNV data for strain frequency inference pipeline

The pipeline to infer strain frequencies was developed based on a previous method put forth in Roodgar et al., 2021 and extended in Wolff et al., 2023.

The strain frequency inference pipeline requires data be formated in the same manner as inputs to the [StrainFinder](https://github.com/cssmillie/StrainFinder) software package (Smillie et al., 2018). To reformat data like so, navigate to the [`Wasney-Briscoe/scripts/postprocessing/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/postprocessing/) directory and run:

```
qsub ./create_StrainFinderInput_wrapper.sh
```
`create_StrainFinderInput.py` and its wrapper script will make the following: 
- A directory/subdirectory with the path `/strain_phasing/input` in your project folder (set to `~/` in the [`config.py`](https://github.com/garudlab/Wasney-Briscoe/blob/main/scripts/postprocessing/postprocessing_scripts/config.py). Within that directory:
  - directories for all species processed (defined by the [`species_snps.txt`](https://github.com/garudlab/Wasney-Briscoe/blob/main/metadata/species_snps.txt) file). Within each species folder:
    - `species_id.strainfinder.locations.p`: a pickle file...
    - `species_id.strainfinder.p`: a pickle file...
    - `species_id.strainfinder.samples.p`: a pickle file...

## Step 2. Infer strain frequency from SNV frequencies

To infer strain frequencies, execute the following code in the [`Wasney-Briscoe/scripts/strain_inference/`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/strain_inference) directory for the species of interest:

```
conda activate python_env
python strain_inference.py --species Bacteroides_vulgatus_57955
```

Ensure that the environment with python 3 (python_env) is loaded. If `--species` is not passed to the script, it will by default process Bacteroides_vulgatus_57955.

the `strain_inference.py` script should produce a file called `species_id_strain_frequency.csv` in the direcotry `strain_phasing/strain_clusters/species_id/`, where `species_id` is replaced by the actual ID of the species of interst. This file contains the frequncies of all inferred strains in each sample, as well as the upper and lower bounds of the bootstrapped 95% confidence interval around the strain frequncies in each sample.


