# Creating supplementary figures

## Supplementary Figure 1

Run `figure_S1.R` in [`Wasney-Briscoe/scripts/figures/figure_S1.R`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/figures/figure_S1.R).

This script will output `figure_S1.png` in your home directory (`~/`).

## Supplementary Figure 2B

First, we create a file called `sample_ploid_status.txt`, which contains information on whether a species has one (haploid) or multiple (polyploid) strains in each sample in which it is detected. To generate this file, run `calculate_ploid.py` at [Wasney-Briscoe/scripts/supplementary_figures/calculate_ploid.py](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/supplementary_figures/calculate_ploid.py), using python 2.7 to do so (e.g., in a python 2.7-based conda environment). `calculate_ploid.py` generates the file `sample_ploid_status.txt` in the `Wasney-Briscoe-2024/metadata/` subdirectory.

Next, run `figure_S2B.R` in [`Wasney-Briscoe/scripts/supplementary_figures/figure_S2B.R`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/supplementary_figures/figure_S2B.R).

This script will output `figure_S2B.png` in your home directory (`~/`).

## Supplementary Figure 3

Run `figure_S3.R` in [`Wasney-Briscoe/scripts/supplementary_figures/figure_S3.R`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/supplementary_figures/figure_S3.R).

This script will output `figure_S3.png` in your home directory (`~/`).

## Supplementary Figure 4

Run `figure_S4.R` in [`Wasney-Briscoe/scripts/figures/figure_S4.R`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/supplementary_figures/figure_S4.R).

This script will output `figure_S4.png` in your home directory (`~/`).

## Supplementary Figure 5

To generate the plot a single trajectory representing both strains of _B. uniformis_, execute the [`strain_inference.py`](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/strain_inference/strain_inference.py)) script with the following flags:

```
python strain_inference.py --species Bacteroides_uniformis_57318 --single_trajectory
```

The `--species` flag specifies the species to be plotted (Bacteroides_uniformis_57318), and the `--single_trajectory` creates a figure called `Bacteroides_uniformis_57318_single_strain_trajectory.png` in your home directory (`~/`) in which a single trajectory represents the frequency of one of two strains colonizing the mice. The figure produced by this code is equivalent to supplementary figure 5.

## Figure 6


