# Packages

import sys
sys.path.insert(0, "~/Wasney-Briscoe/scripts/postprocessing/postprocessing/scripts/")

import config
import pandas as pd
import numpy as np

import os

import subprocess

#MIDAS postprocessing scripts
from calculate_intersample_changes import *
import parse_midas_data
import diversity_utils
import core_gene_utils

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

sys.path.insert(0, "~/Wasney-Briscoe-2024/scripts/supplementary_figures/helper_scripts/")

from annotation import *

# Load species

species_list_path = "cat ~/Wasney-Briscoe-2024/metadata/species_snps.txt"
species_list = subprocess.check_output(species_list_path, shell=True, stderr=subprocess.STDOUT).split("\n")
species_list = [species for species in species_list if species != ""]

# Calculate haploid samples

ploid_df = pd.DataFrame(columns = ['species', 'sample', 'ploid'])

species_vec = []
sample_vec = []
ploid_vec = []

for species in species_list:
    samples = diversity_utils.calculate_highcoverage_samples(species)
    haploid_samples = diversity_utils.calculate_haploid_samples(species)
    for sample in samples:
        species_vec.append(species)
        sample_vec.append(sample)
        if sample in haploid_samples:
            ploid_vec.append("Haploid")
        else:
            ploid_vec.append("Polyploid")
    
ploid_df['species'] = species_vec
ploid_df['sample'] = sample_vec
ploid_df['ploid'] = ploid_vec

ploid_df['mouse'] = ploid_df['sample'].apply(lambda sample: extract_mouse_number(sample))
ploid_df['region'] = ploid_df['sample'].apply(lambda sample: extract_region(sample))

# Saving

out_path = "%s%s" % (config.analysis_directory, "metadata/sample_ploid_status.txt")
ploid_df.to_csv(out_path, index = False)
