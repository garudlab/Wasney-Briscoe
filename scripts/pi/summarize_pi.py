import sys
sys.path.insert(0, "~/Wasney-Briscoe/scripts/postprocessing/postprocessing_scripts/")
import config

sys.path.insert(0, "~/Wasney-Briscoe/scripts/supplementary_figures/helper_scripts/")
from annotation import *

mport os

import numpy as np
import pandas as pd

import scipy as sc

cov = 4

single_sample_pi_dir = "%sSinglePi_MinCoverage%s/" % (config.project_folder, str(cov))

single_sample_pi_output = "~/SinglePi_SchloissnigPi_cov%s.csv" % (config.data_directory, str(cov))

single_pi_df = pd.DataFrame(columns = ["species", 
                                       "sample", 
                                       "genomewide_pi", 
                                       "variable_sites", 
                                       "mean_depth", 
                                       "total_loci"])

species_vec = []
sample_vec = []
Genomewide_pi = []
Genomewide_pi_variable_sites = []
Mean_depth = []
n_total_loci = []

for filename in os.listdir(single_sample_pi_dir):
    if "Loci_Stats" in filename:
        continue
    #IDing species
    species = "_".join(filename.split("_")[1:4])
    #IDing sample name
    start_index = filename.find("SampleID1_") + len("SampleID1_")
    end_index = filename.find("_Pi")
    sample = filename[start_index:end_index]
    
    #Loading data
    summary_stats = pd.read_csv("%s%s" % (single_sample_pi_dir, filename), index_col = 0)
    
    #creating vectors
    species_vec.append(species)
    sample_vec.append(sample)
    Genomewide_pi.append(summary_stats['Genomewide_pi'].values[0])
    Genomewide_pi_variable_sites.append(summary_stats['Genomewide_pi_variable_sites'].values[0])
    Mean_depth.append(summary_stats['Mean_depth'].values[0])
    n_total_loci.append(summary_stats['n_total_loci'].values[0])

single_pi_df["species"] = species_vec
single_pi_df["sample"] = sample_vec
single_pi_df["genomewide_pi"] = Genomewide_pi
single_pi_df["variable_sites"] = Genomewide_pi_variable_sites
single_pi_df["mean_depth"] = Mean_depth
single_pi_df["total_loci"] = n_total_loci

#Annotation functions

single_pi_df['mouse'] = single_pi_df['sample'].apply(lambda sample: extract_mouse_number(sample))
single_pi_df['cage'] = single_pi_df['sample'].apply(lambda sample: extract_cage(sample))
single_pi_df['diet'] = single_pi_df['sample'].apply(lambda sample: extract_diet(sample))
single_pi_df['gut_site'] = single_pi_df['sample'].apply(lambda sample: extract_gut_site(sample))
single_pi_df['gut_region'] = single_pi_df['sample'].apply(lambda sample: extract_region(sample))
single_pi_df["good_species"] = [True if species in good_species else False for species in single_pi_df.species]

single_pi_df.to_csv(single_sample_pi_output, index=False)