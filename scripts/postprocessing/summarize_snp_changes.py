# Packages
import sys

sys.path.insert(0, "~/Wasney-Briscoe-2024/scripts/postprocessing/postprocessing_scripts/")

import config
import pandas as pd
import numpy as np

import subprocess

import pickle

import os

import itertools

## MIDAS postprocessing scripts
from calculate_intersample_changes import *
import parse_midas_data
import diversity_utils
import core_gene_utils

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

from annotation import *

# Make evolutionary changes directory if it doesn't already exist

evolutionary_changes_path = os.path.join(config.project_folder, "evolutionary_changes/")
if not os.path.exists(evolutionary_changes_path):
    os.makedirs(evolutionary_changes_path)

# Generating species list

species_list_path = "cat /u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/metadata/species_snps.txt.bz2"
species_list = subprocess.check_output(species_list_path, shell=True, stderr=subprocess.STDOUT).split("\n")
species_list = [species for species in species_list if species != ""]

# Calculating species to consider

good_species = []

for species in species_list:
    haploid_samples = diversity_utils.calculate_haploid_samples(species)
    if len(haploid_samples) > 1:
        good_species.append(species)
    else:
        continue
      
# Calculating which samples have the same strains for all species

different_strain_list = []
different_strain_trajectory = 1e-3
for species in good_species:
    intersample_change_map_extreme = load_intersample_change_map(species)
    haploid_samples = diversity_utils.calculate_haploid_samples(species)
    for sample_pair in intersample_change_map_extreme.keys():
        sample_1 = sample_pair[0]
        sample_2 = sample_pair[1]
        opportunities, pval, mutations, reversions = calculate_mutations_reversions_from_intersample_change_map(intersample_change_map_extreme, sample_1, sample_2, lower_threshold = 0.2, upper_threshold = 0.8)
        if (len(mutations + reversions)/opportunities > different_strain_trajectory) & (sample_1 in haploid_samples) & (sample_2 in haploid_samples):
            different_strain_list.append((species,sample_1,sample_2))
        else:
            continue

########################################################################################################################################################
################################################################ Extracting SNP changes ################################################################
########################################################################################################################################################
sample_metadata_map = parse_midas_data.parse_sample_metadata_map()
samples = sample_metadata_map.keys()

#Initialize df
snp_changes_df = pd.DataFrame(columns = ['species','contig','site_pos','gene','variant_type','sample1', 'sample2', 'alternate_freq_1', 'depth_1', 'alternate_freq_2', 'depth_2', 'gene_description'])
species_vec = []
contig = []
site_pos = []
gene = []
variant_type = []
sample1 = []
sample2 = []
alternate_freq_1 = []
depth_1 = []
alternate_freq_2 = []
depth_2 = []
opportunities_vec = []
gene_description_vec = []

counter = 0
no_of_species = len(good_species)

for species in species_list:

  counter += 1
    
  print "%s%s (%s/%s species)" % ("Processing ", species, counter, no_of_species)
  
  #gene descriptions
  genome_ids = parse_midas_data.get_ref_genome_ids(species)
  non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species)
  gene_descriptions = parse_patric.load_patric_gene_descriptions(genome_ids, non_shared_genes)
  centroid_gene_map = parse_midas_data.load_centroid_gene_map(species)


  #haploid samples
  haploid_samples = diversity_utils.calculate_haploid_samples(species)
  
  #change map
  intersample_change_map = load_intersample_change_map(species)
  
  #Intersample change map subsetting
#     ICMs = {key: intersample_change_map[key]['snps'] for key in intersample_change_map.keys() if len(intersample_change_map[key]['snps'][2]) > 0} #do we actually want this?
  ICMs = {key: intersample_change_map[key]['snps'] for key in intersample_change_map.keys()}
  ICMs = {key: ICMs[key] for key in ICMs.keys() if ((key[0] in haploid_samples) & (key[1] in haploid_samples))}
  
  
  if len(ICMs) != 0:
      for key in ICMs.keys():
          for snp in ICMs[key][2]:
              species_vec.append(species)
              contig.append(snp[1])
              site_pos.append(snp[2])
              gene.append(snp[0])
              variant_type.append(snp[3])
              sample1.append(key[0])
              sample2.append(key[1])
              alternate_freq_1.append(snp[4])
              depth_1.append(snp[5])
              alternate_freq_2.append(snp[6])
              depth_2.append(snp[7])
              if snp[0] in gene_descriptions:
                  gene_description_vec.append(gene_descriptions[snp[0]])
              elif snp[0] in centroid_gene_map:
                  if centroid_gene_map[snp[0]] in gene_descriptions:
                      gene_description_vec.append(gene_descriptions[centroid_gene_map[snp[0]]])
                  else:
                      gene_description_vec.append("")
              else:
                    gene_description_vec.append("")
  else:
      continue

snp_changes_df['species'] = species_vec
snp_changes_df['contig'] = contig
snp_changes_df['site_pos'] = site_pos
snp_changes_df['gene'] = gene
snp_changes_df['variant_type'] = variant_type
snp_changes_df['sample1'] = sample1
snp_changes_df['sample2'] = sample2
snp_changes_df['alternate_freq_1'] = alternate_freq_1
snp_changes_df['depth_1'] = depth_1
snp_changes_df['alternate_freq_2'] = alternate_freq_2
snp_changes_df['depth_2'] = depth_2  
snp_changes_df['gene_description'] = gene_description_vec   

snp_changes_df['mouse_1'] =  [extract_mouse_number(sample) for sample in snp_changes_df['sample1']]
snp_changes_df['mouse_2'] =  [extract_mouse_number(sample) for sample in snp_changes_df['sample2']]

snp_changes_df['location_1'] =  [extract_gut_site(sample) for sample in snp_changes_df['sample1']]
snp_changes_df['location_2'] =  [extract_gut_site(sample) for sample in snp_changes_df['sample2']]

snp_changes_df['region_1'] =  [extract_region(sample) for sample in snp_changes_df['sample1']]
snp_changes_df['region_2'] =  [extract_region(sample) for sample in snp_changes_df['sample2']]

snp_changes_df['diet_1'] =  [extract_diet(sample) for sample in snp_changes_df['sample1']]
snp_changes_df['diet_2'] =  [extract_diet(sample) for sample in snp_changes_df['sample2']]

snp_changes_df['cage_1'] =  [extract_cage(sample) for sample in snp_changes_df['sample1']]
snp_changes_df['cage_2'] =  [extract_cage(sample) for sample in snp_changes_df['sample2']]

strain_orientation = []
for i,snp_change in snp_changes_df.iterrows():
  if ((snp_change['species'],snp_change['sample1'], snp_change['sample2']) in different_strain_list) | ((snp_change['species'],snp_change['sample2'], snp_change['sample1']) in different_strain_list):
      strain_orientation.append("Different strain")
  else:
      strain_orientation.append("Same strain")
    
snp_changes_df['strain_orientation'] = strain_orientation   

snp_changes_path = "%s%s%s" % (config.project_folder, "evolutionary_changes/", "snp_changes.txt.bz2")
snp_changes_df.to_csv(snp_changes_path, sep = ",")

#################################################################################################################################################################
################################################################ Extracting Opportunities #######################################################################
#################################################################################################################################################################

opportunities_df = pd.DataFrame(columns = ['species', 'sample_1', 'sample_2', 'opportunities'])

species_vec = []
sample_1 = []
sample_2 = []
opportunity_vec = []
for species in good_species:
    print "Processing " + species
    intersample_change_map = load_intersample_change_map(species)
    
    haploid_samples_inoc = diversity_utils.calculate_haploid_samples(species)
    haploid_samples = [sample for sample in haploid_samples_inoc if sample != "TL1gDNAshort"]
    mice = set([extract_mouse_number(sample) for sample in haploid_samples])
    sample_pairs = list(itertools.combinations(haploid_samples_inoc, 2))


    for sample_pair in sample_pairs:
        if sample_pair in intersample_change_map:
            opportunities = intersample_change_map[sample_pair]['snps'][0]
        elif (sample_pair[1], sample_pair[0]) in intersample_change_map:
            sample_pair = (sample_pair[1], sample_pair[0])
            opportunities = intersample_change_map[sample_pair]['snps'][0]
        else: 
            print "Error: " + str(sample_pair) + "not found."

        species_vec.append(species)
        sample_1.append(sample_pair[0])
        sample_2.append(sample_pair[1])
        opportunity_vec.append(opportunities)

opportunities_df['species'] = species_vec
opportunities_df['sample_1'] = sample_1
opportunities_df['sample_2'] = sample_2
opportunities_df['opportunities'] = opportunity_vec

#Annotation
opportunities_df['mouse_1'] =  [extract_mouse_number(sample) for sample in opportunities_df['sample_1']]
opportunities_df['mouse_2'] =  [extract_mouse_number(sample) for sample in opportunities_df['sample_2']]

opportunities_df['location_1'] =  [extract_gut_site(sample) for sample in opportunities_df['sample_1']]
opportunities_df['location_2'] =  [extract_gut_site(sample) for sample in opportunities_df['sample_2']]

opportunities_df['region_1'] =  [extract_region(sample) for sample in opportunities_df['sample_1']]
opportunities_df['region_2'] =  [extract_region(sample) for sample in opportunities_df['sample_2']]

opportunities_df['diet_1'] =  [extract_diet(sample) for sample in opportunities_df['sample_1']]
opportunities_df['diet_2'] =  [extract_diet(sample) for sample in opportunities_df['sample_2']]

opportunities_df['cage_1'] =  [extract_cage(sample) for sample in opportunities_df['sample_1']]
opportunities_df['cage_2'] =  [extract_cage(sample) for sample in opportunities_df['sample_2']]


host_orientation = []
gut_orientation = []
inoculum_orientation = []
diet_orientation = []
for i,row in opportunities_df.iterrows():
    if row['mouse_1'] == row['mouse_2']:
        host_orientation.append("Within host")
    elif (row['mouse_1'] == "Inoculum") | (row['mouse_2'] == "Inoculum"):
        host_orientation.append("Between inoculum")
    else:
        host_orientation.append("Between host")
        
    if row['region_1'] == row['region_2']:
        gut_orientation.append("Within region")
    elif (row['mouse_1'] == "Inoculum") | (row['mouse_2'] == "Inoculum"):
        gut_orientation.append("Between inoculum")
    else:
        gut_orientation.append("Upper gut vs. lower gut")
    
    if row["mouse_1"] == "Inoculum":
        inoculum_orientation.append(row["region_2"])
    elif row["mouse_2"] == "Inoculum":
        inoculum_orientation.append(row["region_1"])
    else: 
        inoculum_orientation.append("Within or between mice")
        
    if row["diet_1"] == row["diet_2"]:
        diet_orientation.append("Same diet")
    else:
        diet_orientation.append("Different diet")



opportunities_df['host_orientation'] =  host_orientation
opportunities_df['gut_orientation'] = gut_orientation
opportunities_df['diet_orientation'] = diet_orientation

opportunities_df['inoculum_orientation'] = inoculum_orientation
opportunities_df['cage_orientation'] = opportunities_df.apply(lambda row: "Between inoculum" if row['host_orientation'] == "Between inoculum" else "Same cage" if row['cage_1'] == row['cage_2'] else "Different cage", axis = 1)

strain_orientation = []
for i,opportunity in opportunities_df.iterrows():
    if ((opportunity['species'],opportunity['sample_1'], opportunity['sample_2']) in different_strain_list) | ((opportunity['species'],opportunity['sample_2'], opportunity['sample_1']) in different_strain_list):
        strain_orientation.append("Different strain")
    else:
        strain_orientation.append("Same strain")
    
opportunities_df['strain_orientation'] = strain_orientation   

opportunities_path = "%s%s%s" % (config.project_folder, "evolutionary_changes/", "opportunities.txt.bz2")
opportunities_df.to_csv(opportunities_path, sep = ",")


