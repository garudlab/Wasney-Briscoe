import sys
sys.path.insert(0, "/u/home/m/michaelw/project-ngarud/microbiome_evolution/microbiome_evolution_SHALON/")

import config
import pandas as pd
import numpy as np
from datetime import datetime

import subprocess

import pickle

import os

import itertools

import bz2
import argparse

#MIDAS postprocessing scripts
from calculate_intersample_changes import *
import parse_midas_data
import diversity_utils
import core_gene_utils
import parse_patric

import matplotlib.pyplot as plt
import seaborn as sns

#ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument("--within_host_changes", help="Extract SNVs for within host changes", action="store_true")
parser.add_argument("--chunk_size", type=int, help="max number of nucleotides to load", default=1e10)

args = parser.parse_args()

within_host_changes_condition = args.within_host_changes
chunk_size = args.chunk_size

#Calculating species list
species_list_path = "cat /u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/metadata/species_snps.txt"
species_list = subprocess.check_output(species_list_path, shell=True, stderr=subprocess.STDOUT).split("\n")
species_list = [species for species in species_list if species != ""]

good_species = []

for species in species_list:
    haploid_samples = diversity_utils.calculate_haploid_samples(species)
    if len(haploid_samples) > 1:
        good_species.append(species)
    else:
        continue
        
        
#Loading SNV changes dataframe
evolution_folder = "%s%s" % (config.project_directory, "evolutionary_changes/")
if within_host_changes_condition:
    snv_changes_df_output = "%s%s" % (evolution_folder, "snp_changes_WithinHost.txt.bz2")
else:
    snv_changes_df_output = "%s%s" % (evolution_folder, "snp_changes.txt.bz2")
snp_changes_df = pd.read_csv(snv_changes_df_output, sep = ",")


#Filtering for within-host changes
condition_1 = (snp_changes_df.sample_type_1 == "Capsule") & (snp_changes_df.sample_type_2 == "Capsule")

within_host_changes = snp_changes_df[condition_1].sort_values(by = ['species', 'contig','gene'])


#1. Identify species that display within host changes
within_host_species = within_host_changes.species.unique()

#2. for loop to capture all changes
freq_df_final = pd.DataFrame()

for species in within_host_species:
    #3. Identify subjects
    subjects = list(within_host_changes[(within_host_changes['species'] == species)].subject.unique())
    sys.stderr.write("Beginning to process %s\n" % (species))
    for subject in subjects:
        #4. extract SNV changes associated with that host
        sys.stderr.write("Processing subject %s in %s.\n" % (subject, species))
        #For the subject in the current loop, pull out relevant samples and create a metadata table
        subject = str(subject)
        subject_integer = int(subject)
        ##Extract haploid subject samples
        subject_samples = parse_midas_data.parse_subject_sample_map()[subject].keys()
        haploid_samples = diversity_utils.calculate_haploid_samples(species, min_coverage=10)
        subject_samples = [sample for sample in subject_samples if sample in haploid_samples]
        ##Make a metadataframe
        metadata_map = parse_midas_data.parse_sample_metadata_map()
        subject_sample_metadata = pd.DataFrame()
        subject_sample_metadata['sample'] = subject_samples
        subject_sample_metadata['sample_type'] = subject_sample_metadata['sample'].apply(lambda sample: metadata_map[sample][9])
        subject_sample_metadata['sample_set'] = subject_sample_metadata['sample'].apply(lambda sample: metadata_map[sample][8])
        subject_sample_metadata['location'] = subject_sample_metadata['sample'].apply(lambda sample: metadata_map[sample][10])
        subject_sample_metadata['swallow_date'] = subject_sample_metadata['sample'].apply(lambda sample: '' if metadata_map[sample][3] == '' else datetime.strptime(metadata_map[sample][3], "%Y-%m-%dT%H:%M:%SZ").date())
        subject_sample_metadata['collection_date'] = subject_sample_metadata['sample'].apply(lambda sample: metadata_map[sample][6])
        subject_sample_metadata['subject'] = subject
        subject_sample_metadata['species'] = species

        #Extract contig and locus information 
        loci_of_interest = set(within_host_changes[(within_host_changes['subject'] == subject_integer) & (within_host_changes['species'] == species)].sort_values(by = ['site_pos'])[['contig', 'site_pos']].drop_duplicates().to_records(index=False).tolist())
        no_of_loci = len(loci_of_interest)

        #Pull out allele frequencies of variants of interest in subject samples


        ## Load the bz2 file
        data_directory = config.data_directory
        snp_file =  bz2.BZ2File("%ssnps/%s/annotated_snps.txt.bz2" % (data_directory, species),"r")

        ## Setting up the loop
        num_sites_processed = 0
        num_extracted_sites = 0
        line_number = -1
        final_line_number = -1
        initial_line_number = -100
        previous_gene_name = ""
        gene_name = ""
        chunk_size = chunk_size

        ## Header info
        line = snp_file.readline() 
        items = line.split()[1:]    
        samples_in_file = sample_utils.parse_merged_sample_names(items)

        ## Sample indices
        desired_sample_idxs = []
        for sample in subject_samples:
            desired_sample_idxs.append( numpy.nonzero(samples_in_file==sample)[0][0] )
        desired_sample_idxs = numpy.array(desired_sample_idxs)    
        desired_samples = samples_in_file[desired_sample_idxs]

        ## Gene descriptions
        genome_ids = parse_midas_data.get_ref_genome_ids(species)
        non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species)
        gene_descriptions = parse_patric.load_patric_gene_descriptions(genome_ids, non_shared_genes)
        centroid_gene_map = parse_midas_data.load_centroid_gene_map(species)

        ## Initializing
        chrom_vec = []
        location_vec = []
        gene_descriptions_vec = []
        loci_extracted = 0


        for line in snp_file: 

            line_number += 1

            previous_gene_name = gene_name

            if line_number%1000000==0:
                sys.stderr.write("%dk sites processed...\n" % (line_number/1000)) 

            items = line.split()
            ## Load information about site
            info_items = items[0].split("|")
            chromosome = info_items[0]
            location = long(info_items[1])
            gene_name = info_items[2]
            variant_type = info_items[3]
            location_tuple = (chromosome, int(location))
            ## If it's not one of the sites of interest, move to the next line

            if line_number >= chunk_size and gene_name!=previous_gene_name:
                ## We are done for now!
                final_line_number = line_number
                sys.stderr.write("Breaking at line " + str(final_line_number) + "\n")
                break

            if location_tuple not in loci_of_interest:
                continue
            else:
                sys.stderr.write("EXTRACTING: " + str(location_tuple) + "\n")
                chrom_vec.append(chromosome)
                location_vec.append(location)
                loci_extracted += 1

            ## Getting alts and depths for sites of interest
            alts = []
            depths = []

            for idx in desired_sample_idxs:    
                item = items[1+idx]
                subitems = item.split(",")
                alts.append(float(subitems[0]))
                depths.append(float(subitems[1]))

            if loci_extracted == 1:
                alts_total = numpy.array(alts)
                depths_total = numpy.array(depths)
            else:
                alts_total = numpy.vstack((alts_total, numpy.array(alts)))
                depths_total = numpy.vstack((depths_total, numpy.array(depths)))

            if gene_name in gene_descriptions:
                gene_descriptions_vec.append(gene_descriptions[gene_name])
            elif gene_name in centroid_gene_map:
                if centroid_gene_map[gene_name] in gene_descriptions:
                    gene_descriptions_vec.append(gene_descriptions[centroid_gene_map[gene_name]])
            else:
                gene_descriptions_vec.append("")


            num_sites_processed += 1

            if num_sites_processed == no_of_loci:
                sys.stderr.write("Successfully extracted all sites.\n")
                break




        snp_file.close()

        ## creating dataframe with alts, depths, and freqs for each sample
        if no_of_loci == 1:
            alts_df = pd.DataFrame(data = [alts_total], columns = desired_samples)
            depths_df = pd.DataFrame(data = [depths_total], columns = desired_samples)
        else:
            alts_df = pd.DataFrame(data = alts_total, columns = desired_samples)
            depths_df = pd.DataFrame(data = depths_total, columns = desired_samples)
        alts_df['contig'] = chrom_vec
        alts_df['site_pos'] = location_vec
        alts_df['gene_descriptions'] = gene_descriptions_vec
        depths_df['contig'] = chrom_vec
        depths_df['site_pos'] = location_vec
        alts_df = alts_df.melt(var_name='sample', value_name='alt', id_vars = ['contig', 'site_pos', 'gene_descriptions'])
        depths_df = depths_df.melt(var_name='sample', value_name='depth', id_vars = ['contig', 'site_pos'])
        freq_df = pd.merge(alts_df, depths_df, on=['sample', 'site_pos', 'contig'])
        freq_df['allele_frequency'] = freq_df['alt']/freq_df['depth']
        ## Merge with subject_sample map
        freq_df = pd.merge(freq_df, subject_sample_metadata, on = 'sample')

        freq_df_final = freq_df_final.append(freq_df)

        sys.stderr.write("Done with subject %s!\n\n" % (subject))
    
    sys.stderr.write("######### Done with %s!\n\n" % (species))
    
# SAVING
evolution_folder = "%s%s" % (config.project_directory, "evolutionary_changes/")
if not os.path.exists(evolution_folder):
    # Create the directory if it doesn't exist
    os.makedirs(evolution_folder)
    print("Directory '{}' created successfully.".format(evolution_folder))
else:
    print("Directory '{}' already exists.".format(evolution_folder))
snv_frequencies_output = "%s%s" % (evolution_folder, "snv_frequencies.txt.bz2")
freq_df_final.to_csv(snv_frequencies_output, index = False, sep = "\t")


