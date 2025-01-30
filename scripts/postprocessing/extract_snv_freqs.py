import pandas as pd
import numpy as np
import bz2

import sys
sys.path.insert(0,"~/Wasney-Briscoe/scripts/postprocessing/postprocessing_scripts/")

import config

#MIDAS postprocessing scripts
from calculate_intersample_changes import *
import parse_midas_data
import diversity_utils
import core_gene_utils
import parse_patric

sys.path.insert(0,"~/Wasney-Briscoe/scripts/postprocessing/postprocessing_scripts/")

from annotation import *

import matplotlib.pyplot as plt

# All changes
all_changes_path = "%sevolutionary_changes/snp_changes.txt.bz2" % (config.project_folder)
all_changes = pd.read_csv(all_changes_path, sep = ",", index_col = 0)[['species', 'contig','site_pos']].drop_duplicates()

only_haploid = False

# species = "Ruminococcus_sp_58571"
min_coverage = 20

full_freq_df = pd.DataFrame()


for species in all_changes.species.unique():
# for species in ['Parabacteroides_distasonis_56985']:
    sys.stderr.write("Processing " + species + "\n")
    
    #Loci of interest
    loci_of_interest = set(all_changes[all_changes.species == species].sort_values(by = ['site_pos'])[['contig', 'site_pos']].drop_duplicates().to_records(index=False).tolist()
    )

    ##Extract haploid samples
    haploid_samples = diversity_utils.calculate_haploid_samples(species, min_coverage=min_coverage)
    high_coverage_samples = diversity_utils.calculate_highcoverage_samples(species, min_coverage=min_coverage)

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
    chunk_size = 1e10

    ## Header info
    line = snp_file.readline() 
    items = line.split()[1:]    
    samples_in_file = sample_utils.parse_merged_sample_names(items)

    # Sample indices
    desired_sample_idxs = []
    if only_haploid:
        for sample in haploid_samples:
            desired_sample_idxs.append( numpy.nonzero(samples_in_file==sample)[0][0] )
    else:
        for sample in high_coverage_samples:
            desired_sample_idxs.append( numpy.nonzero(samples_in_file==sample)[0][0] )
    
    desired_sample_idxs = numpy.array(desired_sample_idxs)    
    desired_samples = samples_in_file[desired_sample_idxs]

    ### gene annotation tools
    genome_ids = parse_midas_data.get_ref_genome_ids(species)
    non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species)
    gene_descriptions = parse_patric.load_patric_gene_descriptions(genome_ids, non_shared_genes)
    centroid_gene_map = parse_midas_data.load_centroid_gene_map(species)

    ## Initializing
    chrom_vec = []
    location_vec = []
    gene_vec = []
    gene_description_vec = []
    loci_extracted = 0

    for line in snp_file: 

        line_number += 1

        previous_gene_name = gene_name

    #     if line_number%100000==0:
    #         sys.stderr.write("%dk sites processed...\n" % (line_number/1000)) 

        items = line.split()
        # Load information about site
        info_items = items[0].split("|")
        chromosome = info_items[0]
        location = long(info_items[1])
        gene_name = info_items[2]
        variant_type = info_items[3]
        location_tuple = (chromosome, int(location))
        #If it's not one of the sites of interest, move to the next line

    #     if (chromosome == "CAHL01000016") & (location == 6905):
    #         print "encountered locus\n"

        if line_number >= chunk_size and gene_name!=previous_gene_name:
            # We are done for now!
            final_line_number = line_number
            sys.stderr.write("Breaking at line " + str(final_line_number) + "\n")
            break

        if location_tuple not in loci_of_interest:
            continue
        else:
            sys.stderr.write("EXTRACTING: " + str(location_tuple) + "\n")
            chrom_vec.append(chromosome)
            location_vec.append(location)
            gene_vec.append(gene_name)
            if gene_name in gene_descriptions:
                gene_description_vec.append(gene_descriptions[gene_name])
            elif gene_name in centroid_gene_map:
                if centroid_gene_map[gene_name] in gene_descriptions:
                    gene_description_vec.append(gene_descriptions[centroid_gene_map[gene_name]])
            else:
                gene_description_vec.append("")        

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


        num_sites_processed += 1




    snp_file.close()
    
    if (len(alts_total.shape) == 1):
        alts_total = np.reshape(alts_total, (1, -1))
        depths_total = np.reshape(depths_total, (1, -1))

    ## creating dataframe with alts, depths, and freqs for each sample
    alts_df = pd.DataFrame(data = alts_total, columns = desired_samples)
    depths_df = pd.DataFrame(data = depths_total, columns = desired_samples)
    alts_df['contig'] = chrom_vec
    alts_df['site_pos'] = location_vec
    alts_df['gene'] = gene_vec
    alts_df['species'] = species
    alts_df['gene_description'] = gene_description_vec

    depths_df['contig'] = chrom_vec
    depths_df['site_pos'] = location_vec
    depths_df['gene'] = gene_vec
    depths_df['species'] = species
    depths_df['gene_description'] = gene_description_vec


    alts_df = alts_df.melt(var_name='sample', value_name='alt', id_vars = ['species','gene','gene_description','contig', 'site_pos'])
    depths_df = depths_df.melt(var_name='sample', value_name='depth', id_vars = ['species','gene','gene_description','contig', 'site_pos'])
    freq_df = pd.merge(alts_df, depths_df, on=['species','sample', 'gene','gene_description','site_pos', 'contig'])
    freq_df['allele_frequency'] = freq_df['alt']/freq_df['depth']
    
    # Append freq_df to full_freq_df
    full_freq_df = full_freq_df.append(freq_df, ignore_index=True)

    sys.stderr.write("Annotating dataframe with sample data\n")

# Annotating sample data

full_freq_df['mouse'] = full_freq_df['sample'].apply(lambda sample: extract_mouse_number(sample))
full_freq_df['diet'] = full_freq_df['sample'].apply(lambda sample: extract_diet(sample))
full_freq_df['cage'] = full_freq_df['sample'].apply(lambda sample: extract_cage(sample))
full_freq_df['gut_site'] = full_freq_df['sample'].apply(lambda sample: extract_gut_site(sample))
full_freq_df['gut_region'] = full_freq_df['sample'].apply(lambda sample: extract_region(sample))

if only_haploid:
    output_path = "%s/evolutionary_changes/SNV_freqs_haploid.txt.bz2" % (config.project_folder, str(min_coverage))
else:
    output_path = "%s/evolutionary_changes/SNV_freqs.txt.bz2" % (config.project_folder)
full_freq_df.to_csv(output_path, sep = "\t")


