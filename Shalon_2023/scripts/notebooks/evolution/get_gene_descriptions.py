import numpy as np
import pandas as pd
import os
import statsmodels.api as sm
from statsmodels.formula.api import ols
import string
import itertools
from datetime import datetime
import random


# plotting functions

import matplotlib.pyplot as plt
import seaborn as sns


# predefined functions
os.sys.path.append('/u/project/ngarud/michaelw/microbiome_evolution/py3.8/microbiome_evolution_SHALON/postprocessing_scripts/')
from parse_midas_data import *
from diversity_utils import *
from calculate_intersample_changes import *
from core_gene_utils import *
from parse_patric import *
import parse_patric

# load species list
species_list_path = "/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/metadata/species_snps.txt"

with open(species_list_path, 'r') as file:
    species_list = file.readlines()
species_list = [species.strip() for species in species_list]


# Load maps
sample_metadata_map = parse_sample_metadata_map()
sample_list = list(sample_metadata_map.keys())
subject_sample_map = parse_subject_sample_map()


# path to output
haploid_summary_path = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/evolutionary_changes/QP_status.csv"
generate = False
save = True
if  os.path.exists(haploid_summary_path) and (not generate):
    print("Loading preexisting QP summary file")
    samples_df = pd.read_csv(haploid_summary_path, sep = ",")

    haploid_samples = samples_df[samples_df['haploid']]

else:

    print("Creating QP summary file from scratch.")
    # load species list
    species_list_path = "/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/metadata/species_snps.txt"

    with open(species_list_path, 'r') as file:
        species_list = file.readlines()
    species_list = [species.strip() for species in species_list]

    metadata_map = parse_midas_data.parse_sample_metadata_map()
    samples_df = pd.DataFrame(columns = ['accession_id', 'species', 'haploid', 'subject_id'])

    accession_id_vec = []
    species_vec = []
    haploid_vec = []
    subject_id_vec = []

    print("\nCalculating lineage structure.")

    counter = 0
    no_of_species = len(species_list)
    for species in species_list:
        counter += 1
        print("Processing " + species + " ({}/{})".format(counter, no_of_species))
        high_coverage_samples = diversity_utils.calculate_highcoverage_samples(species)
        haploid_samples = diversity_utils.calculate_haploid_samples(species, use_HMP_freqs = True)
        
        haploid_boolean_vec = [True if sample in haploid_samples else False for sample in high_coverage_samples]
        subject_id_vec_temp = [metadata_map[sample][0] for sample in high_coverage_samples]
        
        species_vec = species_vec + [species]*len(high_coverage_samples)
        accession_id_vec = accession_id_vec + list(high_coverage_samples)
        haploid_vec = haploid_vec + haploid_boolean_vec
        subject_id_vec = subject_id_vec + subject_id_vec_temp

    samples_df['accession_id'] = accession_id_vec
    samples_df['species'] = species_vec
    samples_df['haploid'] = haploid_vec
    samples_df['subject_id'] = subject_id_vec

    if save:
        samples_df.to_csv(haploid_summary_path, sep = ",", index = False)



# path to output
change_summary_path = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/evolutionary_changes/change_summary_Full.csv"
generate = False
save = True

if os.path.exists(change_summary_path) and (not generate):
    print("Loading preexisting change summary file")
    change_summary_df = pd.read_csv(change_summary_path, sep = ",")
    
else:

    haploid_samples = samples_df[samples_df['haploid']]
        
    change_summary_array = []
    haploid_species = haploid_samples.species.unique()

    print("\nCalculating change rate.")

    for i,species in enumerate(haploid_species):
        print("Processing %s (%d/%d)" % (species, i+1, len(haploid_species)))
        intersample_change_map = load_intersample_change_map(species)
        samples = haploid_samples[haploid_samples['species'] == species]['accession_id'].unique()
        sample_pairs = list(itertools.combinations(samples, 2))
        for sample_pair in sample_pairs:
            if (sample_pair[1], sample_pair[0]) in intersample_change_map:
                sample_pair = (sample_pair[1], sample_pair[0])

            if sample_pair not in intersample_change_map:
                print("%s not in intersample change map for %s" % (str(sample_pair), species))
                continue

            opportunities = intersample_change_map[sample_pair]['snps'][0]
            gene_opportunities = intersample_change_map[sample_pair]['genes'][0]
            snv_changes = len(intersample_change_map[sample_pair]['snps'][2])
            gene_changes = len(intersample_change_map[sample_pair]['genes'][2])
            rate_of_change = snv_changes/opportunities
            gene_rate_of_change = gene_changes/gene_opportunities

            change_summary_array.append([species, sample_pair[0], sample_pair[1], snv_changes, opportunities, rate_of_change, gene_changes, gene_opportunities, gene_rate_of_change])
            

    change_summary_df = pd.DataFrame(change_summary_array, columns=['species', 'sample_1', 'sample_2', 'snv_changes', 'opportunities', 'rate_of_change', 'gene_changes', 'gene_opportunities', 'gene_rate_of_change'])

    # annotate
    print("\n\nAnnotating.")
    change_summary_df['subject_1'] = change_summary_df['sample_1'].apply(lambda x: sample_metadata_map[x][0])
    change_summary_df['subject_2'] = change_summary_df['sample_2'].apply(lambda x: sample_metadata_map[x][0])
    change_summary_df['device_type_1'] = change_summary_df['sample_1'].apply(lambda x: sample_metadata_map[x][2])
    change_summary_df['device_type_2'] = change_summary_df['sample_2'].apply(lambda x: sample_metadata_map[x][2])
    timestamp_format = '%Y-%m-%dT%H:%M:%SZ'
    change_summary_df['day_1'] = change_summary_df[['sample_1','device_type_1']].apply(lambda row: datetime.strptime(sample_metadata_map[row['sample_1']][5], timestamp_format).strftime('%Y-%m-%d') if row['device_type_1'] == "Stool" or row['device_type_1'] == "Saliva" else datetime.strptime(sample_metadata_map[row['sample_1']][3], timestamp_format).strftime('%Y-%m-%d'), axis = 1)
    change_summary_df['day_2'] = change_summary_df[['sample_2','device_type_2']].apply(lambda row: datetime.strptime(sample_metadata_map[row['sample_2']][5], timestamp_format).strftime('%Y-%m-%d') if row['device_type_2'] == "Stool" or row['device_type_2'] == "Saliva" else datetime.strptime(sample_metadata_map[row['sample_2']][3], timestamp_format).strftime('%Y-%m-%d'), axis = 1)
    change_summary_df['time_1'] = change_summary_df[['sample_1','device_type_1']].apply(lambda row: datetime.strptime(sample_metadata_map[row['sample_1']][5], timestamp_format).strftime('%H:%M:%S') if row['device_type_1'] == "Stool" or row['device_type_1'] == "Saliva" else datetime.strptime(sample_metadata_map[row['sample_1']][3], timestamp_format).strftime('%H:%M:%S'), axis = 1)
    change_summary_df['time_2'] = change_summary_df[['sample_2','device_type_2']].apply(lambda row: datetime.strptime(sample_metadata_map[row['sample_2']][5], timestamp_format).strftime('%H:%M:%S') if row['device_type_2'] == "Stool" or row['device_type_2'] == "Saliva" else datetime.strptime(sample_metadata_map[row['sample_2']][3], timestamp_format).strftime('%H:%M:%S'), axis = 1)

    # Label as within timepoint or between time point
    change_summary_df['timepoint_orientation'] = change_summary_df.apply(lambda row: "Within timepoint" if row['day_1'] == row['day_2'] and row['time_1'] == row['time_2'] else "Between timepoint", axis = 1)
    change_summary_df['device_orientation'] = change_summary_df.apply(lambda row: "Same device" if row['device_type_1'] == row['device_type_2'] else "Different device", axis = 1)
    change_summary_df['datetime_1'] = pd.to_datetime(change_summary_df['day_1'] + ' ' + change_summary_df['time_1'])
    change_summary_df['datetime_2'] = pd.to_datetime(change_summary_df['day_2'] + ' ' + change_summary_df['time_2'])
    change_summary_df['time_difference_hours'] = (change_summary_df['datetime_2'] - change_summary_df['datetime_1']).dt.total_seconds() / 3600

    if save:
        change_summary_df.to_csv(change_summary_path, sep = ",", index=False)
        

    
    # Filter for <= 20 SNV changes
change_summary_df_full = change_summary_df.copy()
change_summary_df = change_summary_df[change_summary_df['snv_changes'] <= 20]
# Filter for within host
change_summary_df = change_summary_df[change_summary_df['subject_1'] == change_summary_df['subject_2']]
# Filter for between capsule
change_summary_df = change_summary_df[(change_summary_df['device_type_1'] != "Stool") &
                                      (change_summary_df['device_type_1'] != "Saliva") &
                                      (change_summary_df['device_type_2'] != "Stool") &
                                      (change_summary_df['device_type_2'] != "Saliva")]
# Filter for not the same device in the same timepoint (which only happens in subject 1)
change_summary_df = change_summary_df[~((change_summary_df['timepoint_orientation'] == "Within timepoint") &
                                      (change_summary_df['device_orientation'] == "Same device"))]
# Reset index
change_summary_df = change_summary_df.reset_index(drop=True)


# species-subject pair list
species_list = change_summary_df[change_summary_df.snv_changes > 0][['species','sample_1','sample_2']].drop_duplicates().reset_index(drop = True)

# Getting gene descriptions
species_vec = []
contig_vec = []
site_pos_vec = []
contig_site_pos_vec = []
variant_type_vec = []
gene_id_vec = []
gene_description_vec = []

last_species = ""


for i,row in species_list.iterrows():
    
    species = row['species']
    sample_1 = row['sample_1']
    sample_2 = row['sample_2']
    
    sys.stderr.write("Processing %s in %s - %s sample pair (%d / %d species)\n" % (species, sample_1, sample_2, i+1, len(species_list)))

    if species != last_species:
        last_species = species
        intersample_change_map = load_intersample_change_map(species)
    
    if (sample_1,sample_2) in intersample_change_map:
        sample_pair = (sample_1,sample_2)
        snvs = intersample_change_map[sample_pair]['snps'][2]
    elif (sample_2,sample_1) in intersample_change_map:
        sample_pair = (sample_2,sample_1)
        snvs = intersample_change_map[sample_pair]['snps'][2]
    else:
        sys.stderr.write("%s not in intersample change map for %s\n" % (str(sample_pair), species))
        continue


    for snv in snvs:
        gene_id = snv[0]
        contig = snv[1]
        site_pos = snv[2]
        variant_type = snv[3]
        contig_site_pos = contig + "|" + str(site_pos)
        if contig_site_pos not in contig_site_pos_vec:

            species_vec.append(species)
            contig_vec.append(contig)
            site_pos_vec.append(site_pos)
            contig_site_pos_vec.append(contig_site_pos)
            gene_id_vec.append(gene_id)
            variant_type_vec.append(variant_type)

            genome_ids = parse_midas_data.get_ref_genome_ids(species)
            non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species)
            gene_descriptions = parse_patric.load_patric_gene_descriptions(genome_ids, non_shared_genes)
            centroid_gene_map = parse_midas_data.load_centroid_gene_map(species)

            if gene_id in gene_descriptions:
                gene_description_vec.append(gene_descriptions[gene_id])
            elif gene_id in centroid_gene_map:
                if centroid_gene_map[gene_id] in gene_descriptions:
                    gene_description_vec.append(gene_descriptions[centroid_gene_map[gene_id]])
            else:
                gene_description_vec.append("")

snv_gene_descriptions_df = pd.DataFrame({
    "species": species_vec,
    "contig": contig_vec,
    "site_pos": site_pos_vec,
    "contig_site_pos": contig_site_pos_vec,
    "gene_id": gene_id_vec,
    "variant_type": variant_type_vec,
    "gene_description": gene_description_vec,
})

# saving
out_path = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/evolutionary_changes/SNV_gene_descriptions.tsv"
snv_gene_descriptions_df.to_csv(out_path, sep = "\t", index = False)

