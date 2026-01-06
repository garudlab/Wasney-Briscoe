import numpy
import sys
import os.path 
import pandas as pd
from collections import Counter

import config


###############################################################################
#
# Set up default source and output directories
#
###############################################################################

midasdb_metadata_directory = "%s%s" % (config.midas_directory,"metadata.tsv")
species_list_directory = "%s%s" % (config.metadata_directory, "species_snps.txt")

###############################################################################
#
# Load species list
#
###############################################################################

def load_species_list():
    with open(species_list_directory, "r") as file:
        species_list = [line.strip() for line in file if line.strip()]

    return(species_list)


###############################################################################
#
# Load species maps
#
###############################################################################

def parse_species_code_maps(map = "all"):

    midasdb_metadata_path = "%s%s" % (config.midas_directory,"metadata.tsv") 

    species_code_name_map = {}

    species_code_list = load_species_list()

    with open(midasdb_metadata_path, "r") as file:
        header = file.readline().strip().split("\t") 
        species_id_idx = header.index("species_id")
        gtdb_species_idx = header.index("species")

        for line in file:
            items = line.strip().split("\t")
            species_id = items[species_id_idx].strip()
            if species_id not in species_code_list:
                continue
            species_name = items[gtdb_species_idx].strip().strip("s__")
            species_code_name_map[species_id] = species_name

    # Find duplicates using collections.Counter
    name_values = list(species_code_name_map.values())
    duplicates = [item for item, count in Counter(name_values).items() if count > 1]

    # Update maps based on duplicates
    if len(duplicates) == 0:
        species_name_code_map = {v: k for k, v in species_code_name_map.items()}
    else:
        species_code_name_map = {
            k: f"{v} ({k})" if v in duplicates else v
            for k, v in species_code_name_map.items()
        }
        species_name_code_map = {v: k for k, v in species_code_name_map.items()}
    
    if map == "all":
        return(species_code_name_map, species_name_code_map)
    elif map == "code-name":
        return(species_code_name_map)
    elif map == "name-code":
        return(species_name_code_map)
    else:
        raise ValueError("map options incude 'all' (default), 'code-name', or 'name-code'.")

###############################################################################
#
# Load reference id
#
###############################################################################

def load_reference_genome_id(species_name):
    genomes_df_path = "%s%s" % (config.midas_directory, "genomes.tsv")
    genomes_df = pd.read_csv(genomes_df_path, sep = "\t")
    genomes_df["species"] = genomes_df["species"].astype(str)
    reference_genome_id = genomes_df.loc[(genomes_df.species == species_name) & (genomes_df.genome_is_representative == 1), "genome"].values[0]

    return(reference_genome_id)
