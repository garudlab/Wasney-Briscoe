import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob

import pandas as pd
import core_gene_utils
import parse_midas_data
import parse_patric
import config
import sys

def return_gene_descriptions(species):
    

    genome_ids=parse_midas_data.get_ref_genome_ids(species)
    non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species)
    non_shared_genes = [n.decode("utf-8") for n in non_shared_genes]
    gene_descriptions = parse_patric.load_patric_gene_descriptions(genome_ids, non_shared_genes)
    all_gene_cats = list(set(gene_descriptions.values()))
    gene_descriptions = pd.Series(gene_descriptions)

    return(gene_descriptions)

## returns annotations 
def read_sites(species):
    
    snps_directory = "/u/project/ngarud/Garud_lab/HumanizedMouse/merged_midas_output/snps"
    
    df_sites = pd.read_csv(f"{snps_directory}/{species}/snps_info.txt.bz2",sep="\t",index_col=0,na_values="NaN")
   
    df_sites["contig"] = [d.split("|")[0] for d in df_sites.index]
    df_sites.index = [d.split("|")[1] for d in df_sites.index]
    
    df_sites["gene_id"] = df_sites["gene_id"].fillna("non coding")
    gene_ids = df_sites["gene_id"].values

    gene_breaks = [0]

    gc = gene_ids[0]
    unq_genes = [gc]
    unq_cont = [df_sites["contig"][0]]

    for i,g in enumerate(gene_ids):
        if g is not gc:
            gene_breaks.append(i)
            gc = g
            unq_genes.append(gc)
            unq_cont.append(df_sites["contig"][i])

    gene_breaks = np.array(gene_breaks)       
    gene_lengths = gene_breaks[1:] - gene_breaks[:-1] 

    df_sites.index.set_names("site_pos",inplace=True)
    
    df_sites.set_index('gene_id', append=True, inplace=True)
    df_sites.set_index('contig', append=True, inplace=True)

    df_sites = df_sites.reorder_levels(["contig",'gene_id', 'site_pos'])
    
    level_to_change = 2
    df_sites.index = df_sites.index.set_levels(df_sites.index.levels[level_to_change].astype(int), level=level_to_change)

    return(df_sites)

if __name__ == '__main__':

    good_species = config.good_species
    
    for species in good_species:
        
        sys.stderr.write(f"Processing {species}\n\n")
        
       # try:
            
        gene_descriptions = return_gene_descriptions(species)
        snps_info = read_sites(species)
        gene_descriptions.to_pickle(f"/u/project/ngarud/rwolff/mouse_sites/gene_descriptions/{species}_gene_descriptions.pkl")
        snps_info.to_pickle(f"/u/project/ngarud/rwolff/mouse_sites/snps_info/{species}_snps_info.pkl")

        sys.stderr.write(f"\t{species} completed\n\n")

#         except:
            
#             sys.stderr.write(f"\tError: {species}\n\n")
            