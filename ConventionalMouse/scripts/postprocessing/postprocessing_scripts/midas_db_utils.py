import config
import gzip
import os
import os.path
import species_utils

###############################################################################
#
# Loads list of genes in the reference genome used by MIDAS for a given species
#
###############################################################################
def load_reference_genes(desired_species_name): # MW 07/14/2025: updated to handle MIDASdb-gtdb
    reference_genome_id = species_utils.load_reference_genome_id(desired_species_name) # MW 09/03/2025: loads the reference genome annotations (used to load random genome)
    features_file = open("%sgene_annotations/%s/%s/%s.genes" % (config.midas_directory, desired_species_name, reference_genome_id, reference_genome_id), 'r') 
    
    features_file.readline() # header
    reference_genes = []
    for line in features_file:
        items = line.split() 
        gene_name = items[0].strip()
        reference_genes.append(gene_name)
    features_file.close()    
    
    return set(reference_genes)

def get_pangenome_map(species_name):
    
    gene_info_filename = '%span_genomes/%s/gene_info.txt.gz' % (config.midas_directory, species_name)
    file = gzip.open(gene_info_filename, 'r')
    file.readline() # header
    
    pangenome_map = {}
    
    for line in file:
        items = line.decode("utf-8").split("\t") # MW 07/09/2025: decode bytes like object
        gene_id = items[0].strip()
        genome_id = items[1].strip()
        centroid_99 = items[2].strip()
        centroid_95 = items[3].strip()
        
        if genome_id not in pangenome_map:
            pangenome_map[genome_id] = {}
        
        pangenome_map[genome_id][gene_id] = (centroid_99, centroid_95)
        
    file.close()
    return pangenome_map
    
def get_number_of_genomes(species_name):
    
    return len(get_pangenome_map(species_name))

def parse_species_list():
    
    species_directories = os.listdir(config.midas_directory+"/pan_genomes")
    
    species_names = []
    for potential_species_name in species_directories:
        if not potential_species_name.startswith('.'):
            species_names.append(potential_species_name)
    
    return species_names
    

####
#
# The gene_ids in the pangenome list are the centroids of gene clusters.
# Sometimes the gene in the reference genome is not chosen as the centroid.
# This function creates a map between pangenome_centroids and genes in 
# reference genome (if it exists)
#
###
def load_centroid_gene_map(desired_species_name=None): # MW 07/09/2025: updated to work with MIDASdb-gtdb
    
    if desired_species_name==None:
        import parse_midas_data
        desired_speciess = parse_midas_data.parse_good_species_list()
    else:
        desired_speciess = [desired_species_name]
    
    for desired_species_name in desired_speciess:
        # First load reference genes
        reference_genes = load_reference_genes(desired_species_name)
    
        gene_info_file = open("%spangenomes/%s/clusters_99_info.tsv" % (config.midas_directory, desired_species_name), 'r') # MW 07/08/2025: changed file path
    
        gene_info_file.readline() # header
    
        centroid_gene_map = {}
    
        for line in gene_info_file:
        
            items = line.split("\t") 
            gene_id = items[0].strip()
            centroid_id = items[1].strip() # MW 07/09/2025: changed index from 3 to 1 (difference in dataframe structure)
        
            if centroid_id not in centroid_gene_map:
                centroid_gene_map[centroid_id] = centroid_id
            
            if (gene_id in reference_genes) and (centroid_id not in reference_genes):
                centroid_gene_map[centroid_id] = gene_id
            
        
        gene_info_file.close()
    
    return centroid_gene_map
    
    
def parse_midas_shared_genes(desired_species):
    
    midas_shared_genes = set()
    
    # get list 
    centroid_gene_map = load_centroid_gene_map(desired_species)
    
    midas_db_shared_gene_filename = (config.midas_directory+"cross_species_centroids.txt.gz") 
    file = gzip.open(midas_db_shared_gene_filename,"r")
    for line in file:
        items = line.decode("utf-8").split() # MW 07/09/2025: decode bytes like object
        big_centroid = items[0]
        midas_shared_genes.add(big_centroid.strip())
        other_centroids = items[1].split(",") 
        for centroid in other_centroids:
            stripped_centroid = centroid.strip()
            if centroid in centroid_gene_map:
                midas_shared_genes.add(centroid_gene_map[stripped_centroid])
            
    return midas_shared_genes
               
    
if __name__=='__main__':
    
    import parse_midas_data
    good_species_list = parse_midas_data.parse_good_species_list()
    for species_name in good_species_list:
        print(species_name)
        print(get_number_of_genomes(species_name))
    