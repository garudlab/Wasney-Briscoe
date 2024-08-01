import pandas as pd
import core_gene_utils
import parse_midas_data
import parse_patric

def return_gene_descriptions(species):
    

    genome_ids=parse_midas_data.get_ref_genome_ids(species)
    non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species)
    non_shared_genes = [n.decode("utf-8") for n in non_shared_genes]
    gene_descriptions = parse_patric.load_patric_gene_descriptions(genome_ids, non_shared_genes)
    all_gene_cats = list(set(gene_descriptions.values()))
    gene_descriptions = pd.Series(gene_descriptions)

    return(gene_descriptions)
