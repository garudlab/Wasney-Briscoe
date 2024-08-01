import numpy as np
import pandas as pd
import config
import sys
import parse_midas_data
import itertools
import diversity_utils
import glob
import matplotlib.pyplot as plt
import calculate_temporal_changes as ct
import diversity_utils 
from numba import njit,jit
import parse_HMP_data
import figure_utils as fu
import os

good_species = config.good_species
analysis_dir = config.analysis_directory
min_depth = config.min_depth

def read_sites(species):
    
    df_sites = pd.read_csv(f"{config.snps_directory}/{species}/snps_info.txt.bz2",sep="\t",index_col=0,na_values="NaN")

    df_sites["contig"] = [d.split("|")[0] for d in df_sites.index]
    df_sites.index = [d.split("|")[1] for d in df_sites.index]

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
    
    return(df_sites,gene_lengths,unq_genes,unq_cont)

if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--species',
                        help="Species to generate haplotypes for",
                        type=str)

    parser.add_argument('--min_depth',
                        help="Minimum depth to use to call a variant at a site in a sample",
                        type=int,
                        default=20)

    parser.add_argument('--min_MAF',
                        help="Minimum allele frequency to be considered a polymorphic site. Note: allele frequency is determined as percentage of total number of reads across all hosts supporting an allele. Large fluctuations in allele frequency across hosts can bias this.",
                        type=float,
                        default=0.01)

#     parser.add_argument('--no_singletons',
#                     help="If true, do not consider any site only p . Default is false.",
#                     action='store_true'
#                     type=bool)

    parser.add_argument('--f_star',
                    help="Within host allele frequency to call consensus",
                    type=float,
                    default=0.2)    
    
    args = parser.parse_args()
    
    species = args.species    
    min_depth = args.min_depth
    min_MAF = args.min_MAF
#    no_singletons = args.no_singletons
    f_star = args.f_star

    sys.stderr.write(f"Processing {species} \n")
    
    os.makedirs(f"{config.haplotype_directory}/{species}",exist_ok=True)
    
    snps_dir = "%ssnps/%s" % (config.data_directory,species)
    snps_summary = pd.read_csv("%s/snps_summary.txt" % snps_dir,sep="\t",index_col=0)
    L = snps_summary["covered_bases"]

    samples_host = list(pd.read_csv("%s/snps_depth.txt.bz2" % snps_dir,sep="\t",index_col=0, nrows=0))
    samples_tuples = list(itertools.combinations(samples_host, 2)) + [(s1,s1) for s1 in samples_host]


    reader=True

    haploid_samples = diversity_utils.calculate_haploid_samples(species)

    subject_map = parse_HMP_data.parse_subject_sample_map()
    sample_map = {}
    for key, item in subject_map.items():
        for subkey in item:
            sample_map[subkey] = key

    sample_map = pd.DataFrame(pd.Series(sample_map),columns=["host_id"])

    ## initialize chunk readers for sample depth and allele frequency 
    df_depth_reader = pd.read_csv("%s/snps_depth.txt.bz2" % snps_dir,sep="\t",index_col=0, iterator=True,low_memory=False)
    df_refreq_reader = pd.read_csv("%s/snps_ref_freq.txt.bz2" % snps_dir,sep="\t",index_col=0, iterator=True,low_memory=False)

    ## reads header. chunking can now proceed on data in files
    df_depth_header = df_depth_reader.get_chunk(0)
    df_refreq_header = df_refreq_reader.get_chunk(0)

    ## read snps_info file, which contains cohort-wide information about allele frequency, 
    ## as well as information about sites (contig/gene, S/NS etc)
    df_sites,gene_lengths,unq_genes,unq_cont = read_sites(species)
    
    ## frequency of each nucleotide at each site
    atcg = df_sites["allele_props"].dropna().str.split("|")
    atcg = pd.DataFrame([[elem[2:] for elem in e] for e in atcg],index=atcg.index,columns=["a","c","t","g"])
    atcg = atcg.astype(float)
    
    ## pulls S/NS status for each nucleotide mutation at each site relative to ref state
    syn_non = df_sites["snps"].loc[df_sites.index]
    syn_non = syn_non.dropna().str.split("|")
    syn_non = syn_non.loc[[len(s) > 1 for s in syn_non]]
    syn_non = pd.DataFrame([[elem[2:] for elem in e] for e in syn_non],index=syn_non.index,columns=["a","c","t","g"])
 
    ## counts number of alleles at a site w/ frequency > min minor allele frequency
    ## essentially, what counts as a "polymorphic" site
    realized = syn_non.where(atcg > min_MAF)
    num_alleles = realized.T.count()
    ## retrieve sites with exactly two alleles 
    realized_bi = realized.loc[num_alleles.loc[num_alleles == 2].index]
   
    ## Assess if bi-allelic site is SYN (AA preserving) or NS (AA changing) for each variant
    realized_bi_syn = 1*(realized_bi == "SYN")
    realized_bi_ns = 1*(realized_bi == "NS")
    
    ## if exactly one variant is SYN, and the other NS, tag site as NS. If both are S, tag site as S
    ns = (realized_bi_syn.loc[realized_bi_syn.T.sum() == 1]).index
    s = (realized_bi_syn.loc[realized_bi_syn.T.sum() == 2]).index
    
    ## in output file, refer to bi-allelic sites w/ non-syn mutation as "1D" and w/ syn as "4D".
    ## the reason for this designation is that downstream applications in the LD work 
    ## use a computationally efficient numeric code system to sort pairwise LD's  
    ## 1x1 = 1: NS/NS comparison
    ## 4x1 = 1x4 = 4: NS/S comparison (commutativity of multiplication is the key to this system)
    ## 4x4 = 16: S/S comparison
    df_sites.loc[ns,"site_type"] = "1D"
    df_sites.loc[s,"site_type"] = "4D"
     
    ## remove mono-, tri- and quad-allelic sites from further consideration
    df_sites = df_sites.loc[realized_bi.index]
        
    ## remove tri- and quad-allelic sites 
    # realized_tri = realized.loc[num_alleles.loc[num_alleles == 3].index]
    # realized_quad = realized.loc[num_alleles.loc[num_alleles == 4].index]    
    #df_sites = df_sites.drop(realized_tri.index)
    # df_sites = df_sites.drop(realized_quad.index)

    sys.stderr.write(f"Looping over genes \n")
                             
    for i,chunk_size in enumerate(gene_lengths):

        ## read next chunk_size number of lines 
        df_depth = df_depth_reader.get_chunk(chunk_size)
        df_refreq = df_refreq_reader.get_chunk(chunk_size)

        ## if we're in a coding region, create a haplotype
        if unq_genes[i] is not np.nan:

            #sys.stderr.write(f"Processing {unq_genes[i]} \n")

            df_depth.columns = [d[:-1] if d[-1] == "c" else d for d in df_depth.columns]
            df_refreq.columns = [d[:-1] if d[-1] == "c" else d for d in df_refreq.columns]

            df_depth = df_depth[haploid_samples]
            df_refreq = df_refreq[haploid_samples]

            ## initialize haplotype dataframe as alternate allele freq in each sample
            df_haplotypes = 1 - df_refreq.copy()

            ## treat sites with less than min_depth coverage as missing data
            df_haplotypes = df_haplotypes.mask(df_depth < min_depth) 

            ## if alternate is consensus, mark sample as 1. if ref is consensus, mark as 0.
            df_haplotypes = df_haplotypes.mask(df_haplotypes >= .8,1)
            df_haplotypes = df_haplotypes.mask(df_haplotypes <= .2,0)

            ## treat intermediate frequency variants as missing data
            df_depth = df_depth.mask(np.logical_and(df_haplotypes < .8,df_haplotypes > .2))
            df_haplotypes = df_haplotypes.mask(np.logical_and(df_haplotypes < .8,df_haplotypes > .2)) 

            df_haplotypes.index = [d.split("|")[1] for d in df_haplotypes.index]
            df_haplotypes.index.set_names("site_pos",inplace=True)
            df_haplotypes["gene_id"] = df_haplotypes.shape[0]*[unq_genes[i]]
            df_haplotypes["contig"] = df_haplotypes.shape[0]*[unq_cont[i]]
            df_haplotypes.set_index('gene_id', append=True, inplace=True) 
            df_haplotypes.set_index('contig', append=True, inplace=True) 
            df_haplotypes = df_haplotypes.reorder_levels(["contig",'gene_id', 'site_pos'])

            ## cut down to sites marked as bi-allelic
            df_haplotypes = df_haplotypes.loc[[d for d in df_haplotypes.index if d in df_sites.index]]          

            df_haplotypes["site_type"] = df_sites.loc[df_haplotypes.index]["site_type"]
            df_haplotypes.set_index('site_type', append=True, inplace=True)   
            
            df_haplotypes = df_haplotypes.astype(float).astype('Int64')
            
            ## haplotype data format: 
            ## index column 1: contig
            ## index column 2: gene
            ## index column 3: site position (i.e. position of site on contig)
            ## index column 4: site type (i.e. syn vs non-syn, tagged as 4D (syn) or 1D (non-syn)
            ## data: sites x allowed samples: 0 if sample is consensus for ref, 1 if consensus for alt
            df_haplotypes = df_haplotypes.reorder_levels(["contig",'gene_id', 'site_pos','site_type'])      

            ## write gene haplotype to file 
            df_haplotypes.to_csv(f"{config.haplotype_directory}/{species}/{unq_genes[i]}_haplotypes.csv")

        if i%50 == 0:
            sys.stderr.write(f"\n \n \n \n \n \n \n {np.around(100*i/len(gene_lengths),3)}% complete \n \n \n \n \n \n \n")

