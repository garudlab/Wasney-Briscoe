### Packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib
# matplotlib.rc('text', usetex=True)
# plt.rc('text', usetex=True)
# plt.rc('text.latex', preamble=r'\usepackage{amsmath}') 
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
# from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
from mpl_toolkits.axes_grid1.inset_locator import (
    inset_axes, InsetPosition, mark_inset
)


import scipy.stats
import figure_utils as fu

from numba import njit 

import matplotlib.pyplot as plt
import pickle
import pandas as pd
import config
import numpy
import random as rand

from random import randint,sample
from math import log

import sys
import os 
from scipy.spatial.distance import pdist,squareform

# microbiome functions
sys.path.insert(0, "/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/postprocessing/postprocessing_scripts/")
from parse_midas_data import *


hap_cmap = ListedColormap(['grey', 'red', 'black', 'black','blue'], 'indexed')

### Functions

## useful utility for quickly returning the upper triangle of a 2-d array as a 1-d array
def take_triu(df):
    
    N = df.shape[0]
    p=np.triu_indices(N,k=1)
    
    return(df[p])

## returns annotations 
def read_sites(species):
    
    snps_directory = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/merged_data/snps"
    
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

def read_haplotypes(species,good_samples=None,clade_control=False):
    
    output_dir = "/u/scratch/r/rwolff/LD/HMP"
    haplotype_directory = "%s/haplotypes" % output_dir
    
    hap_files = glob.glob(f"{haplotype_directory}/{species}/*_haplotypes.csv")
    
    idx_cols = ["contig","gene_id","site_pos","site_type"]
    all_haps = []

    for hap in hap_files:
        
        if good_samples is not None and clade_control==False:
            df = pd.read_csv(hap,index_col=[0,1,2,3],usecols = idx_cols + list(good_samples))
            
        elif good_samples is None and clade_control != False:
            df = pd.read_csv(hap,index_col=[0,1,2,3],usecols = idx_cols + list(clade_utils.load_largest_clade(species)))
                  
        else:
            df = pd.read_csv(hap,index_col=[0,1,2,3])

        all_haps.append(df)

    df = pd.concat(all_haps)
    
    return(df)
    
def get_cmap(n, name='Set3_r'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

### use of numba/njit below *substantially* increases computational efficiency

## calculate distances for forward polarization 
@njit
def D_mat_fun1(num,F,D,D_mat):   

    for k in range(num - 1):
        
        O = np.zeros(num)
        
        di = D[k]
        fi = F[k]
        
        for i in range(num - k - 1):

            j = i + k + 1

            fj = F[j]
            dj = D[j]

            O[j] = 2*np.nanmean((di + dj)*((fi - fj)**2)/((fi + fj)*(1 - fi + 1 - fj)))        
        
        D_mat[k] = O
    
    return D_mat

## calculate distances for reverse polarization
@njit
def D_mat_fun2(num,F,D,D_mat_in):   

    for k in range(num - 1):
        
        O = np.zeros(num)       
        di = D[k]
        
        fi = 1-F[k]
        
        for i in range(num - k - 1):

            j = i + k + 1
            
            fj = F[j]
            dj = D[j]

            O[j] = 2*np.nanmean((di + dj)*((fi - fj)**2)/((fi + fj)*(1 - fi + 1 - fj)))        
        
        D_mat_in[k] = O
    
    return D_mat_in

def return_clus(D_mat_close,Fs_sub, co_cluster_pct=0.25):
    D_mat_close_sorted_sum = D_mat_close.sum().sort_values()
    desired_idx = D_mat_close_sorted_sum.index[-1]
    clus_idxs = D_mat_close.loc[D_mat_close[desired_idx]].index
    
    ### only return indices which are co-clustered w/ at least .25 of other points
    idxtrue = (D_mat_close.loc[clus_idxs,clus_idxs].T.mean() > co_cluster_pct)
    idxtrue = idxtrue[idxtrue].index
    clus_idxs = idxtrue
    clus = Fs_sub.loc[clus_idxs]
     
    return clus,clus_idxs

def return_clus_sub25(D_mat_close,Fs_sub):
    D_mat_close_sorted_sum = D_mat_close.sum().sort_values()
    desired_idx = D_mat_close_sorted_sum.index[-1]
    clus_idxs = D_mat_close.loc[D_mat_close[desired_idx]].index
    
    ### only return indices which are co-clustered w/ at least .25 of other points
    idxtrue = (D_mat_close.loc[clus_idxs,clus_idxs].T.mean() > .25)
    idxtrue = idxtrue[idxtrue].index
    clus_idxs = idxtrue
    clus = Fs_sub.loc[clus_idxs]
     
    return clus,clus_idxs

def drop_clus_idxs(D_mat_close,clus_idxs):
    D_mat_close_out = D_mat_close.drop(clus_idxs).drop(clus_idxs,axis=1)
    return D_mat_close_out

def polarize_clus(clus,clus_idxs,D_mat_1,D_mat_2):
    
    ## polarize whole cluster based on polarization of first cluster element
    clus_to_pol = clus_idxs[np.where(D_mat_1.loc[clus_idxs[:1],clus_idxs] >= D_mat_2.loc[clus_idxs[:1],clus_idxs])[1]]
    pol_clus = 1 - clus.loc[clus_to_pol]
    clus_non_pol = clus_idxs[np.where(D_mat_1.loc[clus_idxs[:1],clus_idxs] < D_mat_2.loc[clus_idxs[:1],clus_idxs])[1]]
    non_pol_plus = clus.loc[clus_non_pol]
    clus_pol = pd.concat([pol_clus,non_pol_plus],ignore_index=True)
    
    return(clus_pol)

@njit
def symmetrize(D_mat):
    for i in range(D_mat.shape[0]-1):
        for j in range(i,D_mat.shape[0]):
            D_mat[j][i] = D_mat[i][j]
    return(D_mat)

## non-looped version of distance calculation
def calc_dis(di,dj,fi,fj):
    
    return(2*np.nanmean((di + dj)*((fi - fj)**2)/((fi + fj)*(1 - fi + 1 - fj)),dtype='float32'))
    
def return_FAD(species,min_coverage=10,min_sample_coverage=5,poly_cov_frac=1/5, calculate_poly_cov_frac=False, read_support = True):
    #MW 07/01/24: read_support indicates that you need at least 2 reads supporting a polymporphism


    strainfinder_dir = "%sinput" % (config.strain_phasing_directory)

    snp_alignment = pd.read_pickle("%s/%s/%s.strainfinder.p" %  (strainfinder_dir ,species, species))
    samples = pd.read_pickle("%s/%s/%s.strainfinder.samples.p" % (strainfinder_dir ,species, species))
    snp_locations = pd.read_pickle("%s/%s/%s.strainfinder.locations.p" % (strainfinder_dir,species,species))
    sample_metadata_map = parse_sample_metadata_map()

    cluster_As = []
    cluster_Ds = []
    for snp_idx in range(0,snp_alignment.shape[1]):

        Ds = snp_alignment[:,snp_idx,:].sum(axis=1)
        As = snp_alignment[:,snp_idx,0]
        As = np.reshape(As, (1,len(As)))
        Ds = np.reshape(Ds, (1,len(Ds)))

        cluster_As.append(As[0])
        cluster_Ds.append(Ds[0])

    cluster_As = np.array(cluster_As)
    cluster_Ds = np.array(cluster_Ds)

    As = pd.DataFrame(cluster_As,columns=samples,index=snp_locations)
    Ds = pd.DataFrame(cluster_Ds,columns=samples,index=snp_locations)

    F = As/Ds

    Ass = As
    Dss = Ds.loc[Ass.index]
    Ass = Ass.mask(Dss < min_coverage)
    
    Fs = Ass/Dss

    samps = Dss.mean() > min_sample_coverage
    samps = samps[samps].index
    Dss = Dss[samps]
    Ass = Ass[samps]
    Fs = Fs[samps]
    
    if calculate_poly_cov_frac: #MW: added this to calculate poly_cov_frac based on the number of samples present after coverage filter
        poly_cov_frac = 1/len(samps)
    
    #MW EDIT: I filter out sites that are fixed for EITHER the reference or the alternate allele
    if read_support:
        sys.stderr.write('Using a read support of 4 for each polymorphism.\n')
        Fs = Fs.loc[Fs.mask(Ass < 4).mask(Ass > Dss-4).mask(Dss < min_coverage).notna().T.sum() > int(Fs.shape[1]*poly_cov_frac)]        
    else:
        sys.stderr.write('Not using read support.\n')
        Fs = Fs.loc[Fs.mask(Ass == 0).mask(Ass == Dss).mask(Dss < min_coverage).notna().T.sum() > int(Fs.shape[1]*poly_cov_frac)]
    #MW EDIT: This is where i would disaggregate coverage requirements and polymorphisms
    #Fs = Fs.loc[Fs.mask(Ass == 0).mask(Ass == Dss).mask(Dss < min_coverage).notna().T.sum() > int(Fs.shape[1]*poly_cov_frac)]

    
    Ass = Ass.loc[Fs.index]
    Dss = Dss.loc[Fs.index]

    mnum = ["Mouse %s" % str(sample_metadata_map[f]["mouse"]) for f in Fs.columns]
    msite = [sample_metadata_map[f]["location"] for f in Fs.columns]
    mcage = ["Cage %s" % str(sample_metadata_map[f]["cage"]) for f in Fs.columns]
    Fs = Fs.T
    Fs["mouse_number"] = mnum
    Fs.set_index('mouse_number', append=True, inplace=True)
    Fs["region"] = msite
    Fs.set_index('region', append=True, inplace=True)
    Fs["cage"] = mcage
    Fs.set_index('cage', append=True, inplace=True)

    Fs.index.names = ["sample","mouse_number","region","cage"]
    Fs = Fs.reorder_levels(["cage","mouse_number","region","sample"])
    Fs = Fs.T

    Ass.columns = Fs.columns
    Dss.columns = Fs.columns

    Fs = Fs.sort_index(level="mouse_number",axis=1)
    Dss = Dss.sort_index(level="mouse_number",axis=1)
    Ass = Ass.sort_index(level="mouse_number",axis=1)

    snv_idx = pd.MultiIndex.from_tuples(Fs.index,names=["contig","site_pos","ref/alt"])
    Fs.index = snv_idx
    Fs = Fs.droplevel("ref/alt")

    Ass.index = Fs.index
    Dss.index = Fs.index

    Fs=Fs.sort_index(level=["contig",'site_pos'])
    Ass=Ass.sort_index(level=["contig",'site_pos'])
    Dss=Dss.sort_index(level=["contig",'site_pos'])

    C_list = np.unique(Fs.index.get_level_values("contig"))

    all_site_pos = []
    offset = 0
    for C in C_list:

        all_site_pos.extend(Fs.loc[C].index + offset)

        offset += Fs.loc[C].index[-1]

    Fs['all_site_pos'] = all_site_pos
    Fs.set_index('all_site_pos', append=True, inplace=True)
    Ass['all_site_pos'] = all_site_pos
    Ass.set_index('all_site_pos', append=True, inplace=True)
    Dss['all_site_pos'] = all_site_pos
    Dss.set_index('all_site_pos', append=True, inplace=True)

    return(Fs,Ass,Dss)


## defines an order for each level type
## M1 --> M6
## Upper gut --> lower gut
## Control --> Guar gum
## Co-housing treatment 1 --> co-housing treatment 3
order_dict = {"Mouse 1":0,"Mouse 2":1,"Mouse 3":2,"Mouse 4":3,"Mouse 5":4,"Mouse 6":5,"Mouse 7":6,"Mouse 8":7,
              'Duodenum': 0, 'Jejunum': 1, 'Ileum': 2,"Cecum":3,"Colon":4,
              "Cage 1":0,"Cage 2":1,"Cage 3":2} 

order_dict_reverse = {"Mouse 1":1,"Mouse 2":2,"Mouse 3":3,"Mouse 4":4,"Mouse 5":5,"Mouse 6":6,"Mouse 7":7,"Mouse 8":8,
                      'Duodenum': 1, 'Jejunum': 2, 'Ileum': 3,"Cecum":4,"Colon":5,
                      "Cage 1":1,"Cage 2":2,"Cage 3":3} 

## function sorts our multiindex of frequencies according to whatever key we specify, e.g. mouse_number
## while maintaining the order of subsequent levels according to order_dict
def reorder_sort(df,first_idx,order_dict=order_dict):
    
    reorder=list(df.index.names)
    reorder.remove(first_idx)
    reorder.insert(0, first_idx)
    return df.reorder_levels(reorder).sort_index(key=lambda x: x.map(order_dict))


##########################################################################################
#
#     MW ADDITIONS
#
##########################################################################################

### Annotation functions

def extract_gut_site(sample_name):
    if sample_name[2:4] == "Co":
        return("Colon")
    elif sample_name[2:4] == "Ce":
        return("Cecum")
    elif sample_name[2:3] == "I":
        return("Ileum")
    elif sample_name[2:3] == "J":
        return("Jejunum")
    elif sample_name[2:3] == "D":
        return("Duodenum")
    elif sample_name == "TL1gDNAshort":
        return("Inoculum")
    else:
        raise ValueError('Incorrect sample_name: cannot extract gut site')

def extract_mouse_number(sample_name):
    mouse_number = sample_name[1:2]
    if mouse_number == "L":
        mouse_number = "Inoculum"
    
    return(mouse_number)

def extract_diet(sample_name):
    mouse_num = extract_mouse_number(sample_name)
    if mouse_num in ['1','2','3']:
        return("Control diet")
    elif mouse_num in ['4','5','6','7','8']:
        return("Guar gum diet")
    elif mouse_num == 'Inoculum':
        return("Human diet")
    else:
        raise ValueError('Incorrect sample_name: cannot extract diet')
        
def extract_cage(sample_name):
    mouse_num = extract_mouse_number(sample_name)
    
    if mouse_num in ['1','2','3']:
        return("Cage 1")
    elif mouse_num in ['4','5']:
        return("Cage 2")
    elif mouse_num in ['6','7', '8']:
        return("Cage 3")
    elif mouse_num == 'Inoculum':
        return("Inoculum")
    else:
        raise ValueError('Incorrect sample_name: cannot extract diet')
    
def abbreviate_gut_site(gut_site, blank_inoc = False):
    if (gut_site != "Cecum") & (gut_site != "Colon"):
        abbreviated_gut_site = gut_site[0:1]
    else:
        abbreviated_gut_site = gut_site[0:2]
    
    return(abbreviated_gut_site)

#####################################################################################
############################## STRAIN PLOTTING ######################################
#####################################################################################

        
def plot_strains(species):
    
    ### Directories ###
    strainfinder_dir = "%sinput" % (config.strain_phasing_directory)

    #Raw cluster
    raw_cluster_path = "%s%s" % (config.strain_phasing_directory, "strain_clusters/")
    species_raw_cluster_dir = "%s%s/" % (raw_cluster_path, species)
    if not os.path.exists(species_raw_cluster_dir):
        os.makedirs(species_raw_cluster_dir)
        print("Cluster directory created successfully!")
    else:
        print("Cluster directory already exists.")
    species_raw_cluster_path = "%s%s%s" % (species_raw_cluster_dir, species, "_RawCluster.pckl")

    #Centroid and polarized cluster paths
    species_centroid_cluster_path = "%s%s%s" % (species_raw_cluster_dir, species, "_ClusterCentroid.pckl")
    species_polarized_cluster_path = "%s%s%s" % (species_raw_cluster_dir, species, "_PolarizedCluster.pckl")

    #Evolutionary SNPs
    evo_snvs_directory = "%sstrain_phasing/snp_changes/%s/" % (config.project_directory, species)
    inoculum_to_mouse_changes_path = "%s%s%s" % (evo_snvs_directory, species,"_inoculum_mouse_changes.pckl")
    all_snp_changes_path = "%s%s%s" % (evo_snvs_directory, species,"_all_snp_changes.pckl")

    #MIDAS data
    annotated_snps_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species)

    ### Parameters ###
    min_cluster_size = 1000
    
    min_cluster_fraction = 1/10
    
    max_num_snvs = 20000
    
    ## minimum coverage to consider allele frequency at a site for purposes of clustering
    min_coverage = 20
    
    ## minimum average sample coverage at polymorphic sites (e.g. sites in the A/D matrices)
    min_sample_coverage = 5
    
    poly_cov_frac = 1/5
    
    n_clusters = 100
    
    #Minimum number of snvs per sample
    min_num_snvs_per_sample = 100
    
    ### Fs and other dataframes ###
    Fs,Ass,Dss = return_FAD(species, min_coverage=min_coverage, min_sample_coverage=min_sample_coverage, poly_cov_frac = poly_cov_frac, calculate_poly_cov_frac=True)
    
    ### filter out samples without an adequate number of SNVs when all is said and done ###
    sample_with_adequate_snv_count = ~((~np.isnan(Fs)).sum() < min_num_snvs_per_sample)

    Fs = Fs.loc[:,sample_with_adequate_snv_count]
    Ass = Ass.loc[:,sample_with_adequate_snv_count]
    Dss = Dss.loc[:,sample_with_adequate_snv_count]
    
    #### Copying the Inoculum ###
    Fs_inoculum = Fs.loc[:,(Fs.columns.get_level_values('sample') == 'TL1gDNAshort')]
    Ass_inoculum = Ass.loc[:,(Ass.columns.get_level_values('sample') == 'TL1gDNAshort')]
    Dss_inoculum = Dss.loc[:,(Dss.columns.get_level_values('sample') == 'TL1gDNAshort')]
    
    ### Loading polarized sites ###
    
    final_clusters = pd.read_pickle(species_polarized_cluster_path)
    final_f = pd.read_pickle(species_centroid_cluster_path)
    
    ### More ordering utilities ###
    mnum = list(set(Fs.T.index.get_level_values("mouse_number")))
    msite = list(set(Fs.T.index.get_level_values("region")))
    mdiet = list(set(Fs.T.index.get_level_values("diet")))
    mcage = list(set(Fs.T.index.get_level_values("cage")))

    mnum_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"mouse_number").index.get_level_values("mouse_number") == m).ravel() for m in list(set(mnum))}

    msite_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"region").index.get_level_values("region") == m).ravel() for m in list(set(msite))}

    mdiet_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"diet").index.get_level_values("diet") == m).ravel() for m in list(set(mdiet))}

    mcage_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"cage").index.get_level_values("cage") == m).ravel() for m in list(set(mcage))}

    all_sample_dics = {"diet":mdiet_sample_dic,"region":msite_sample_dic,"mouse_number":mnum_sample_dic,"cage":mcage_sample_dic}

    ######################################################################################################
    ############################################## PLOTTING ##############################################
    ######################################################################################################
    
def plot_strains(species, fig, ax):
    
    ### Directories ###
    strainfinder_dir = "%sinput" % (config.strain_phasing_directory)

    #Raw cluster
    raw_cluster_path = "%s%s" % (config.strain_phasing_directory, "strain_clusters/")
    species_raw_cluster_dir = "%s%s/" % (raw_cluster_path, species)
    if not os.path.exists(species_raw_cluster_dir):
        os.makedirs(species_raw_cluster_dir)
        print("Cluster directory created successfully!")
    else:
        print("Cluster directory already exists.")
    species_raw_cluster_path = "%s%s%s" % (species_raw_cluster_dir, species, "_RawCluster.pckl")

    #Centroid and polarized cluster paths
    species_centroid_cluster_path = "%s%s%s" % (species_raw_cluster_dir, species, "_ClusterCentroid.pckl")
    species_polarized_cluster_path = "%s%s%s" % (species_raw_cluster_dir, species, "_PolarizedCluster.pckl")
    
    #FINAL Fs path
    final_Fs_path = "%s%s%s" % (species_raw_cluster_dir, species, "_final_Fs.pckl")

    #Evolutionary SNPs
    evo_snvs_directory = "%sstrain_phasing/snp_changes/%s/" % (config.project_directory, species)
    inoculum_to_mouse_changes_path = "%s%s%s" % (evo_snvs_directory, species,"_inoculum_mouse_changes.pckl")
    all_snp_changes_path = "%s%s%s" % (evo_snvs_directory, species,"_all_snp_changes.pckl")

    #MIDAS data
    annotated_snps_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species)

    ### Parameters ###
    min_cluster_size = 1000
    
    min_cluster_fraction = 1/10
    
    max_num_snvs = 20000
    
    ## minimum coverage to consider allele frequency at a site for purposes of clustering
    min_coverage = 20
    
    ## minimum average sample coverage at polymorphic sites (e.g. sites in the A/D matrices)
    min_sample_coverage = 5
    
    poly_cov_frac = 1/5
    
    n_clusters = 100
    
    #Minimum number of snvs per sample
    min_num_snvs_per_sample = 100
    
    ### Loading polarized sites ###
    
    final_clusters = pd.read_pickle(species_polarized_cluster_path)
    final_f = pd.read_pickle(species_centroid_cluster_path)
    Fs = pd.read_pickle(final_Fs_path)
    
    ### More ordering utilities ###
    mnum = list(set(Fs.T.index.get_level_values("mouse_number")))
    msite = list(set(Fs.T.index.get_level_values("region")))
    mdiet = list(set(Fs.T.index.get_level_values("diet")))
    mcage = list(set(Fs.T.index.get_level_values("cage")))

    mnum_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"mouse_number").index.get_level_values("mouse_number") == m).ravel() for m in list(set(mnum))}

    msite_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"region").index.get_level_values("region") == m).ravel() for m in list(set(msite))}

    mdiet_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"diet").index.get_level_values("diet") == m).ravel() for m in list(set(mdiet))}

    mcage_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"cage").index.get_level_values("cage") == m).ravel() for m in list(set(mcage))}

    all_sample_dics = {"diet":mdiet_sample_dic,"region":msite_sample_dic,"mouse_number":mnum_sample_dic,"cage":mcage_sample_dic}

    ######################################################################################################
    ############################################## PLOTTING ##############################################
    ######################################################################################################
    
    key_to_sort = "mouse_number"
    secondary_key = "region"

    cmap_cage = []
    colors_library = [(1, 1, 0.8, 0.25), #creamy yellow
                      (0.529, 0.808, 0.922, 0.25), #sky blue
                      (0.529, 0.808, 0.922, 0.25)]
#                       (0.529, 0.808, 0.722, 0.25), #green blue
#                       (0.294, 0.427, 0.804, 0.15)] #purple blue
    if "Cage 1" in mcage_sample_dic:
    #     cmap_cage.append(colors_library(0))
        cmap_cage.append(colors_library[0])
    if "Cage 2" in mcage_sample_dic:
        cmap_cage.append(colors_library[1])

    if "Cage 3" in mcage_sample_dic:
        cmap_cage.append(colors_library[2])
    
    cmap_clus = get_cmap(5,name="Set3")

    #fig,ax = plt.subplots(figsize=(16,8))
    ax.set_title(fu.get_pretty_species_name(species),size=30)

    for i,f in enumerate(final_f):

        ff = reorder_sort(f,key_to_sort).sort_index(key=lambda x: x.map(order_dict))
        ax.plot(ff.values,zorder=100,lw=6,color=cmap_clus(i),label=f"Strain {i+1}");
        ax.plot(ff.values,zorder=80,lw=7,color="k");
        if len(final_clusters) != 0:
            ff_c = reorder_sort(final_clusters[i].T,key_to_sort).T.sort_index(key=lambda x: x.map(order_dict),axis=1)
            ax.plot(ff_c.sample(min(ff_c.shape[0],10000)).T.values,color=cmap_clus(i),alpha=.01)
        else:
            for i in np.arange(len([Fs])):
                ff_c = reorder_sort((1-Fs).T,key_to_sort).T.sort_index(key=lambda x: x.map(order_dict),axis=1)
                ax.plot(ff_c.sample(min(ff_c.shape[0],10000)).T.values,color=cmap_clus(i),alpha=.01)

    major_x = []
    minor_x = []
    labels = []

    second_xlabels = list(ff.index.get_level_values(secondary_key))


    #Making the vertical lines and labels
    i = 0
    for key, item in all_sample_dics[key_to_sort].items():

        xmin = item.min() 
        xmax = item.max()

        for e in item:
            if key == "Inoculum":
                ax.axvline(e,color="k",zorder=0,alpha=1)
            else:
                ax.axvline(e,color="k",zorder=0,alpha=.5)   

            ax.text(e, -0.1, abbreviate_gut_site(second_xlabels[e], blank_inoc = True), ha='center',va='top', clip_on=False,size=15, rotation=0)

        if xmin != xmax:
            major_x.extend([xmin,(xmax + xmin)/2,xmax])
            minor_x.append((xmax + xmin)/2)
            labels.extend(["",key,""])
        else:
            major_x.append(xmin)
            minor_x.append(xmax)
            labels.extend([key])   

        i+=1

        if ("Inoculum" not in all_sample_dics[key_to_sort]) & (max(all_sample_dics[key_to_sort]) == key):
            ax.vlines(item[0] - 0.5, 0, -0.15, color='black', lw=0.8, clip_on=False, transform=ax.get_xaxis_transform())
            ax.vlines(item[-1] + 0.5, 0, -0.15, color='black', lw=2, clip_on=False, transform=ax.get_xaxis_transform())

        elif key == "Inoculum":
            ax.vlines(item[0] - 0.5, 0, -0.15, color='black', lw=0.8, clip_on=False, transform=ax.get_xaxis_transform())
        else:
            ax.vlines(item[0] - 0.5, 0, -0.15, color='black', lw=0.8, clip_on=False, transform=ax.get_xaxis_transform())
            ax.vlines(item[-1] + 0.5, 0, -0.15, color='black', lw=0.8, clip_on=False, transform=ax.get_xaxis_transform())

    #Making the vertical colors
    key_to_sort = "cage"
    i = 0    
    for key, item in sorted(all_sample_dics[key_to_sort].items()):

        xmin = item.min() 
        xmax = item.max()

        if (key == max(all_sample_dics[key_to_sort])):
            ax.vlines(item[0] - 0.5, 0, -0.15, color='black', lw=2, clip_on=False, transform=ax.get_xaxis_transform())
            if "Inoculum" != max(all_sample_dics[key_to_sort]):
                ax.vlines(item[-1] + 0.5, 0, -0.15, color='black', lw=2, clip_on=False, transform=ax.get_xaxis_transform())
                ax.axvspan(xmin - .25,xmax+.25,color=cmap_cage[i]) #,alpha=.2
            else:
                ax.axvspan(xmin - .25,xmax+.25,alpha=.2,color='black') 

            i+=1
        else:
            ax.axvspan(xmin - .1,xmax+.1,color=cmap_cage[i]) #alpha=.2,
            ax.vlines(item[0] - 0.5, 0, -0.15, color='black', lw=2, clip_on=False, transform=ax.get_xaxis_transform())
            ax.vlines(item[-1] + 0.5, 0, -0.15, color='black', lw=2, clip_on=False, transform=ax.get_xaxis_transform())

            i+=1


    ax.set_xticks(major_x)
    ax.set_xticks(minor_x, minor = True)
    ax.set_xticklabels([label if label != "Inoculum" else "" for label in labels]);

    ax.axhline(0,color="grey")
    ax.axhline(1,color="grey")

#     ax.tick_params(axis = 'x', which = 'major', length=0,labelsize = 20,pad=45)
    ax.tick_params(axis = 'x', which = 'major', length=0,labelsize = 20, pad = 30)
    ax.tick_params(axis = 'x', which = 'minor', length = 10,labelsize = 0)

#     ax.set_ylabel("Strain frequency",size=30)
    ax.set_ylim([-0.05,1.05]);
    if "Inoculum" in labels:
        ax.set_xlim([-1,max(major_x)+0.25]);


#     ax.set_xlabel("Sample", size = 30)


#     fig.legend(prop={"size":20});

    ############ ADDING SFS #################
    # Create a new axis on the right side
    if ("Inoculum" in labels):

        ax_hist = ax.inset_axes([1, 0, 0.1, 1])

        if len(final_clusters) ==0:
            hist_data = Fs["Inoculum"].dropna().values
            ax_hist.hist(hist_data, alpha = 0.5, bins = 25, color = cmap_clus(0), orientation = "horizontal")
        else:   
            for i,final_cluster in enumerate(final_clusters):
                hist_data = final_cluster["Inoculum"].dropna().values
                ax_hist.hist(hist_data, alpha = 0.5, bins = 25, color = cmap_clus(i), orientation = "horizontal")



        ax_hist.invert_yaxis()
        ax_hist.set_yticks([])
        ax_hist.set_xticks([])
        ax_hist.set_frame_on(False)
        ax_hist.invert_yaxis()

        # Set the title for the histogram
        ax_hist.set_xlabel("Inoculum", fontsize = 20, labelpad=17.5)
    

    
