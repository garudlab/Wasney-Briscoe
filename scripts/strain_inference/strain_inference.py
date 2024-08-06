import sys
### ADD PATH HERE
sys.path.insert(0, "~/Wasney-Briscoe-2024/scripts/postprocessing/postprocessing_scripts/")
import config

### Packages
import pandas as pd
import numpy as np
import scipy.stats

import os 
from scipy.spatial.distance import pdist,squareform
import figure_utils as fu
from numba import njit 
import pickle
import random as rand
from random import randint,sample,choices
from math import log

### Plotting packages
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}') 
import seaborn as sns
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import matplotlib.cm as cm
plasma_cmap = cm.get_cmap('plasma')
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from strain_phasing_functions import *


#MIDAS data
annotated_snps_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species)

## Meta-parameters: experiment with these—no hard and fast rules!

## minimum number of SNVs which need to be clustered together in order to qualify as a "strain"
## if we didn't cap max_num_snvs, then min_cluster_size would be O(10^4), based on the typical evolutionary
## distance between strains
min_cluster_size = 1000
# min_cluster_size = 50

## minimum fraction of sites which pass our coverage threshold which must be in a cluster in order for it to qualify 
## as a strain
## basically, the idea is that as the initial number of sites we pass in gets bigger, we want to incrwease the min_cluster_size
## here, we say that 10% of all variable sites must be in a cluster in order for it to be considered a "strain"
## this will largely be redundant w/ min_cluster_size, but adds some more functionality to play with
min_cluster_fraction = 1/10
#min_cluster_fraction = 0

## For computational efficiency, we can downsample the SNVs we actually perform strain phasing on
## should still give us the same strain trajectory 
## clustering 20k SNVs takes ~90 seconds. 
max_num_snvs = 20000

## distance threshold to be considered linked—lower means trajectories have to be more 
## similar, higher means less similar, to be in a cluster
    
max_d = 3.5
# max_d = 6 #USE THIS FOR B. WEXLERAE
#max_d = 4 #USE THIS FOR B. UNIFORMIS, A. hadrus
#max_d = 5 #Coprococcus_sp_62244


## minimum coverage to consider allele frequency at a site for purposes of clustering
min_coverage = 10 #MW: 20 until 6/30/2024
# min_coverage = 5

## minimum average sample coverage at polymorphic sites (e.g. sites in the A/D matrices)
min_sample_coverage = 5
# min_sample_coverage = 0


## polymorphic & covered fraction: what percentage of samples does a site need 
## with coverage > min_coverage and polymorphic to be included in downstream analyses? 
## NOTE: we may want to disaggregate coverage and polymorphic-ness so as to not lose evolutionary snvs
## but for strain clustering purposes, I think we should focus on SNVs that are actually polymorphic
## in a good number of samples
poly_cov_frac = 1/5 #

## Number of clusters to calculate
n_clusters = 100
#n_clusters = 500

#Minimum number of snvs per sample
min_num_snvs_per_sample = 100

#Loading the frequency, alternate read, and total depth dataframes
Fs,Ass,Dss = return_FAD(species, min_coverage=min_coverage, 
                        min_sample_coverage=min_sample_coverage, 
                        poly_cov_frac = poly_cov_frac, 
                        calculate_poly_cov_frac=True, 
                        read_support = True) 

#Copying the Inoculum
Fs_inoculum = Fs.loc[:,(Fs.columns.get_level_values('sample') == 'TL1gDNAshort')]
Ass_inoculum = Ass.loc[:,(Ass.columns.get_level_values('sample') == 'TL1gDNAshort')]
Dss_inoculum = Dss.loc[:,(Dss.columns.get_level_values('sample') == 'TL1gDNAshort')]



%%time

fss = Ass.values/(Dss.values + (Dss.values == 0)) #This is so it doesn't produce a na (division by 0)

cluster_As = Ass.values
cluster_Ds = Dss.values
cluster_fs = cluster_As/(cluster_Ds + (cluster_Ds == 0))

## for compatibility in case of threshold number of SNVs
num = min(max_num_snvs,Fs.shape[0])

i_list = Dss.T.mean().sort_values(ascending=False).index[:num]

sys.stderr.write("Processing %s SNVs" % num)

## simply shuffles indices if no threshold is specified
#i_list = sample(range(Fs.shape[0]),num)
i_list_idx = Fs.loc[i_list].index

Ass_sub = Ass.loc[i_list_idx]
Dss_sub = Dss.loc[i_list_idx]
Fs_sub = Fs.loc[i_list_idx]

fss_sub = Ass_sub.values/(Dss_sub.values + (Dss_sub.values == 0))

cluster_As_sub = Ass_sub.values
cluster_Ds_sub = Dss_sub.values
cluster_fs_sub = cluster_As_sub/(cluster_Ds_sub + (cluster_Ds_sub == 0))

D_mat = np.zeros([num,num])
D_mat_1 = D_mat_fun1(num,fss_sub,cluster_Ds_sub,D_mat)
D_mat = np.zeros([num,num]) 
D_mat_2 = D_mat_fun2(num,fss_sub,cluster_Ds_sub,D_mat)

D_mat = np.fmin(D_mat_1,D_mat_2) #I believe this is filling in the minimum of the two polarizations
D_mat = symmetrize(D_mat)

D_mat_1 = pd.DataFrame(D_mat_1,index=Fs_sub.index,columns=Fs_sub.index)
D_mat_2 = pd.DataFrame(D_mat_2,index=Fs_sub.index,columns=Fs_sub.index)

D_mat_close = pd.DataFrame(D_mat < max_d) 

D_mat_close.index = Fs_sub.index
D_mat_close.columns = Fs_sub.index


## extracts up to 100 clusters
## in practice all SNVs should fall into one of a fairly small number of clusters
## really should re-write this with a while loop but this works for now
## the idea is that we exhaust all clusters—there should only be a small number of them ultimately

###Idea with while loop:
##### While there are still variants out there, have it try to be clusterings

all_clus_pol = []
all_clus_idx = []
all_clus_A = []
all_clus_D = []

all_clus_F = []

if species == "Bacteroides_uniformis_57318": #MW 7/17/2024
    
    clus_1,clus_idxs_1 = return_clus(D_mat_close,Fs_sub)
    clus_pol_1 = polarize_clus(clus_1,clus_idxs_1,D_mat_1,D_mat_2)
    clus_pol_1.index = clus_idxs_1
    D_mat_close = drop_clus_idxs(D_mat_close,clus_idxs_1)
    
    clus_2,clus_idxs_2 = return_clus(D_mat_close,Fs_sub)
    clus_pol_2 = polarize_clus(clus_2,clus_idxs_2,D_mat_1,D_mat_2)
    clus_pol_2.index = clus_idxs_2
    D_mat_close = drop_clus_idxs(D_mat_close,clus_idxs_2)
    
    clus = pd.concat([clus_1, clus_2], axis = 0)
    clus_idxs = clus_idxs_1.append(clus_idxs_2)
    clus_pol = pd.concat([clus_pol_1, clus_pol_2], axis = 0)
    
    all_clus_D.append(Dss.loc[clus.index].mean().values)
    all_clus_pol.append(clus_pol)
    all_clus_A.append(clus_pol.mean()*all_clus_D[-1])
    all_clus_F.append(clus_pol.mean())
    
    print(clus_pol.shape[0])

    
else:

    for i in range(n_clusters):

        try:

            clus,clus_idxs = return_clus(D_mat_close,Fs_sub)
    #         clus,clus_idxs = return_clus(D_mat_close,Fs_sub, co_cluster_pct=0.5) #Finding points that cluster with 25% other points. That's a cluster.
                                                             #We would modify this function to get smaller clusters...
            clus_pol = polarize_clus(clus,clus_idxs,D_mat_1,D_mat_2)
            clus_pol.index = clus_idxs
            D_mat_close = drop_clus_idxs(D_mat_close,clus_idxs)

            if clus_pol.shape[0] > min_cluster_size and clus_pol.shape[0] > Fs.shape[0]*min_cluster_fraction:

                all_clus_D.append(Dss.loc[clus.index].mean().values)
                all_clus_pol.append(clus_pol)
                all_clus_A.append(clus_pol.mean()*all_clus_D[-1])
                all_clus_F.append(clus_pol.mean())

                print(clus_pol.shape[0])

        except:
            pass

## now, choosing a representative SNV from each cluster, and finding all other sites (not just limited to the 20k)
## which are consistent w/ being linked to it

final_clusters = []

all_aligned_sites = []

for i in range(len(all_clus_D)):
    
    sys.stderr.write(f'\n\nCluster {i+1}\n')
    ancD = all_clus_D[i]
    ancF = all_clus_F[i]

    dss = Dss.values
    fss = Fs.values
    
    disAnc_forward = []
    disAnc_backward = []

    for j in range(Dss.shape[0]):
        disAnc_forward.append(calc_dis(ancD,dss[j],ancF,fss[j]))
        disAnc_backward.append(calc_dis(ancD,dss[j],ancF,1-fss[j]))
        if j % 1000 == 0:
            sys.stderr.write(f"\n\t{np.around(100*j/Dss.shape[0],3)}% finished")
    
    disAnc = [min(els) for els in zip(disAnc_forward, disAnc_backward)]
    disAnc = np.array(disAnc)
    aligned_sites = Fs.loc[disAnc < max_d].index
    f_dist =  pd.DataFrame(np.array([disAnc_forward,disAnc_backward]).T,index=Fs.index)
    pols = f_dist.T.idxmin() > 0
    
    aligned_sites = [a for a in aligned_sites if a not in all_aligned_sites]
    
    pols = pols.loc[aligned_sites]
    re_polarize = pols.loc[pols].index
    
    all_aligned_sites.extend(aligned_sites)
    
    Fs_cluster = Fs.loc[aligned_sites]
    
    Fs_cluster.loc[re_polarize] = 1 - Fs_cluster.loc[re_polarize]
        
    final_clusters.append(Fs_cluster)

#Creating polarized clusters
no_cluster = False
multiple_inoculum_strains = False
if (len(final_clusters) == 0):
    sys.stderr.write("No clusters detected.\n")
    no_cluster = True
    
    Fs_nc,Ass_nc,Dss_nc = return_FAD(species, 
                            min_coverage=min_coverage,
                            min_sample_coverage=min_sample_coverage,
                            poly_cov_frac = 0, 
                            calculate_poly_cov_frac=False, 
                            read_support = False)
    
    if ("Inoculum" not in Fs_nc):
        sys.stderr.write("Inoculum not in data.\n")
        if Fs_nc.mean().mean() < 0.5:
            df_final_f = 1 - Fs_nc.mean().T
            df_final_f.loc[:,:] = 1
            final_f = []
            final_f.append(df_final_f)
            final_clusters = []
            final_clusters.append(1 - Fs_nc)
        else:
            df_final_f = Fs_nc.mean().T
            df_final_f.loc[:,:] = 1
            final_f = []
            final_f.append(df_final_f)
            final_clusters = []
            final_clusters.append(Fs_nc)
        
    elif ((((Fs_nc > 0.0) | (Fs_nc < 1.0)) & (Ass_nc > 4)).sum()["Inoculum"].values[0] > min_cluster_size):
        #If the number of supported polymoprhic SNVs in the inoculum > 1000, infer two strains in the inoculum
        sys.stderr.write("Multiple strains detected in inoculum.\n")
        multiple_inoculum_strains = True
        Fs_nc = Fs_nc.loc[:,((Dss_nc > min_coverage).sum() > min_cluster_size)]
        Ass_nc = Ass_nc.loc[:,((Dss_nc > min_coverage).sum() > min_cluster_size)]
        Dss_nc = Dss_nc.loc[:,((Dss_nc > min_coverage).sum() > min_cluster_size)] #Number of good coverage sites exceeds 1000
        
        if Fs_nc.mean().mean() < 0.5:
            #if alleles are polarized such that the alt is < 0.5, flip the polarization
            #extract mean using only polymorphic sites in the inoculum
            Fs_nc = Fs_nc.loc[Fs_nc['Inoculum'][(Fs_nc['Inoculum'] > 0) & (Dss_nc['Inoculum'] > 4)].dropna().index]
            Ass_nc = Ass_nc.loc[Fs_nc['Inoculum'][(Fs_nc['Inoculum'] > 0) & (Dss_nc['Inoculum'] > 4)].dropna().index]
            Dss_nc = Dss_nc.loc[Fs_nc['Inoculum'][(Fs_nc['Inoculum'] > 0) & (Dss_nc['Inoculum'] > 4)].dropna().index]
            df_final_f = 1 - Fs_nc.mean().T
            columns_to_update = pd.DataFrame(df_final_f).T.columns.difference([('Inoculum', 'Inoculum', 'Human diet', 'Inoculum', 'TL1gDNAshort')])
            df_final_f = pd.DataFrame(df_final_f).T
            df_final_f.loc[:,columns_to_update] = 1
            df_final_f = df_final_f.mean()
            final_f = []
            final_f.append(df_final_f)
            final_f.append(1 - df_final_f)
            final_clusters = []
            final_clusters.append(1 - Fs_nc)
            final_clusters.append(Fs_nc)
            
        else:
            Fs_nc = Fs_nc.loc[Fs_nc['Inoculum'][(Fs_nc['Inoculum'] < 1) & (Dss_nc['Inoculum'] > 4)].dropna().index]
            Ass_nc = Ass_nc.loc[Fs_nc['Inoculum'][(Fs_nc['Inoculum'] < 1) & (Dss_nc['Inoculum'] > 4)].dropna().index]
            Dss_nc = Dss_nc.loc[Fs_nc['Inoculum'][(Fs_nc['Inoculum'] < 1) & (Dss_nc['Inoculum'] > 4)].dropna().index]
            df_final_f = Fs_nc.mean().T
            columns_to_update = pd.DataFrame(df_final_f).T.columns.difference([('Inoculum', 'Inoculum', 'Human diet', 'Inoculum', 'TL1gDNAshort')])
            df_final_f = pd.DataFrame(df_final_f).T
            df_final_f.loc[:,columns_to_update] = 1
            df_final_f = df_final_f.mean()
            final_f = []
            final_f.append(df_final_f)
            final_f.append(1-df_final_f)
            final_clusters = []
            final_clusters.append(Fs_nc)
            final_clusters.append(1 - Fs_nc)
            
    else:
        sys.stderr.write("A single strain detected in inoculum.\n")
        Fs_nc = Fs_nc.loc[:,((Dss_nc > min_coverage).sum() > min_cluster_size)]
        Ass_nc = Ass_nc.loc[:,((Dss_nc > min_coverage).sum() > min_cluster_size)]
        Dss_nc = Dss_nc.loc[:,((Dss_nc > min_coverage).sum() > min_cluster_size)] #Number of good coverage sites exceeds 1000
        df_final_f = Fs_nc.mean().T
        df_final_f.loc[:,:] = 1
        final_f = []
        final_f.append(df_final_f)
        if (Fs_nc.mean().mean() < 0.5): 
            #if alleles are polarized such that the alt is < 0.5, flip the polarization
            final_clusters = []
            final_clusters.append(1 - Fs_nc)
        else:
            final_clusters = []
            final_clusters.append(Fs_nc)
else:
    sys.stderr.write("Multiple strains detected.\n")
    ## If only a single cluster is detected, add a second "cluster" which is simply 1 minus the allele frequencies
    ## in the first cluster
    ## aids in visualization for people not familiar with this kind of clustering
    if len(final_clusters) == 1:
        final_clusters.append(1-final_clusters[0])
    
    ## add cluster centroids
    final_f = []
    for cluster in final_clusters:
        final_f.append(cluster.mean())
    df_final_f = pd.DataFrame(final_f)

    ## now, polarize clusters so that the sum of squareds of the centroids to 1 is minimized
    ## the idea here is that accurate strain frequencies should sum to 1
    polarize = True
    
    pol_d2 = {}

    for i in range(df_final_f.shape[0]):
        df_final_f_temp = df_final_f.copy() #Makes a copy of the centroids
        df_final_f_temp.iloc[i] = 1 - df_final_f_temp.iloc[i] #gets the polarized version of ONE of the centroids.
        pol_d2[i] =  ((1 - df_final_f_temp.sum())**2).sum()   #Get the across centroids for all samples (should be close to 1), 
                                                                  #subtract this from 1, and square. Sum all those values
                                                                  #Ideally, this value is really close to 0. 
                                                                  #Add this value to the dictionary.

        pol_d2 = pd.Series(pol_d2)                                #Make the dictionary a series 

        if pol_d2.min() < ((1 - df_final_f.sum())**2).sum(): #If any of the above repolarizations actually made the overall sum of centroids closer to 1, repolarize.
            clus_to_re_pol = pol_d2.idxmin()
            final_f[clus_to_re_pol] = 1 - final_f[clus_to_re_pol]
            final_clusters[clus_to_re_pol] = 1 - final_clusters[clus_to_re_pol]
            df_final_f = pd.DataFrame(final_f)  
            Fs_inoculum.loc[final_clusters[clus_to_re_pol].index] = 1 -  Fs_inoculum.loc[final_clusters[clus_to_re_pol].index]#polarize inoculum accordingly

    

#Filtering out samples without adequate coverage
good_indices = []

for i,cluster in enumerate(final_clusters):
    if i == 0:
        good_samples = (len(final_clusters[i]) - np.isnan(final_clusters[i]).sum(axis = 0) > min_num_snvs_per_sample).values
    else:
        new_good_samples = (len(final_clusters[i]) - np.isnan(final_clusters[i]).sum(axis = 0) > min_num_snvs_per_sample).values
        good_samples = good_samples & new_good_samples

for i,cluster in enumerate(final_clusters): 
    final_clusters[i] = final_clusters[i].T.loc[good_samples].T
    final_f[i] = final_f[i].T.loc[good_samples]
if no_cluster:
    Fs = Fs_nc.T.loc[good_samples].T
    Ass = Ass_nc.T.loc[good_samples].T
    Dss = Dss_nc.T.loc[good_samples].T
else:
    Fs = Fs.T.loc[good_samples].T
    Ass = Ass.T.loc[good_samples].T
    Dss = Dss.T.loc[good_samples].T

#Filter all all na columns if there are any - THis is redundant
if len(final_clusters) > 0:
    for i,cluster in enumerate(final_clusters):
        if i == 0:
            mask = ~np.isnan(cluster).all(axis = 0)
        final_clusters[i] = cluster.loc[:,mask]
        final_f[i] = final_f[i][mask]
    Fs = Fs.loc[:,mask]

#Creating output for R plotting (with bootstrapped 95% confidence intervals
bootstrap_ci = True
boostrap_k = 100
bootstrap_N = 1000
for i in np.arange(len(final_f)): 
    if i == 0:
        final_f_all_strains = pd.DataFrame(final_f[i]).reset_index().rename(columns = {0: "freq"})
        final_f_all_strains['species'] = [species]*len(final_f_all_strains)
        final_f_all_strains['strain'] = [i + 1]*len(final_f_all_strains)
        final_f_all_strains = final_f_all_strains[['species','strain','mouse_number', 'region', 'diet', 'cage', 'sample', 'freq']]
        iqr_75 = pd.DataFrame(final_clusters[i].quantile(0.75)).reset_index()
        iqr_25 = pd.DataFrame(final_clusters[i].quantile(0.25)).reset_index()
        final_f_all_strains = final_f_all_strains.merge(iqr_25).merge(iqr_75)
        if bootstrap_ci:
            upper_ci_vec = []
            lower_ci_vec = []
            
            for sample_i in range(final_clusters[i].shape[1]):
                
                sample_name = final_clusters[i].iloc[:,sample_i].name[4]
                if  (no_cluster) & (multiple_inoculum_strains) & (sample_name != "TL1gDNAshort"):
                    upper_ci_vec.append(np.nan)
                    lower_ci_vec.append(np.nan)
                
                elif (no_cluster) & (not multiple_inoculum_strains):
                    upper_ci_vec.append(np.nan)
                    lower_ci_vec.append(np.nan)
                    
                else:
                    snv_freqs = final_clusters[i].iloc[:,sample_i].dropna().to_list()
                    mean_freq_vec = []
                    for n in range(bootstrap_N):
                        sampled_snv_freqs = choices(snv_freqs, k = boostrap_k)
                        mean_freq = np.mean(sampled_snv_freqs)
                        mean_freq_vec.append(mean_freq)
                    upper_ci = np.quantile(mean_freq_vec, 0.975)
                    lower_ci = np.quantile(mean_freq_vec, 0.025)
                    upper_ci_vec.append(upper_ci)
                    lower_ci_vec.append(lower_ci)
                
            

    else:

        final_f_all_strains_temp = pd.DataFrame(final_f[i]).reset_index().rename(columns = {0: "freq"})
        final_f_all_strains_temp['species'] = [species]*len(final_f_all_strains_temp)
        final_f_all_strains_temp['strain'] = [i + 1]*len(final_f_all_strains_temp)
        final_f_all_strains_temp = final_f_all_strains_temp[['species','strain','mouse_number', 'region', 'diet', 'cage', 'sample', 'freq']]
        iqr_75 = pd.DataFrame(final_clusters[i].quantile(0.75)).reset_index()
        iqr_25 = pd.DataFrame(final_clusters[i].quantile(0.25)).reset_index()
        final_f_all_strains_temp = final_f_all_strains_temp.merge(iqr_25).merge(iqr_75)
        final_f_all_strains = final_f_all_strains.append(final_f_all_strains_temp, ignore_index = True)
        if bootstrap_ci:
            
            for sample_i in range(final_clusters[i].shape[1]): #if we did not infer multiple strains, but there may have been multiple strains in the inoculum, just calculate CI for inoculum
                sample_name = final_clusters[i].iloc[:,sample_i].name[4]
                if  (no_cluster) & (multiple_inoculum_strains) & (sample_name != "TL1gDNAshort"):
                    upper_ci_vec.append(np.nan)
                    lower_ci_vec.append(np.nan)
               
                elif (no_cluster) & (not multiple_inoculum_strains):
                    upper_ci_vec.append(np.nan)
                    lower_ci_vec.append(np.nan)
                    
                else:
                    snv_freqs = final_clusters[i].iloc[:,sample_i].dropna().to_list()
                    mean_freq_vec = []
                    for n in range(bootstrap_N):
                        sampled_snv_freqs = choices(snv_freqs, k = boostrap_k)
                        mean_freq = np.mean(sampled_snv_freqs)
                        mean_freq_vec.append(mean_freq)
                    upper_ci = np.quantile(mean_freq_vec, 0.975)
                    lower_ci = np.quantile(mean_freq_vec, 0.025)
                    upper_ci_vec.append(upper_ci)
                    lower_ci_vec.append(lower_ci)
        
        
    final_f_all_strains['upper_ci'] = upper_ci_vec
    final_f_all_strains['lower_ci'] = lower_ci_vec
#Renaming
final_f_all_strains = final_f_all_strains.rename(columns = {0.25: "quantile_25", 0.75: "quantile_75"})

#Saving 
### Make strain_cluster directory if necessary
output_dir = "%sstrain_phasing/strain_clusters/%s/" % (config.project_folder, species)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
### Saving
output_file = "%sstrain_phasing/strain_clusters/%s/%s_strain_frequency.csv" % (config.project_folder, species, species)
final_f_all_strains.to_csv(output_file, sep = "\t", index = False)


