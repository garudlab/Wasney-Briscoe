### SET UP ####

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob

import bz2

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}') 
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import scipy.stats
import figure_utils as fu
from return_gene_descriptions import return_gene_descriptions

from numba import njit 

import matplotlib.pyplot as plt
import pickle
import pandas as pd
import config
import numpy
import random as rand

from random import randint,sample
from math import log

import matplotlib.cm as cm
plasma_cmap = cm.get_cmap('plasma')
import matplotlib.colors as mcolors

import sys
import os 
from scipy.spatial.distance import pdist,squareform

import matplotlib.gridspec as gridspec

from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
from matplotlib.colors import ListedColormap

hap_cmap = ListedColormap(['grey', 'red', 'black', 'black','blue'], 'indexed')

from strain_phasing_functions import *

#### SCRIPT ####
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--species", type=str, help="Loads the species to be processed.", default="Bacteroides_vulgatus_57955")
args = parser.parse_args()

species=args.species

# DIRECTORIES
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

## Meta-parameters: experiment with these—no hard and fast rules!

## minimum number of SNVs which need to be clustered together in order to qualify as a "strain"
## if we didn't cap max_num_snvs, then min_cluster_size would be O(10^4), based on the typical evolutionary
## distance between strains
min_cluster_size = 1000
#min_cluster_size = 100

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

if species == "Blautia_wexlerae_56130":
    max_d = 3.75
elif species == "Bacteroides_uniformis_57318":
    max_d = 4
elif species == "Eubacterium_hallii_61477":
    max_d = 4.25
else:
    max_d = 3.5


## minimum coverage to consider allele frequency at a site for purposes of clustering
min_coverage = 20

## minimum average sample coverage at polymorphic sites (e.g. sites in the A/D matrices)
min_sample_coverage = 5

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

Fs,Ass,Dss = return_FAD(species, min_coverage=min_coverage, min_sample_coverage=min_sample_coverage, poly_cov_frac = poly_cov_frac, calculate_poly_cov_frac=True)

#filter out samples without an adequate number of SNVs when all is said and done
sample_with_adequate_snv_count = ~((~np.isnan(Fs)).sum() < min_num_snvs_per_sample)

Fs = Fs.loc[:,sample_with_adequate_snv_count]
Ass = Ass.loc[:,sample_with_adequate_snv_count]
Dss = Dss.loc[:,sample_with_adequate_snv_count]


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


for i in range(n_clusters):
    
    try:
        
        clus,clus_idxs = return_clus(D_mat_close,Fs_sub) #Finding points that cluster with 25% other points. That's a cluster.
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

if "Inoculum" in Fs:

    inoculum_AF = Fs.Inoculum.dropna().values.ravel() 
    percent_polymorphic_inoculum = ((inoculum_AF >= 0.2) & (inoculum_AF <= 0.8)).sum()/len(inoculum_AF)

if len(final_clusters) == 0:
    df_final_f = pd.DataFrame(Fs.mean()).T
    df_final_f.loc[:,:] = 1
    final_f = []
    final_f.append(df_final_f.mean())

elif (((final_clusters[0].mean() < 0.05).sum()/len(final_clusters[0].columns) > 0.75) | ((final_clusters[0].mean() > 0.95).sum()/len(final_clusters[0].columns) > 0.75)):   #Need to edit this so it's selecting species with multiple strains in the inoculum 
# elif (((final_clusters[0].mean() < 0.05).sum() == len(final_clusters[0].columns) - 1) | ((final_clusters[0].mean() > 0.95).sum() == len(final_clusters[0].columns) - 1)):
    if len(final_clusters) == 1:
        final_clusters.append(1-final_clusters[0])
        final_f = []
        for cluster in final_clusters:
            final_f.append(cluster.mean())
    
elif (len(final_clusters) == 1) & ((final_clusters[0].mean() < 0.05).sum()/len(final_clusters[0].columns) > 0.75):
    df_final_f = pd.DataFrame(Fs.mean()).T
    df_final_f.loc[:,:] = 1
    final_f = []
    final_f.append(df_final_f.mean())
    final_clusters[0] = 1 - final_clusters[0]

else:
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

    if polarize:

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

            
#Filter all all na columns if there are any
if len(final_clusters) > 0:
    for i,cluster in enumerate(final_clusters):
        if i == 0:
            mask = ~np.isnan(cluster).all(axis = 0)
        final_clusters[i] = cluster.loc[:,mask]
        final_f[i] = final_f[i][mask]
    Fs = Fs.loc[:,mask]

## More ordering utilities
mnum = list(set(Fs.T.index.get_level_values("mouse_number")))
msite = list(set(Fs.T.index.get_level_values("region")))
mdiet = list(set(Fs.T.index.get_level_values("diet")))
mcage = list(set(Fs.T.index.get_level_values("cage")))

mnum_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"mouse_number").index.get_level_values("mouse_number") == m).ravel() for m in list(set(mnum))}

msite_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"region").index.get_level_values("region") == m).ravel() for m in list(set(msite))}

mdiet_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"diet").index.get_level_values("diet") == m).ravel() for m in list(set(mdiet))}

mcage_sample_dic = {m:np.argwhere(reorder_sort(Fs.T,"cage").index.get_level_values("cage") == m).ravel() for m in list(set(mcage))}

all_sample_dics = {"diet":mdiet_sample_dic,"region":msite_sample_dic,"mouse_number":mnum_sample_dic,"cage":mcage_sample_dic}

#### PLOTTING #######
key_to_sort = "mouse_number"
secondary_key = "region"

cmap = get_cmap(len(list(set(mnum))))
# cmap_clus = get_cmap(len(list(set(mnum))),name="Set3")
cmap_clus = get_cmap(5,name="Set3")

fig,ax = plt.subplots(figsize=(16,8))
fig.suptitle(fu.get_pretty_species_name(species),size=30)

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
for key, item in all_sample_dics[key_to_sort].items():
    xmin = item.min() 
    xmax = item.max()
        
    if (key == max(all_sample_dics[key_to_sort])):
        ax.vlines(item[0] - 0.5, 0, -0.15, color='black', lw=2, clip_on=False, transform=ax.get_xaxis_transform())
        if "Inoculum" != max(all_sample_dics[key_to_sort]):
            ax.vlines(item[-1] + 0.5, 0, -0.15, color='black', lw=2, clip_on=False, transform=ax.get_xaxis_transform())
            ax.axvspan(xmin - .25,xmax+.25,alpha=.2,color=cmap(i))
        else:
            ax.axvspan(xmin - .25,xmax+.25,alpha=.2,color='black')

        i+=1
    else:
        ax.axvspan(xmin - .1,xmax+.1,alpha=.2,color=cmap(i))
        ax.vlines(item[0] - 0.5, 0, -0.15, color='black', lw=2, clip_on=False, transform=ax.get_xaxis_transform())
        ax.vlines(item[-1] + 0.5, 0, -0.15, color='black', lw=2, clip_on=False, transform=ax.get_xaxis_transform())

        i+=1


ax.set_xticks(major_x)
ax.set_xticks(minor_x, minor = True)
ax.set_xticklabels([label if label != "Inoculum" else "" for label in labels]);

ax.axhline(0,color="grey")
ax.axhline(1,color="grey")

ax.tick_params(axis = 'x', which = 'major', length=0,labelsize = 20,pad=45)
ax.tick_params(axis = 'x', which = 'minor', length = 10,labelsize = 0)
    
ax.set_ylabel("Strain frequency",size=30)
ax.set_ylim([-0.05,1.05]);
if "Inoculum" in labels:
    ax.set_xlim([-1,max(major_x)+0.25]);
    

ax.set_xlabel("Sample", size = 30)


fig.legend(prop={"size":20});
        
############ ADDING SFS #################
# Create a new axis on the right side
if ("Inoculum" in labels):

    ax_hist = ax.inset_axes([1, 0, 0.1, 1])
    
    if len(final_clusters) == 0:
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
#plt.tight_layout()

######## SAVING
if (len(final_f) > 0) & (len(Fs.columns) > 0):

    #Figure directory
    figure_directory = "%s%s%s" % (config.figure_directory, "strain_phasing/", species)
    if not os.path.exists(figure_directory):
        os.makedirs(figure_directory)
        print("Figure directory created successfully!")
    else:
        print("Figure directory already exists.")

    figure_path = "%s/%s_snp_changes.png" % (figure_directory, species)
    fig.savefig(figure_path, facecolor='white', transparent=False, dpi=300, bbox_inches='tight')





