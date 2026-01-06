# Setting working directory
import sys
sys.path.insert(0, "/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/scripts/strain_phasing/")

# Normal Libraries
import pandas as pd
import numpy as np
import pickle
import scipy.stats
import random as rand
from random import randint,sample,choices
from math import log
import os 
from datetime import datetime

# Plotting libraries
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib import cm
import figure_utils as fu
plasma_cmap = cm.get_cmap('plasma')

# config
import config

# predefined functions
from strain_phasing_functions import *
from microbiome_evolution_functions import *

# Variables

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-s", "--species", type=str, help="The list of species to use", default = "Bacteroides_vulgatus_57955")

args = parser.parse_args()

species = args.species

sys.stderr.write("Species: %s \n\n" % (species))

strainfinder_dir = "%sinput" % (config.strain_phasing_directory)

#Raw cluster
raw_cluster_path = "%s%s" % (config.strain_phasing_directory, "strain_clusters_ReadSupport/")
species_raw_cluster_dir = "%s%s/" % (raw_cluster_path, species)
if not os.path.exists(species_raw_cluster_dir):
    os.makedirs(species_raw_cluster_dir)
    print("Cluster directory created successfully!")
else:
    print("Cluster directory already exists.")

strain_phasing_figures_dir = "%s%s" % (config.figure_directory, "strain_phasing/")
if not os.path.exists(strain_phasing_figures_dir):
    os.makedirs(strain_phasing_figures_dir)
    print("Figure directory created successfully!")
else:
    print("Figure directory already exists.")

# Loading metadata parsers

sample_metadata_map = parse_sample_metadata_map()
subject_sample_map = parse_subject_sample_map()

subjects = subject_sample_map.keys()

# Clustering parameters
## Meta-parameters: experiment with these—no hard and fast rules!

## minimum number of SNVs which need to be clustered together in order to qualify as a "strain"
min_cluster_size = 1000

## minimum fraction of sites which pass our coverage threshold which must be in a cluster in order for it to qualify 
## as a strain
min_cluster_fraction = 1/10

## For computational efficiency, we can downsample the SNVs we actually perform strain phasing on
max_num_snvs = 20000

## distance threshold to be considered linked—lower means trajectories have to be more   
max_d = 3.5

## minimum coverage to consider allele frequency at a site for purposes of clustering
min_coverage = 10 
# min_coverage = 5

## minimum average sample coverage at polymorphic sites (e.g. sites in the A/D matrices)
# min_sample_coverage = 5
min_sample_coverage = 5

## polymorphic & covered fraction: what percentage of samples does a site need 
## with coverage > min_coverage and polymorphic to be included in downstream analyses? 
poly_cov_frac = 1/5 #

## Number of clusters to calculate
n_clusters = 100

#Minimum number of snvs per sample
# min_num_snvs_per_sample = 100
min_num_snvs_per_sample = 100


for subject_id in subjects:

    sys.stderr.write("\nProcessing subject %s \n\n" % (subject_id))

    samples_of_interest = list(subject_sample_map[subject_id].keys())

    Fs,Ass,Dss = return_FAD(species, min_coverage=min_coverage, 
                            min_sample_coverage=min_sample_coverage, 
                            poly_cov_frac = poly_cov_frac, 
                            calculate_poly_cov_frac=False, 
                            read_support = True,
                            subject_id=subject_id) 
    
    #filter out samples without an adequate number of SNVs when all is said and done
    sample_with_adequate_snv_count = ~((~np.isnan(Fs)).sum() < min_num_snvs_per_sample)

    Fs = Fs.loc[:,sample_with_adequate_snv_count]
    Ass = Ass.loc[:,sample_with_adequate_snv_count]
    Dss = Dss.loc[:,sample_with_adequate_snv_count]

    # Creating distance matrix for all SNVs

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

    # If there are no clusters, skip to next subject
    if len(final_clusters) == 0:
        sys.stderr.write("\nNo %s strain clusters detected in subject %s \n" % (species, subject_id))
        continue
        
    #SAVING RAW FILE
    species_raw_cluster_path = "%s%s%s%s%s" % (species_raw_cluster_dir, species, "_subject_", subject_id, "_RawCluster.pckl")

    pickle_object = open(species_raw_cluster_path, "wb")
    pickle.dump(final_clusters, pickle_object)
    pickle_object.close()

    # Creating polarized clusters

    if len(final_clusters) > 1:
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

    # Filtering out sample that don't have adquate SNVs

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

    if final_clusters[0].shape[1] <= 1:
        sys.stderr.write("\nOnly a single sample detected in subject %s \n" % (subject_id))
        continue

    
    # Plotting

    subset = True
    subset_value = 1000

    if subset:
        high_coverage_snv_idxs = Dss.median(axis = 1).sort_values(ascending = False).index

    # Create the line plot
    fig, ax = plt.subplots(figsize=(20, 8))

    for i,f in enumerate(final_f):
        if subset:
            high_coverage_snv_idxs_strain = high_coverage_snv_idxs.intersection(final_clusters[i].index)[:subset_value]
            sns.lineplot(data = final_clusters[i].loc[high_coverage_snv_idxs_strain].T.values, ax = ax, palette=[["red", "blue", "green"][i]]*subset_value, alpha = 0.075, dashes = False, legend = False)
        else:
            sns.lineplot(data = final_clusters[i].T.values, ax = ax, palette=[["red", "blue", "green"][i]]*final_clusters[i].shape[0], alpha = 0.075, dashes = False, legend = False)
        # sns.lineplot(data = f.values, ax = ax, color = ["red", "blue", "green"][i], linewidth = 4)
        sns.lineplot(data = pd.DataFrame(f), x = "sample", y = 0,  ax = ax, color = "black", linewidth = 9)
        sns.lineplot(data = pd.DataFrame(f), x = "sample", y = 0,  ax = ax, color = ["red", "blue", "green"][i], linewidth = 7)

        

    # Creating x axis minor ticks and extracting locations for vspan; creating major tick labels
    major_ticks = []
    major_tick_labels = []
    minor_ticks = []
    minor_tick_labels = []
    time_point = ""
    x_ticks_loc = ax.get_xticks()
    vspan_counter = 0
    vspan_vec = []
    for i, column in enumerate(final_clusters[0].columns):
        # Major ticks: timepoint
        new_time_point = "\n\n" +column[1] + ",\n" + column[2]
        if (time_point != new_time_point) & (i != len(final_clusters[0].columns) - 1):
            time_point = new_time_point
            # major_ticks.append(x_ticks_loc[i])
            major_tick_labels.append(time_point)

            # add vspan
            vspan_counter += 1

            if vspan_counter == 1:
                xmin = 0
            else: 
                xmax = x_ticks_loc[i] - 0.5
                vspan_vec.append([xmin, xmax])
                xmin = xmax
        elif (time_point != new_time_point) & (i == len(final_clusters[0].columns) - 1):  
            time_point = new_time_point    
            major_tick_labels.append(time_point)
            
            xmax = x_ticks_loc[i] - 0.5
            vspan_vec.append([xmin, xmax])
            xmin = xmax
            xmax = x_ticks_loc[i]
            vspan_vec.append([xmin, xmax])
        elif (i == len(final_clusters[0].columns) - 1):
            xmax = x_ticks_loc[i]
            vspan_vec.append([xmin, xmax])

        # Minor ticks: sample type
        sample_type = column[0]
        minor_ticks.append(x_ticks_loc[i])
        minor_tick_labels.append(sample_type)

    # adding vspan and creating minor ticks
    for i,v in enumerate(vspan_vec):
        # adding vspan
        if i % 2 == 1:
            ax.axvspan(v[0],v[1],alpha=.2,color='grey') 
        
        if (i == 0) & (v[1] == 0.5):
            major_ticks.append(0)
        elif (i == len(vspan_vec) - 1):
            major_ticks.append(np.mean([v[0] + 0.5, v[1]]))
        else:
            major_ticks.append(np.mean(v))


    ax.set_xticks(major_ticks)
    ax.set_xticklabels(major_tick_labels)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticklabels(minor_tick_labels, minor=True)

    ax.xaxis.remove_overlapping_locs = False

    plt.tick_params(axis='x',which='major',bottom=False,left=False,top=False) 

    # Titles

    ax.set_title("%s%s%s%s" % ("Strains of ", species, " in subject ", subject_id), size = 20)
    ax.set_xlabel("Sample", size = 20)
    ax.set_ylabel("Frequency", size = 20)

            
    # Legend

    legend_elements = []

    for i in np.arange(len(final_f)):
        legend_elements.append(Line2D([0], [0], color=["red", "blue", "green"][i], lw=8, label='Strain %s' % (i+1)))

    ax.legend(handles=legend_elements, fontsize = 16)


    # Saving
    plt.tight_layout()


    out_file = "%s%s%s%s%s" % (strain_phasing_figures_dir, species, "_subject_", subject_id, "_StrainFreq_ReadSupport.png")
    fig.savefig(out_file, dpi = 300)









        


        



