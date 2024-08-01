import sys
### ADD PATH HERE

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

# Directories
strainfinder_dir = "%sinput" % (config.strain_phasing_directory)

### Raw cluster
raw_cluster_path = "%s%s" % (config.strain_phasing_directory, "strain_clusters/")
species_raw_cluster_dir = "%s%s/" % (raw_cluster_path, species)
if not os.path.exists(species_raw_cluster_dir):
    os.makedirs(species_raw_cluster_dir)
    print("Cluster directory created successfully!")
else:
    print("Cluster directory already exists.")
species_raw_cluster_path = "%s%s%s" % (species_raw_cluster_dir, species, "_RawCluster.pckl")

### Centroid and polarized cluster paths
species_centroid_cluster_path = "%s%s%s" % (species_raw_cluster_dir, species, "_ClusterCentroid.pckl")
species_polarized_cluster_path = "%s%s%s" % (species_raw_cluster_dir, species, "_PolarizedCluster.pckl")

### FINAL Fs path
final_Fs_path = "%s%s%s" % (species_raw_cluster_dir, species, "_final_Fs.pckl")

#MIDAS data
annotated_snps_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species)

