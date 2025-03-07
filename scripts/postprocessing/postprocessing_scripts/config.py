###############################################################################
#
# Set up default source and output directories
#
###############################################################################
import os.path 
from math import log10

#CHANGE FOR EACH PROJECT
data_directory = os.path.expanduser("~/merged_data/")
project_folder = os.path.expanduser("~/")
HMP_data_directory = os.path.expanduser("~/Wasney-Briscoe/scripts/postprocessing/")
analysis_directory = os.path.expanduser("~/Wasney-Briscoe/") #change this for your project
metadata_directory = os.path.expanduser("~/Wasney-Briscoe/metadata/")
accessions = os.path.expanduser("~/Wasney-Briscoe/scripts/accessions.txt")

#output
figure_directory = os.path.expanduser("~/Wasney-Briscoe/figures/")

#STATIC
scripts_directory = os.path.expanduser("~/Wasney-Briscoe/scripts/postprocessing/")
patric_directory = os.path.expanduser("~/PATRIC/")
midas_directory = os.path.expanduser("~/midas_db_v1.2/")

# We use this one to debug because it was the first one we looked at
debug_species_name = 'Bacteroides_uniformis_57318'

good_species_min_coverage = 10
good_species_min_prevalence = 2 #MW: changing this from 10 to 2

min_median_coverage = 20
min_site_coverage = 20 #This is the minimum coverage a site needs to have to be considered in a change.

consensus_lower_threshold = 0.2
consensus_upper_threshold = 0.8

# consensus_lower_threshold = 0.4
# consensus_upper_threshold = 0.6

fixation_min_change = consensus_upper_threshold-consensus_lower_threshold
fixation_log10_depth_ratio_threshold = log10(3)

threshold_within_between_fraction = 0.1
threshold_pi = 1e-03

min_opportunities = 100000

modification_difference_threshold = 20
replacement_difference_threshold = 500

twin_modification_difference_threshold = 1000
twin_replacement_difference_threshold = 1000

gainloss_max_absent_copynum = 0.05
gainloss_min_normal_copynum = 0.6
gainloss_max_normal_copynum = 1.2

core_genome_min_copynum = 0.3
core_genome_max_copynum = 3 # BG: should we use a maximum for "core genome"? I'm going to go w/ yes for now
core_genome_min_prevalence = 0.9
shared_genome_min_copynum = 3

# Default parameters for pipe snps
# (Initial filtering for snps, done during postprocessing)
pipe_snps_min_samples=1 #changed this from 4 to 1
pipe_snps_min_nonzero_median_coverage=5
pipe_snps_lower_depth_factor=0.3
pipe_snps_upper_depth_factor=3

parse_snps_min_freq = 0.05

between_host_min_sample_size = 33
between_host_ld_min_sample_size = 10
within_host_min_sample_size = 3
within_host_min_haploid_sample_size = 10

between_low_divergence_threshold = 2e-04

# Comment this out
from parse_midas_data import *
# and uncomment this
#from parse_simulated_data import *
# for isolate data
