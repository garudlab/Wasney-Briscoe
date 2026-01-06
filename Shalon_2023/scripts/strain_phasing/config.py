###############################################################################
#
# Set up default source and output directories
#
###############################################################################
import os.path 
from math import log10

data_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/merged_data/")
project_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/")
snps_directory = os.path.expanduser("%s/snps" % data_directory)
analysis_directory = os.path.expanduser("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/")
scripts_directory = os.path.expanduser("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/scripts/strain_phasing/")
strain_phasing_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/strain_phasing/")
figure_directory = os.path.expanduser("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/figures/")
strain_cluster_directory =  os.path.expanduser("/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/strain_phasing/")
metadata_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/metadata/")

patric_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/software/PATRIC/")
midas_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/midas_db_v1.2/")
# with open('good_species.txt') as f:
#     file = f.readlines()
# good_species = [f.strip() for f in file]
    
debug_species_name = 'Bacteroides_vulgatus_57955'

min_depth = 20

min_SNV_gene = 5
short_range_thresh = 300

good_species_min_coverage = 10
good_species_min_prevalence = 10

min_median_coverage = 20

consensus_lower_threshold = 0.2
consensus_upper_threshold = 0.8
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
core_genome_max_copynum = 3 
core_genome_min_prevalence = 0.9
shared_genome_min_copynum = 3

# Default parameters for pipe snps
# (Initial filtering for snps, done during postprocessing)
pipe_snps_min_samples=4
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
#from parse_HMP_data import *
# and uncomment this
#from parse_simulated_data import *
# for isolate data
