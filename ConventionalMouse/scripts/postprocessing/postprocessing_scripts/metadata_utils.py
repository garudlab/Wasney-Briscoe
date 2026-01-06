import numpy
import sys
import os.path 

import config

###############################################################################
#
# Set up default source and output directories
#
###############################################################################

metadata_directory = "%s%s" % (config.metadata_directory, "metadata_files/") #MW: added 07/30/24

###############################################################################
#
# Parse metadata
#
###############################################################################


##### DICTIONARY OF (1) PATHS TO METADATA FILES (2) VARIABLES INCLUDED IN METADATA
metadata_variable_list = ["sample_id","subject_id","sample_info","batch","sex","original_sample_number","time_point"]









