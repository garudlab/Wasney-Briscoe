import sys
import os
import bz2
import parse_midas_data

if len(sys.argv) > 1:
    species_name=sys.argv[1]
else:
    species_name=parse_midas_data.debug_species_name

sys.stderr.write("Calculating pvalues for %s...\n" % species_name)

output_filename = "%ssnps/%s/annotated_snps.txt.bz2" % (parse_midas_data.data_directory, species_name)

os.system('python %spostprocessing_scripts/pipe_midas_data.py %s | %spostprocessing_scripts/error_filtering_cpp_code/annotate_pvalue --disabled | bzip2 -c > %s' % (parse_midas_data.scripts_directory, species_name, parse_midas_data.scripts_directory, output_filename) )  #changed the directory of annotate_pvalue (original script had it in main microbiome_evolution folder
 
