import sys
import bz2
import lz4.frame
import numpy
import parse_midas_data

if len(sys.argv) > 1:
    species_name=sys.argv[1]
else:
    species_name=parse_midas_data.debug_species_name

sys.stderr.write("Calculating marker gene coverage for %s...\n" % species_name)

#####
#
# Creates a species-specific marker_coverage.txt.bz2 in the species snp directory
#
#####

depth_file = lz4.frame.open("%ssnps/%s/%s.snps_depth.tsv.lz4" % (parse_midas_data.data_directory, species_name, species_name),"rt", encoding="utf-8") # MW 07/09/2025: decoding bytes like objects in file
    
# get list of samples to use
depth_line = depth_file.readline()
depth_items = depth_line.split()
output_samples = depth_items[1:]
depth_file.close()
    
coverage_file = open("%sspecies/species_marker_coverage.tsv" % (parse_midas_data.data_directory),"r") # MW 07/14/2025: this is no longer a bzipped file, so altering line. In addition, there are two marker coverage files, species_marker_coverage.tsv (mean of the 15 SCGs) and species_marker_median_coverage.tsv median of the 15 SCGs). I chose to use the mean file.

output_coverage_file = bz2.open("%ssnps/%s/marker_coverage.txt.bz2" % (parse_midas_data.data_directory, species_name),"wt", encoding = "utf-8") # MW 07/09/2025: encoding as bytes like objects in file
sys.stderr.write("Writing output to %s" % ("%ssnps/%s/marker_coverage.txt.bz2\n" % (parse_midas_data.data_directory, species_name)))

# get header line
line = coverage_file.readline()

# get list of samples
items = line.split()
samples = items[1:]

# get the indexes of the samples in the coverage file corresponding to the depth file.
output_idxs = numpy.array([samples.index(sample) for sample in output_samples])
    
# write header lines for output file
output_coverage_file.write("\t".join([items[0]]+output_samples))
    
total_depths_defined = False
while True:
            
    # load next lines
    line = coverage_file.readline()
    
    # quit if file has ended
    if line=="":
        break
        
    items = line.split()
        
    current_species = items[0]
    depths = numpy.array([float(item) for item in items[1:]])[output_idxs]
        
    if current_species==species_name:
        # write output
        output_coverage_file.write("\n")
        output_coverage_file.write("\t".join([current_species]+["%g" % d for d in depths]))
        
    #initialize the total depths field
    if not total_depths_defined:
        total_depths = numpy.zeros_like(depths)
        total_depths_defined=True

    # total depths reports the depth across all species. 
    total_depths += depths
        
    
output_coverage_file.write("\n")
output_coverage_file.write("\t".join(["Total"] + ["%g" % d for d in total_depths]))
    
coverage_file.close()
output_coverage_file.close()
    
# Done!
# no return value
