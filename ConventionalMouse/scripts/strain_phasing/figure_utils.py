
import sys
sys.path.insert(0, "/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/postprocessing/postprocessing_scripts/")
import species_utils

def get_pretty_species_name(species_name, include_number=False):
    
    species_code_map = species_utils.parse_species_code_maps()[0]
    pretty_name = species_code_map[species_name]
        
    return pretty_name
    
def get_abbreviated_species_name(species_name):
    
    items = species_name.split("_")
    
    pretty_name = "%s. %s" % (items[0][0], items[1])
        
    return pretty_name