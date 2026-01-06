import numpy
import sys
import config
import gzip
import os.path
import os
import midas_db_utils
import parse_midas_data #MW: Added 09/21/23

core_genes_directory = ("%score_genes/" % config.data_directory)
external_core_genes_directory = ("%score_genes/external/" % config.data_directory)

default_external_shared_gene_filename = (external_core_genes_directory+"shared_genes.txt.gz")
default_external_core_gene_filename = (external_core_genes_directory+"core_genes.txt.gz")
default_external_stringent_core_gene_filename = (external_core_genes_directory+"core_genes_stringent.txt.gz")
default_external_gene_freq_template = (external_core_genes_directory+"%s_gene_freqs.txt.gz")

default_shared_gene_filename = (core_genes_directory+"shared_genes.txt.gz")
default_core_gene_filename = (core_genes_directory+"core_genes.txt.gz")
default_stringent_core_gene_filename = (core_genes_directory+"core_genes_stringent.txt.gz")
default_gene_freq_template = (core_genes_directory+"%s_gene_freqs.txt.gz")

#HMP data #MW added on 11/2/23
shared_gene_HMP_filename = (core_genes_directory+"shared_genes_HMP.txt.gz")
core_gene_HMP_filename = (core_genes_directory+"core_genes_HMP.txt.gz")
stringent_core_gene_HMP_filename = (core_genes_directory+"core_genes_stringent_HMP.txt.gz")


def parse_core_genes(desired_species_name="", core_gene_filename=default_core_gene_filename, external_core_gene_filename=default_external_core_gene_filename, external_filtering=True):
    
    core_genes = set()
    core_gene_file = gzip.open(core_gene_filename,"rt", encoding="utf-8") 
    for line in core_gene_file:
        
        items = line.split(":")
        if len(items)<2:
            continue
            
        species_name = items[0].strip()
        gene_names = [subitem.strip() for subitem in items[1].split(",")]
        
        if (species_name==desired_species_name) or (desired_species_name==""):
            core_genes.update(gene_names)
            
    core_gene_file.close() 
    
    external_core_genes = set()
    if os.path.isfile(external_core_gene_filename):
        
        external_core_gene_file = gzip.open(external_core_gene_filename,"rt", encoding="utf-8") 
        
        for line in external_core_gene_file:
        
            items = line.split(":")
            if len(items)<2:
                continue
            
            species_name = items[0].strip()
            gene_names = [subitem.strip() for subitem in items[1].split(",")]
        
            if (species_name==desired_species_name) or (desired_species_name==""):
                external_core_genes.update(gene_names)
            
        external_core_gene_file.close() 
    
    if external_filtering and len(external_core_genes)>0:
        # some externally provided core genes
        core_genes = (core_genes & external_core_genes)
        
    return core_genes
    
    
def parse_shared_genes(desired_species_name="", shared_gene_filename=default_shared_gene_filename, external_shared_gene_filename=default_external_shared_gene_filename, external_filtering=True):
    
    shared_genes = set()
    shared_gene_file = gzip.open(shared_gene_filename,"rt", encoding="utf-8") 
    for line in shared_gene_file:
        items = line.split(":")
        if len(items)<2:
            continue
            
        species_name = items[0].strip()
        gene_names_str = items[1].strip()
        if gene_names_str.startswith('N/A'): # Wasn't enough pangenome data to detect shared genes
            gene_names = []
        else:
            gene_names = [subitem.strip() for subitem in gene_names_str.split(",")]
        
        if (species_name==desired_species_name) or (desired_species_name==""):
            shared_genes.update(gene_names)
            
    shared_gene_file.close() 
    
    external_shared_genes = set()
    if os.path.isfile(external_shared_gene_filename):
        
        external_shared_gene_file = gzip.open(external_shared_gene_filename,"rt", encoding="utf-8") 
            
        for line in external_shared_gene_file:
        
            items = line.split(":")
            if len(items)<2:
                continue
            
            species_name = items[0].strip()
            gene_names_str = items[1].strip()
            if gene_names_str.startswith('N/A'): # Wasn't enough pangenome data to detect shared genes
                gene_names = []
            else:
                gene_names = [subitem.strip() for subitem in gene_names_str.split(",")]
        
            if (species_name==desired_species_name) or (desired_species_name==""):
                external_shared_genes.update(gene_names)
            
        external_shared_gene_file.close() 
    
    if external_filtering and len(external_shared_genes)>0:
        # some externally provided core genes
        shared_genes = (shared_genes | external_shared_genes)
        
    return shared_genes

def parse_non_shared_reference_genes(desired_species_name="", shared_gene_filename=default_shared_gene_filename, external_shared_gene_filename=default_external_shared_gene_filename, external_filtering=True):
    import parse_midas_data
    shared_genes = parse_shared_genes(desired_species_name, shared_gene_filename, external_shared_gene_filename, external_filtering)
    reference_genes = parse_midas_data.load_reference_genes(desired_species_name)
    non_shared_reference_genes = set(reference_genes)-shared_genes
    return non_shared_reference_genes
    
def get_good_pangenome_samples(species_name, marker_coverages, gene_copynum_matrix): #MW 09/03/2025: added species_name as an argument

    cmin = config.core_genome_min_copynum
    cmax = config.core_genome_max_copynum  

    # Load reference genes
    num_reference_genes = len(parse_midas_data.load_reference_genes(species_name))
    
    num_present_genes = (gene_copynum_matrix>cmin).sum(axis=0)
    num_high_genes = (gene_copynum_matrix>cmax).sum(axis=0)
    
    good_sample_idxs = (num_present_genes>0.3*num_reference_genes)*(num_high_genes<0.3*num_present_genes)
    return good_sample_idxs

# def get_good_pangenome_samples(marker_coverages, gene_copynum_matrix, species_name): #THIS VERSION ALTERED 09/21/23; go back to old version when done

#     cmin = config.core_genome_min_copynum
#     cmax = config.core_genome_max_copynum  

#     # Load reference genes
#     num_reference_genes = len(parse_midas_data.load_reference_genes(species_name))
    
#     num_present_genes = (gene_copynum_matrix>cmin).sum(axis=0)
#     num_high_genes = (gene_copynum_matrix>cmax).sum(axis=0)
    
#     good_sample_idxs = (num_present_genes>0.3*num_reference_genes)*(num_high_genes<0.3*num_present_genes)
#     return good_sample_idxs


def parse_gene_freqs(desired_species_name, use_external=False):
    
    if use_external:
        filename_template = default_external_gene_freq_template
    else:
        filename_template = default_gene_freq_template
    
    
    filename = filename_template % (desired_species_name)
    if not os.path.isfile(filename):
        return None
        
    file = gzip.open(filename,"rt", encoding="utf-8") 
    gene_freq_map = {}
    for line in file:
        items = line.split()
        gene_name = items[0]
        f = float(items[1])
        gene_freq_map[gene_name] = f
    file.close()
    
    return gene_freq_map

#MW: This combines the results of parse_core_genes() and parse_non_shared_reference_genes()
def parse_personal_core_genes(desired_species_name):
    unfiltered_core_genes = parse_core_genes(desired_species_name, external_filtering=False)
    non_shared_genes = parse_non_shared_reference_genes(desired_species_name,external_filtering=True)
    
    personal_core_genes = (unfiltered_core_genes & non_shared_genes)
    
    #if len(personal_core_genes)==0:
        #sys.stderr.write("Warning: zero overlap between %d core genes and %d non_shared genes!\n" % (len(unfiltered_core_genes), len(non_shared_genes)))
    
    return personal_core_genes    
            
    
# Actually calculate the core genes
if __name__=='__main__':
    
    #PARSE ARGUMENTS
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--use_HMP", help="Determine whether to add to existing HMP results", action="store_true")
    
    args = parser.parse_args()

    use_HMP = args.use_HMP
    
    import parse_midas_data
    
    os.system('mkdir -p %s' % core_genes_directory)
    os.system('mkdir -p %s' % external_core_genes_directory)
    
    # pangenome_species = parse_midas_data.parse_good_species_list()
    pangenome_species = parse_midas_data.parse_species_list() # MW 08/29/2025: loaded all species in the species_snps.txt file
 
    cmin = config.core_genome_min_copynum
    cmax = config.core_genome_max_copynum  
    shared_cmin = config.shared_genome_min_copynum
    
    min_good_fraction = config.core_genome_min_prevalence
    min_coverage = 5 # (for assessing core genome, we'll use a lower coverage value than when we look at real changes)
    
    #Loading files to write to
    output_filename = default_core_gene_filename 
    output_file = gzip.open(output_filename, mode = "wt", encoding = "utf-8") # MW 07/09/2025: changing gzip.GzipFile to gzip open, mode to "wt", and added encoding = "utf-8"

    stringent_output_filename = default_stringent_core_gene_filename
    stringent_output_file = gzip.open(stringent_output_filename,mode = "wt", encoding="utf-8") # MW 07/09/2025: changing gzip.GzipFile to gzip open, mode to "wt", and added encoding = "utf-8"
    
    shared_output_file = gzip.open(default_shared_gene_filename, mode = "wt", encoding="utf-8") # MW 07/09/2025: changing gzip.GzipFile to gzip open, mode to "wt", and added encoding = "utf-8"
    
    for species_name in pangenome_species:
      
        # Load reference genes
        sys.stderr.write("Loading genes on reference genome..\n")
        reference_genes = midas_db_utils.load_reference_genes(species_name)
        sys.stderr.write("Done!\n")

        # Load reference genes
        sys.stderr.write("Loading shared genes from midas db..\n")
        # midas_shared_genes = midas_db_utils.parse_midas_shared_genes(species_name)  # MW 07/08/2025: no cross_species_centroids.txt.gz file            
        sys.stderr.write("Done!\n")


        bad_pangenome_data = False
                        
        # Load gene coverage information for species_name
        sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
        gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
        sys.stderr.write("Done!\n")  
        
        if len(marker_coverages)==0:
            bad_pangenome_data = True
        else:        
            
            high_coverage_idxs = (marker_coverages>=min_coverage)

            if high_coverage_idxs.sum() < 0.5:
                bad_pangenome_data = True

        if bad_pangenome_data:
            # Just use reference genes
            sys.stderr.write("Bad pangenome data for %s!\n" % species_name)
            # shared_gene_names = sorted(midas_shared_genes) # MW 07/08/2025: no cross_species_centroids.txt.gz, so not using this
            # core_gene_names = sorted(reference_genes - midas_shared_genes)
            core_gene_names = sorted(reference_genes)
            # stringent_gene_names = sorted(reference_genes - midas_shared_genes)
            stringent_gene_names = sorted(reference_genes)
      
        else:    
        
            gene_names = numpy.array(gene_names)
            gene_samples = gene_samples[high_coverage_idxs]
            marker_coverages = marker_coverages[high_coverage_idxs]
            gene_depth_matrix = gene_depth_matrix[:,high_coverage_idxs] 
            gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))

            good_sample_idxs = get_good_pangenome_samples(species_name, marker_coverages, gene_copynum_matrix) # MW 09/03/2025: added species name
            bad_sample_idxs = numpy.logical_not(good_sample_idxs)
                    
            #sys.stderr.write("%d bad samples!\n" % bad_sample_idxs.sum())
            
            gene_samples = gene_samples[good_sample_idxs]
            marker_coverages = marker_coverages[good_sample_idxs]
            gene_copynum_matrix = gene_copynum_matrix[:,good_sample_idxs]
            
            reference_gene_idxs = numpy.array([gene_name in reference_genes for gene_name in gene_names])
            
            # midas_shared_idxs = numpy.array([gene_name in midas_shared_genes for gene_name in gene_names]) # MW 07/09/2025: no cross_species_centroids.txt.gz
            
        
            # These are genes that have coverage >=3x normal in some sample. This are candidates for being linked to another species.
            # (they could also be multi-copy genes, but we can't look at much on these genes anyway, so might as well toss them out)
            metagenome_shared_idxs = ((gene_copynum_matrix>shared_cmin).sum(axis=1)>0.5)
            
            # Now union with those we identified from midas db
            # shared_idxs = numpy.logical_or(metagenome_shared_idxs, midas_shared_idxs) 
            shared_idxs = metagenome_shared_idxs # MW 07/09/2025: no cross_species_centroids.txt.gz
            non_shared_idxs = numpy.logical_not(shared_idxs)
            
            shared_gene_names = gene_names[shared_idxs]
            
            # calculating good genes
            good_idxs = (((gene_copynum_matrix>=cmin)*(gene_copynum_matrix<=cmax)).sum(axis=1)*1.0/len(marker_coverages) >= min_good_fraction) #indexes in which gene copy numbers are between 0.3 and 3 in 90% of samples
            core_gene_idxs = good_idxs*reference_gene_idxs*non_shared_idxs
            core_gene_names = gene_names[core_gene_idxs]
            
            # num_metagenome_and_midas = numpy.logical_and(midas_shared_idxs, metagenome_shared_idxs).sum() # MW 07/10/2025: no cross_species_centroids.txt.gz, and  midas_shared_idxs is created using that
            num_metagenome_and_midas = metagenome_shared_idxs.sum()
            # num_metagenome_only = numpy.logical_and(metagenome_shared_idxs, numpy.logical_not(midas_shared_idxs)).sum()
            num_metagenome_only = metagenome_shared_idxs.sum()
            # num_midas_only = numpy.logical_and(midas_shared_idxs, numpy.logical_not(metagenome_shared_idxs)).sum()
            num_midas_only = "NA"
            num_metagenome_or_midas = shared_idxs.sum()
            num_remaining = non_shared_idxs.sum()
            num_reference_remaining = (non_shared_idxs*reference_gene_idxs).sum()
            num_core = core_gene_idxs.sum()
            
            print("%s %d %d %s %d %d %d %d" % (species_name, num_metagenome_and_midas, num_metagenome_only, num_midas_only, num_metagenome_or_midas, num_remaining, num_reference_remaining, num_core)) # MW 07/10/2025: Made the 4th a string (because it's now NA, see above)
            
            # Measure frequencies and output them
            gene_prevalence_numerators = ((gene_copynum_matrix>=cmin)*(gene_copynum_matrix<=cmax)).sum(axis=1)
            gene_prevalence_denominators = ((gene_copynum_matrix<=cmax).sum(axis=1))
            
            good_prevalence_idxs = (gene_prevalence_numerators>0.5)*(gene_prevalence_denominators>0.5)*non_shared_idxs
            
            gene_prevalence_names = gene_names[good_prevalence_idxs]
            gene_prevalences = gene_prevalence_numerators[good_prevalence_idxs]*1.0/gene_prevalence_denominators[good_prevalence_idxs]
            
            
            
            gene_freq_output_file = gzip.open(default_gene_freq_template % species_name, mode="wt", encoding="utf-8") # MW 07/09/2025: encoding as utf by changing gzip.GzipFile(default_gene_freq_template % species_name, "w") to gzip.open(default_gene_freq_template % species_name, mode="wt", encoding="utf-8")
            for gene_name, f in zip(gene_prevalence_names, gene_prevalences):
                gene_freq_output_file.write("%s %g\n" % (gene_name, f)) #MW: I want to rerun it with this
            gene_freq_output_file.close()
    
            
            # calculating good genes w/ stringent definition (100%)
            bad_idxs = (gene_copynum_matrix<config.gainloss_max_absent_copynum).sum(axis=1) > 0.5 #"absent" in 50% of samples
            good_idxs = numpy.logical_not(bad_idxs)
            stringent_gene_names = gene_names[good_idxs*reference_gene_idxs*non_shared_idxs]
            #sys.stderr.write("%d stringent core genes out of %d\n" % (len(stringent_gene_names), len(gene_names))) 
            
            #MW 11/2/23 addition: Synthesize results with previously determined core/shared genes using HMP data
            
            if use_HMP:
                ### Load HMP data
                core_genes_HMP = parse_core_genes(species_name,core_gene_filename=core_gene_HMP_filename)
                stringent_core_genes_HMP = parse_core_genes(species_name,core_gene_filename=stringent_core_gene_HMP_filename)
                shared_genes_HMP = parse_shared_genes(species_name, shared_gene_HMP_filename, default_external_shared_gene_filename)
                ### Determine if there's anything in the HMP files for that species. 
                    ### If not, just add the genes that were calculated in this script for that species
                ### In the case of core genes, if there is something present in HMP, subtract any
                ### gene we calculate as being shared.
                ### Treat stringent_core_genes the same. 
                ### However, with shared genes, if HMP has nothing for a species, add what was calculated here.
                ### If HMP does have shared genes for that species, add any genes that aren't already present
                if len(core_genes_HMP) > 0:
                    final_core_genes = [gene_name for gene_name in core_genes_HMP if gene_name not in shared_gene_names]
                    output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in final_core_genes])))
                else:
                    output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in core_gene_names])))

                if len(stringent_core_genes_HMP) > 0:
                    final_stringent_core_genes = [gene_name for gene_name in stringent_core_genes_HMP if gene_name not in shared_gene_names]
                    stringent_output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in final_stringent_core_genes])))
                else:
                    stringent_output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in stringent_gene_names])))

                if len(shared_genes_HMP) > 0:
                    final_shared_genes = shared_genes_HMP
                    final_shared_genes.update(set(shared_gene_names))
                    final_shared_genes = numpy.array(list(final_shared_genes))
                    shared_output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in final_shared_genes])))
                else:
                    shared_output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in shared_gene_names])))


            else:
                # Write output to file!
                shared_output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in shared_gene_names])))
                output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in core_gene_names])))
                stringent_output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in stringent_gene_names])))
    
                

    shared_output_file.close()
    output_file.close()
    stringent_output_file.close()

    
    