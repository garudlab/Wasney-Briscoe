import tarfile, bz2, cPickle, os, sys

# MIDAS STuff
import config
import parse_midas_data
import numpy as np
from numpy.random import shuffle
import stats_utils
import diversity_utils
import core_gene_utils

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument("--species", help="Name of specific species to run code on")
parser.add_argument("-o", "--outdir", help="Where to write output file",metavar="DIR")
parser.add_argument("-Lmax", "--downsample", type=int, help="The number of SNPs above which downsampling will occur",default=1e09) #originally 1e8
parser.add_argument("-fstar", "--freq-threshold", type=float, help="Frequency has to exceed to be included",default=0.05)
parser.add_argument("--fraction-covered", type=float, help="Fraction of samples with sufficient coverage",default=0.3)
parser.add_argument("--fraction-var", type=float, help="Fraction of samples with sufficient # snvs",default=0.2)

parser.add_argument("--min-depth", type=float, help="Sufficient coverage",default=5)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
species_name = args.species
outdir = args.outdir + species_name + "/"
max_snps = args.downsample
fstar = args.freq_threshold
fraction_covered = args.fraction_covered 
min_depth = args.min_depth
fraction_var = args.fraction_var

print(fraction_covered)

snps_directory= config.data_directory + "%s/snps/" #MW: changed this to work with normal microbiome_evolution scripts

sys.stderr.write("Loading core genes...\n")
core_genes = core_gene_utils.parse_core_genes(species_name)
personal_core_genes = core_gene_utils.parse_personal_core_genes(species_name)
non_shared_genes = personal_core_genes
shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
sys.stderr.write("Done! %d core genes and %d shared genes and %d non-shared genes\n" % (len(core_genes), len(shared_pangenome_genes), len(non_shared_genes)))
   
snp_file = outdir+species_name+".strainfinder.p"
location_file = outdir+species_name+".strainfinder.locations.p"
sample_file = outdir+species_name+".strainfinder.samples.p"

#Make directories if needed
if os.path.exists(args.outdir):
    if os.path.exists(outdir):
        print "Strainfinder directory for " + species_name + " already exists."
    else: 
        print "Making strainfinder directory for " + species_name + "."
        os.mkdir(outdir)
else:
    os.mkdir(args.outdir)
    os.mkdir(outdir)

#snp_samples = parse_timecourse_data.morteza_samples
snp_samples = np.loadtxt(config.accessions, comments="#", delimiter="\n", dtype = "str", unpack=False)
sample_size = len(snp_samples)
sys.stderr.write("Proceeding with %d temporal samples!\n" % sample_size)

snp_alignment = [] # (construct a # sites x # samples x # bases (4) array)   
snp_locations = []

final_line_number = 0
while final_line_number >= 0:
    
    print "Final line number: " + str(final_line_number)

    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    
    snp_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_polarized_snps(species_name,debug=debug,allowed_variant_types=set(['1D','2D','3D','4D']),chunk_size=chunk_size,allowed_samples=snp_samples, initial_line_number=final_line_number, allowed_genes=non_shared_genes, snps_directory=snps_directory) #MW: this is from ricky's version (find in the modified repo). 
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
   # import cPickle
    #with open('fprau_allele_counts.p', 'wb') as handle:
    #    cPickle.dump(allele_counts_map, handle)
    #    print("pickle written \n")
    #snp_samples = dummy_samples
    #with open('fprau_allele_counts.p', 'rb') as f:
    #    allele_counts_map = cPickle.load(f) 
        
    for gene_name in allele_counts_map.keys():
    
        for var_type in allele_counts_map[gene_name].keys():
            locations = allele_counts_map[gene_name][var_type]['locations']
            allele_counts = allele_counts_map[gene_name][var_type]['alleles']
            polarizations = allele_counts_map[gene_name][var_type]['polarizations']
                
            if len(allele_counts)==0:
                continue
                    
            depths = allele_counts.sum(axis=2)
            freqs = allele_counts[:,:,0]*1.0/(depths+(depths==0))
                
            for snp_idx in xrange(0,len(locations)):
                good_idxs = (depths[snp_idx,:]>min_depth)
                insufficient_coverage = (good_idxs.sum() < fraction_covered*depths.shape[1]) #Does this snv have good coverage (over 5) in greater than or equal to 0.3 of samples?
                
                masked_freqs = freqs[snp_idx][good_idxs]
                #low_frequency = numpy.mean(numpy.logical_or(masked_freqs<fstar,masked_freqs>1-fstar)) > 1 - fraction_var
                low_frequency = ((masked_freqs<fstar).all() or (masked_freqs>(1-fstar)).all()) #for samples that have over 5 coverages at this snv, identify if all have below 0.05 frequency
                #low_frequency = ((numpy.mean(masked_freqs<fstar) > fraction_var) or numpy.mean((masked_freqs>(1-fstar))) > fraction_var)

                if insufficient_coverage or low_frequency: #if either insufficient_coverage or low_frequency, eliminate this snv
                    continue
                 
                four_base_counts = [np.hstack([allele_counts[snp_idx,sample_idx,:]*good_idxs[sample_idx], [0,0]]) for sample_idx in xrange(0,depths.shape[1])]   
                
                snp_alignment.append( four_base_counts )
                #np.append(snp_alignment, four_base_counts)
                
                snp_locations.append( (locations[snp_idx][0], locations[snp_idx][1], polarizations[snp_idx] ) )
                #np.append(snp_locations, (locations[snp_idx][0], locations[snp_idx][1], polarizations[snp_idx] ))

  
    if len(snp_alignment) > max_snps:
        print "Shuffling!"
        shuffle(snp_alignment)
        snp_alignment = snp_alignment[0:max_snps]   
    
    print len(snp_alignment)             
                    
    snp_alignment = np.array(snp_alignment)
    #sys.stderr.write("Original aligmmetn: %s \n" % (str(snp_alignment.shape)))
    snp_alignment = np.swapaxes(snp_alignment,0,1)
    sys.stderr.write("Saving %s alignment file \n" % (str(snp_alignment.shape)))
    cPickle.dump(snp_alignment, open(snp_file, 'wb'))
    cPickle.dump(snp_locations, open(location_file, 'wb'))
    cPickle.dump(snp_samples, open(sample_file, 'wb'))
