# Calculating $\pi$ for bacterial species

## Step 1. Downsampling nucleotide sites with high coverage in mice 7 and 8

This code produces a midas-format snps_ref_freq.txt.bz2 and snps_depth.txt.bz2 where the columns for samples in Mouse 7 and Mouse 8 are updated on a genomic site-by-genomic site basis if the depth of those samples is above the median depth for samples from Mouse 1 through 6. If the median depth for Mouse 1 through 6 is less than 4, then the sample depth in Mouse 7 and 8 are automatically converted to 0 at that site because a site with less than 4 depth will not be considered in $\pi$ calcultion regardless. 

Input: `snps_ref_freq.txt.bz2` and `snps_depth.txt.bz2`
Output: snps_ref_freq.txt.bz2 and snps_depth.txt.bz2 with modifications to columns corresponding to Mouse 7 and 8. These files can be found in a new directory: 

At each nucleotide site, the script [optimized_site_downsampling.py](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/optimized_site_downsampling.pi) does the following
1. Compute median depth of reads for all samples from Mice 1 through 6
2. If median depth < 4, pass to next site. Else, continue
3. This creates a vector of 0s and 1s, representing the total pool of reference and alternative nucleotides for a given sample, respectively.
2. Based on median

To run fo all species with SNP data, submit the [optimized_site_downsampling.py](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/optimized_site_downsampling.pi) script for multiple species simultaneously by executing the following code from the [Wasney-Briscoe/scripts/pi/](https://github.com/garudlab/Wasney-Briscoe/tree/main/scripts/pi/) directory:

```
while IFS= read -r species;
do
  	qsub -cwd -V -N $species -e logs -o logs -l h_data=4G,time=12:00:00 -b y "./job_script_site_downsampling.sh $species" 
done < ~/Wasney-Briscoe/metadata/species_snps.txt
```

## Step 2. Copying files from normal MIDAS output directory to the downsampled directory

## step 3. Inferring $\pi$
