#!/bin/bash
#$ -N bzip_MIDAS_files
#$ -l h_data=32G,h_rt=168:00:00,highp
#$ -e ~/Wasney-Briscoe-2024/scripts/MIDAS/errors/
#$ -o ~/Wasney-Briscoe-2024/scripts/MIDAS/outputs/



############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: scriptTemplate [-m|--merged_path / -h|--help]"
   echo "options:"
   echo "-m|--merged_path      The directory holding the merged midas outputs."
   echo
   echo "-h|--help             Display this help message and exit script."
   echo
}

############################################################
# Setup                                                    #
############################################################

while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done

#Default arguments
merged_output_path=~/merged_data/


#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -m|--merged_path)
        merged_output_path="$2"
        shift # Remove --merged_path from processing
        shift # Remove value from processing
        ;;
    esac
done


############################################################
# execute                                                  #
############################################################

### Species
cd ${merged_output_path}species

if [ -f count_reads.txt ]; then
   echo "compressing count_reads.txt"
   bzip2 count_reads.txt
elif [ -f count_reads.txt.bz2 ]; then
   echo "count_reads.txt.bz2 already exists"
else
   echo "count_reads file doesn't exist"
fi

if [ -f coverage.txt ]; then
   echo "compressing coverage.txt"
   bzip2 coverage.txt
elif [ -f coverage.txt.bz2 ]; then
   echo "coverage.txt.bz2 already exists"
else
   echo "coverage file doesn't exist"
fi

if [ -f relative_abundance.txt ]; then
   echo "compressing relative_abundance.txt"
   bzip2 relative_abundance.txt
elif [ -f relative_abundance.txt.bz2 ]; then
   echo "relative_abundance.txt.bz2 already exists"
else
   echo "relative_abundance file doesn't exist"
fi

if [ -f species_prevalence.txt ]; then
   echo "compressing species_prevalence.txt"
   bzip2 species_prevalence.txt
elif [ -f species_prevalence.txt.bz2 ]; then
   echo "species_prevalence.txt.bz2 already exists"
else
   echo "species_prevalence file doesn't exist"
fi


# ### genes
cd ${merged_output_path}genes

for dir in */; do
   cd ${dir}
   
   echo "Proccessing genes of " $dir

   if [ -f genes_copynum.txt ]; then
      echo "compressing genes_copynum.txt"
      bzip2 genes_copynum.txt
   elif [ -f genes_copynum.txt.bz2 ]; then
      echo "genes_copynum.txt.bz2 already exists"
   else
      echo "genes_copynum file doesn't exist"
   fi

   if [ -f genes_depth.txt ]; then
      echo "compressing genes_depth.txt"
      bzip2 genes_depth.txt
   elif [ -f genes_depth.txt.bz2 ]; then
      echo "genes_depth.txt.bz2 already exists"
   else
      echo "genes_depth file doesn't exist"
   fi

   if [ -f genes_presabs.txt ]; then
      echo "compressing genes_presabs.txt"
      bzip2 genes_presabs.txt
   elif [ -f genes_presabs.txt.bz2 ]; then
      echo "genes_presabs.txt.bz2 already exists"
   else
      echo "genes_presabs file doesn't exist"
   fi

   if [ -f genes_reads.txt ]; then
      echo "compressing genes_reads.txt"
      bzip2 genes_reads.txt
   elif [ -f genes_reads.txt.bz2 ]; then
      echo "genes_reads.txt.bz2 already exists"
   else
      echo "genes_reads file doesn't exist"
   fi

   #if [ -f genes_summary.txt ]; then
   #   echo "compressing genes_summary.txt"
   #   bzip2 genes_summary.txt
   #elif [ -f genes_summary.txt.bz2 ]; then
   #   echo "genes_summary.txt.bz2 already exists"
   #else
   #   echo "genes_summary file doesn't exist"
   #fi

   cd ..

done


### SNPs

cd ${merged_output_path}snps

for dir in */; do
   cd ${dir}
   
   echo "Proccessing snps of " $dir

   if [ -f snps_alt_allele.txt ]; then
      echo "compressing snps_alt_allele.txt"
      bzip2 snps_alt_allele.txt
   elif [ -f snps_alt_allele.txt.bz2 ]; then
      echo "snps_alt_allele.txt.bz2 already exists"
   else
      echo "snps_alt_allele file doesn't exist"
   fi

   if [ -f snps_depth.txt ]; then
      echo "compressing snps_depth.txt"
      bzip2 snps_depth.txt
   elif [ -f snps_depth.txt.bz2 ]; then
      echo "snps_depth.txt.bz2 already exists"
   else
      echo "snps_depth file doesn't exist"
   fi

   if [ -f snps_info.txt ]; then
      echo "compressing snps_info.txt"
      bzip2 snps_info.txt
   elif [ -f snps_info.txt.bz2 ]; then
      echo "snps_info.txt.bz2 already exists"
   else
      echo "snps_info file doesn't exist"
   fi

   if [ -f snps_ref_freq.txt ]; then
      echo "compressing snps_ref_freq.txt"
      bzip2 snps_ref_freq.txt
   elif [ -f snps_ref_freq.txt.bz2 ]; then
      echo "snps_ref_freq.txt.bz2 already exists"
   else
      echo "snps_ref_freq file doesn't exist"
   fi

   cd ..

done
