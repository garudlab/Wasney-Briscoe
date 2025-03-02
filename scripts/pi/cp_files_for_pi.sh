#!/bin/bash
#$ -N cp_files_for_pi
#$ -e ~/Wasney-Briscoe/scripts/pi/errors/
#$ -o ~/Wasney-Briscoe/scripts/pi/outputs/
#$ -cwd
#$ -l h_data=12G
#$ -l time=4:00:00
#$ -t 1-85:1

###########################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-s|--species_list / -h|--help]"
   echo "options:"
   echo "-s|--species_list               Path to the list of species to be processed."
   echo
   echo "-p|--project_directory          Path to directory holding MIDAS output directory (merged_data)."
   echo
   echo "-m|--metric                     Median or mean (default is median)."
   echo
   echo "-h|--help                       Display this help message and exit script."
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
species_list=~/Wasney-Briscoe/scripts/postprocessing/species_snps.txt
project_directory=~/
metric=median

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -s|--species_list)
        species_list="$2"
        shift 
        shift 
        ;;
        -p|--project_directory)
        project_directory="$2"
        shift 
        shift 
        ;;
        -m|--mean)
        mean="$2"
        shift 
        shift 
        ;;
    esac
done


############################################################
# execute                                                  #
############################################################

# Species
i=$((SGE_TASK_ID))

species_name=$(sed -n "$i"p $species_file)

# Set the destination directory based on the metric
if [ "$metric" == "median" ]; then
    merge_directory=${project_directory}SiteDownsampled
elif [ "$metric" == "mean" ]; then
    merge_directory=${project_directory}SiteDownsampledByMean
else
    echo "Error: The metric argument is incorrectly defined. Must be 'mean' or 'median'."
    exit 1
fi

# Copying snp files
cp ${project_directory}merged_data/snps/${species_name}/snps_info.txt.bz2 ${merge_directory}/merged_data/snps/${species_name}/.

# Copying gene files

cp ${project_directory}merged_data/genes/${species_name}/genes_copynum.txt.bz2 ${merge_directory}/merged_data/genes/${species_name}/.
