#!/bin/bash
#$ -N CM-mgbc_merge_species
#$ -e /u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/midas3_mgbc/errors/
#$ -o /u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/midas3_mgbc/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=12G
#$ -l time=6:00:00
#4 -js 1000000000
#$ -l highp


############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-s|--sample_manifest / -m|--midas_output / -f|--fastqs / -h|--help]"
   echo "options:"
   echo "-s|--sample_manifest     Path to a sample manifest file list_of_samples.tsv."
   echo
   echo "-m|--midas_output        Path to the midas_output directory."
   echo
   echo "-f|--fastqs              Path to the fastq directory."
   echo
   echo "-h|--help                Display this help message and exit script."
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
sample_manifest=/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/metadata/list_of_samples.tsv
midas_output=/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc
fastqs=/u/scratch/m/michaelw/ConventionalMouse_fastqs/fastqs_merged

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -s|--sample_manifest)
        sample_manifest="$2"
        shift
        shift
        ;;
        -m|--midas_output)
        midas_output="$2"
        shift
        shift
        ;;
        -f|--fastqs)
        fastqs="$2"
        shift 
        shift 
        ;;
    esac
done

############################################################
# Main program                                             #
############################################################

source /u/local/Modules/default/init/modules.sh
module load miniforge/23.11.0
conda activate midasv3

# # make the sample manifest
# run merge

midas merge_species \
  --samples_list $sample_manifest \
  --min_cov 2 \
  ${midas_output}/merge

