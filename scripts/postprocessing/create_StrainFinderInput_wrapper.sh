#!/bin/bash
#$ -N strainfinder_preprocessing
#$ -l h_data=16G,h_rt=24:00:00,highp
#$ -e ~/Wasney-Briscoe/scripts/postprocessing/errors/
#$ -o ~/Wasney-Briscoe/scripts/postprocessing/outputs/
#$ -t 1-85:1

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-s/--species_list | -o/--outdir | -h/--help]"
   echo "options:"
   echo "-s|--species_list     The file used to select which species to process."
   echo
   echo "-o|--outdir           The strainfinder output directory."
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
species_file=~/Wasney-Briscoe/scripts/postprocessing/species_snps.txt
outdir=~/strain_phasing/input/

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -s|--species_list)
        species_file="$2"
        shift # Remove --species list from processing
        shift # Remove value from processing
        ;;
        -o|--outdir)
        outdir="$2"
        shift # Remove --species list from processing
        shift # Remove value from processing
        ;;
    esac
done

############################################################
# Main program                                             #
############################################################

source /u/local/Modules/default/init/modules.sh
module load anaconda3
source activate python27_env


i=$((SGE_TASK_ID))

species_name=$(sed -n "$i"p $species_file)

mkdir -p $outdir

cd ~/Wasney-Briscoe/scripts/postprocessing/postprocessing_scripts/

python create_StrainFinderInput.py --species $species_name --outdir $outdir

