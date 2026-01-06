#!/bin/bash
#$ -N strainfinder_preprocessing
#$ -l h_data=8G,h_rt=2:00:00,highp
#$ -o /u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/postprocessing/outputs/
#$ -e /u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/postprocessing/errors/
#$ -t 1-117:1
#$ -js 100000000000000000000000000000000000000000000

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-s/--species_list | -d|--scripts_dir | -o/--outdir | -h/--help]"
   echo "options:"
   echo "-s|--species_list     The file used to select which species to process."
   echo
   echo "-d|--scripts_dir     The file used to select which species to process."
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
species_file=/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/metadata/species_snps.txt
outdir=/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/strain_phasing/input/
scripts_dir=/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/postprocessing/postprocessing_scripts/

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -s|--species_list)
        species_file="$2"
        shift # Remove --species list from processing
        shift # Remove value from processing
        ;;
        -d|--scripts_dir)
        scripts_dir="$2"
        shift 
        shift 
        ;;
        -o|--outdir)
        outdir="$2"
        shift
        shift 
        ;;
    esac
done

############################################################
# Main program                                             #
############################################################

source /u/local/Modules/default/init/modules.sh
module load miniforge
conda activate py38


i=$((SGE_TASK_ID))

species_name=$(sed -n "$i"p $species_file)

mkdir -p $outdir

cd $scripts_dir

python create_StrainFinderInput.py --species $species_name --outdir $outdir

