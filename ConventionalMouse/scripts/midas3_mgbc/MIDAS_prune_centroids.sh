#!/bin/bash
#$ -N MIDAS_prune-centroids_step
#$ -e /u/project/ngarud/michaelw/hsiao_data/scripts/MIDASv3/errors/
#$ -o /u/project/ngarud/michaelw/hsiao_data/scripts/MIDASv3/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=12G
#$ -l time=24:00:00
#$ -js 1000000000
#$ -l highp
#$ -t 1-130:1


############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-s|--species_list / -m|--midasdb / -h|--help]"
   echo "options:"
   echo "-s|--species_list        Path to a list of species ids."
   echo
   echo "-m|--midasdb             Path to the MIDASDB."
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
species_list=/u/project/ngarud/michaelw/hsiao_data/metadata/MIDASv3/all_species_list.tsv
midasdb=/u/scratch/m/michaelw/MIDASDB-gtdb

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -s|--species_list)
        species_list="$2"
        shift
        shift
        ;;
        -m|--midasdb)
        midasdb="$2"
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

# read species

readarray species_list < $species_list
species_list=(null ${species_list[@]}) # zero to one start index
species_id=${species_list[$SGE_TASK_ID]}
echo $species_id

# run prune
midas prune_centroids \
    --midasdb_name gtdb \
    --midasdb_dir $midasdb \
    -t 1 \
    --remove_singleton \
    --species $species_id \
    --force \
    --debug 


