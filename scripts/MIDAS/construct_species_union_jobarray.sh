#!/bin/bash
#$ -N MIDAS_species_union_step
#$ -o ~/Wasney-Briscoe-2024/scripts/MIDAS/errors/
#$ -e ~/Wasney-Briscoe-2024/scripts/MIDAS/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=15G
#$ -l time=1:00:00
#$ -l highp
#$ -tc 1 #only 1 task at once
#$ -t 1-41:1


############################################################
# Purpose                                                  #
############################################################

# This script constructs the species union files necessary 
# for MIDAS with multiple samples per host.

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "This program creates species union files necessary for MIDAS processing."
   echo
   echo "Syntax: [-m/--midas_path | -a/--accession_list |  | -h/--help]"
   echo "options:"
   echo "-m|--midas_path      Path to directory including MIDAS processing folders."
   echo
   echo "-a|--accession_list  Path to a list of accessions to run through post-processing."
   echo
   echo "-h|--help            Display this help message and exit script."
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
midas_path=~
accession_list=~/Wasney-Briscoe-2024/scripts/MIDAS/species_union_list.txt


#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -m|--midas_path)
        midas_path="$2"
        shift # Remove --midas_path from processing
        shift # Remove value from processing
        ;;
        -a|--accession_list)
        accession_list="$2"
        shift # Remove --accession_list from processing
        shift # Remove value from processing
        ;;
    esac
done

############################################################
# Main program                                             #
############################################################


#Loading environment
source /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate python_env


# Creating directories, if necessary

if [ -d ${midas_path}species_union/ ]; then
    echo "Species union directory exists"
else
    echo "Making the species union directory."
    mkdir ${midas_path}species_union/ 
fi 

while IFS= read -r accession; do
   if [ -d "${midas_path}species_union/${accession}/" ]; then
      echo "$accession directory exists"
   else
      echo "Making the $accession directory"
      mkdir ${midas_path}species_union/${accession}
   fi
done < <(awk -F ',' 'NR>1 {print $2}' "$accession_list")


#Creating species union
i=0
while read line; do
   i=$((i+1))
   echo "Processing" $line
   if [ $i -eq $SGE_TASK_ID ]; then
      file=$(echo $line | cut -d ',' -f 1)
      echo $file
   fi
done < <(sed '1d' "$accession_list")

python construct_species_union.py $accession_list $file $midas_path

