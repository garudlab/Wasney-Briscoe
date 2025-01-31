#!/bin/bash
#$ -N cp_files_for_pi
#$ -e ~/Wasney-Briscoe/scripts/pi/errors/
#$ -o ~/Wasney-Briscoe/scripts/pi/outputs/
#$ -cwd
#$ -l h_data=12G
#$ -l time=4:00:00
#$ -t 1-85:1

# Change directory

project_directory=~/
species_file=~/Wasney-Briscoe/metadata/species_snps.txt

# Species
i=$((SGE_TASK_ID))

species_name=$(sed -n "$i"p $species_file)

# Making gene directory, if necessary
mkdir -p ${project_directory}/merged_data/genes/${species_name}/

# Copying snp files
cp ~/merged_data/snps/${species_name}/snps_info.txt.bz2 ${project_directory}/merged_data/snps/${species_name}/.

# Copying gene files

cp ~/merged_data/genes/${species_name}/genes_copynum.txt.bz2 ${project_directory}/merged_data/genes/${species_name}/.
