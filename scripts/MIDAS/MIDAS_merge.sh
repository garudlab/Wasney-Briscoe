#!/bin/bash
#$ -N MIDAS_merge
#$ -e ~/Wasney-Briscoe/scripts/MIDAS/errors/
#$ -o ~/Wasney-Briscoe/scripts/MIDAS/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=32G
#$ -l time=72:00:00
#$ -l highp

#Load modules
source /u/local/Modules/default/init/modules.sh
module load anaconda3
module load singularity

#Load python 2.7 conda environment
conda activate python27_env

#Sources
export PYTHONPATH=$PYTHONPATH:~/MIDAS_mod
export PATH=$PATH:~/MIDAS_mod/scripts
export MIDAS_DB=~/midas_db_v1.2

#Making the merged data directories
mkdir -p ~/merged_data/species/
mkdir -p ~/merged_data/genes/
mkdir -p ~/merged_data/snps/

INDIR=~/midas_output/
OUTDIR=~/merged_data/

####Merging 
#1. species 
#2. genes
#3. snps

#### 1. Species
singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif merge_midas.py species ${OUTDIR}species -i $INDIR -t dir > ${OUTDIR}species/merge_species_record.log

#### 2. Genes
singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif merge_midas.py genes ${OUTDIR}genes -i $INDIR -t dir --sample_depth 1 --min_samples 1 > ${OUTDIR}genes/merge_genes_record.log

#### 3. snps
singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif merge_midas.py snps ${OUTDIR}snps -i $INDIR -t dir  --sample_depth 5 --site_depth 3 --min_samples 1 --site_prev 0.0 > ${OUTDIR}snps/merge_snps_record.log
