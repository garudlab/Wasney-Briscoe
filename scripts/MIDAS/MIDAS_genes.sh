#!/bin/bash
#$ -N MIDAS_genes
#$ -e ~/Wasney-Briscoe/scripts/MIDAS/errors/
#$ -o ~/Wasney-Briscoe/scripts/MIDAS/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=48G
#$ -l time=48:00:00
#$ -l highp
#$ -t 1-41:1

#Load modules
source /u/local/Modules/default/init/modules.sh
module load anaconda3
module load singularity

#Load python 2.7 conda environment
conda activate python27_env

export PYTHONPATH=$PYTHONPATH:~/MIDAS_mod
export PATH=$PATH:~/MIDAS_mod/scripts
export MIDAS_DB=~/midas_db_v1.2

#Accessions
readarray accs < ~/Wasney-Briscoe/scripts/accessions.txt
accs=(null ${accs[@]}) # zero to one start index
acc=${accs[$SGE_TASK_ID]}
echo $acc

OUTDIR=~/midas_output/${acc}

fastq1=~/fastqs/${acc}_1.fq.gz
fastq2=~/fastqs/${acc}_2.fq.gz
species_union=~/species_union/${acc}/species_union.txt

### MIDAS GENES
singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif run_midas.py genes $OUTDIR -1 $fastq1 -2 $fastq2 --extra_species_file $species_union --remove_temp
