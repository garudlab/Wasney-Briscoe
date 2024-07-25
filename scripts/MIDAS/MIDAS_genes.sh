#!/bin/bash
#$ -N MIDAS_genes
#$ -e ~/Wasney-Briscoe-2024/scripts/MIDAS/errors/
#$ -o ~/Wasney-Briscoe-2024/scripts/MIDAS/outputs/
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

export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS_mod
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS_mod/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2

#Accessions
readarray accs < ~/Wasney-Briscoe-2024/scripts/accessions.txt
accs=(null ${accs[@]}) # zero to one start index
acc=${accs[$SGE_TASK_ID]}
echo $acc

#Accession paths
readarray accs_paths < ~/Wasney-Briscoe-2024/scripts/accession_paths.txt
accs_paths=(null ${accs_paths[@]}) # zero to one start index
acc_path=${accs_paths[$SGE_TASK_ID]}
echo $acc_path

OUTDIR=~/midas_output/${acc}

fastq1=~/fastqs/${acc_path}_1.fq.gz
fastq2=~/fastqs/${acc_path}_2.fq.gz
species_union=~/species_union/${acc}/species_union.txt

### MIDAS GENES
singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif run_midas.py genes $OUTDIR -1 $fastq1 -2 $fastq2 --extra_species_file $species_union --remove_temp

### MIDAS SNPs
#singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif run_midas.py snps $OUTDIR -1 $fastq1 -2 $fastq2 --extra_species_file $species_union
