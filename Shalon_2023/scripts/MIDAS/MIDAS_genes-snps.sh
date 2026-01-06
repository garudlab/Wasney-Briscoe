#!/bin/bash
#$ -N MIDAS_genes_snps
#$ -e /u/project/ngarud/michaelw/PaulAllen/Shalon_2023/scripts/errors/
#$ -o /u/project/ngarud/michaelw/PaulAllen/Shalon_2023/scripts/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=20G
#$ -l time=72:00:00
#$ -l highp
#$ -tc 60 # Throttle to max 60 tasks at a time
#$ -t 1-302:1

#Load modules
source /u/local/Modules/default/init/modules.sh
module load anaconda3
module load singularity

#Load python 2.7 conda environment
conda activate python27_env


export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS_mod
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS_mod/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2

#Accession
readarray accs < /u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/metadata/SRA_list.txt
accs=(null ${accs[@]}) # zero to one start index
acc=${accs[$SGE_TASK_ID]}

OUTDIR=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/midas_output/${acc}

fastq1=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/fastqs/${acc}/${acc}_1.fastq.bz2
fastq2=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/fastqs/${acc}/${acc}_2.fastq.bz2
species_union=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/species_union/${acc}/species_union.txt

#run_midas.py genes $OUTDIR -1 $fastq1 -2 $fastq2 --extra_species_file $species_union --remove_temp


### MIDAS GENES
singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif run_midas.py genes $OUTDIR -1 $fastq1 -2 $fastq2 --extra_species_file $species_union --remove_temp

### MIDAS SNPs
singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif run_midas.py snps $OUTDIR -1 $fastq1 -2 $fastq2 --extra_species_file $species_union
