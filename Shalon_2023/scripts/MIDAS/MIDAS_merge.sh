#!/bin/bash
#$ -N MIDAS_merge
#$ -e /u/project/ngarud/michaelw/PaulAllen/Shalon_2023/scripts/errors/
#$ -o /u/project/ngarud/michaelw/PaulAllen/Shalon_2023/scripts/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=20G
#$ -l time=168:00:00
#$ -l highp

#Load modules
source /u/local/Modules/default/init/modules.sh
module load anaconda3
module load singularity

#Load python 2.7 conda environment
conda activate python27_env

#Sources
export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2

#Making the merged data directories
if [ -d /u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/merged_data/ ]; then
	echo "merged_data/ directory already exists"
else
	echo "Creating merged_data/"
	mkdir /u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/merged_data/
fi

if [ -d /u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/merged_data/species/ ]; then
	echo "merged_data/species/ directory already exists"
else
	echo "Creating merged_data/species/"
	mkdir /u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/merged_data/species/
fi

if [ -d /u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/merged_data/genes/ ]; then
	echo "merged_data/genes/ directory already exists"
else
	echo "Creating merged_data/genes/"
	mkdir /u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/merged_data/genes/
fi

if [ -d /u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/merged_data/snps/ ]; then
	echo "merged_data/snps/ directory already exists"
else
	echo "Creating merged_data/snps/"
	mkdir /u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/merged_data/snps/
fi


INDIR=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/midas_output/
OUTDIR=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/merged_data/


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
