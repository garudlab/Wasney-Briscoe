#!/bin/bash

source /u/local/Modules/default/init/modules.sh
module load anaconda3
source activate python_env

python optimized_site_downsampling.py --strain $1 --metric median
