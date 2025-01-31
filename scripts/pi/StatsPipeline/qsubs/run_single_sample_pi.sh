#!/bin/bash

source /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate python_env

while read -r arg_1 arg_2 arg_3 arg_4 arg_5 arg_6 arg_7; do
    python ../pop_gen_calculator/pop_gen_calculator.py --sample_id $arg_1 --sample_id2 "" --end_index $arg_2 --strain $arg_3 --study $arg_4 --calc_type $arg_5 --statistic $arg_6 --cov_min $arg_7

done < inputs/data_$SGE_TASK_ID.in






