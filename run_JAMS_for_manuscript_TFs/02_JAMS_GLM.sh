#!/bin/bash
#SBATCH --account=def-hsn
#SBATCH --mincpus=1
#SBATCH --array=1-260
#SBATCH --time=01:30:00
#SBATCH --mem=10000M
#SBATCH --job-name=run_JAMS
#SBATCH --output=%x_job_%j.out

line_number=1 # Which experiment
ref=./ref/experiment_IDs.txt

line=$( sed -n "${line_number}p" ${ref} )
EXPERIMENT_ID=$(echo $line | awk '{print $1}')

IN_DATA_DIR=./data/JAMS_input/${EXPERIMENT_ID}
OUT_DIR=./data/JAMS_output
mkdir -p ${OUT_DIR}

../JAMS --task GLM \
	--experiment ${EXPERIMENT_ID} \
	--flanking 20 \
	--data_dir ${IN_DATA_DIR} \
	--output_dir ${OUT_DIR}
