#!/bin/bash

which R

EXPERIMENT_ID=CTCF_HEK293_GSM2026781_smallest
DATA=./data/CTCF_demo/

DATA_DIR=${DATA}/02_formatted_data/smallest_demo
OUT_DIR=${DATA}/03_output

./JAMS --task GLM \
       --experiment ${EXPERIMENT_ID} \
       --flanking 20 \
       --data_dir ${DATA_DIR} \
       --output_dir ${OUT_DIR}
