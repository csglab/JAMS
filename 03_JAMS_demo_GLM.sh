#!/bin/bash

which R

EXPERIMENT_ID=CTCF_HEK293_GSM2026781_small
DATA=~/repos/tools/vanilla_JAMS/JAMS/data/CTCF_demo

DATA_DIR=${DATA}/02_formatted_data/small_demo
OUT_DIR=${DATA}/03_output

set -o xtrace
JAMS --task BUILD \
     --experiment ${EXPERIMENT_ID} \
     --flanking 20 \
     --data_dir ${DATA_DIR} \
     --output_dir ${OUT_DIR}
