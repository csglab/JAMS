#!/bin/bash

## Check software requirements
which bwtool bedtools python3 R

EXPERIMENT_ID=CTCF_HEK293_GSM2026781
DATA=./data/CTCF_demo/

FILE_DIR=${DATA}/01_input_files
DATA_DIR=${DATA}/02_formatted_data
OUT_DIR=${DATA}/03_output


PEAKS=${FILE_DIR}/GSM2026781_hg38_CTCF_input_Hughes_7_hg38_pval_1e-02_peaks.xls
METH=${FILE_DIR}/HEK293_hg38_MET.bw
UNMETH=${FILE_DIR}/HEK293_hg38_UNMET.bw
DNA_ACC=${FILE_DIR}/DNA_ACC_HEK293_HG38_ENCFF148BGE.bw
PFM=${FILE_DIR}/CTCF_HUMAN.pfm.txt
CHR_SIZES=${FILE_DIR}/hg38.chrom.sizes
MASK_REGIONS=${FILE_DIR}/blacklisted_and_repeats_sorted.bed
CONTROL_BAM=${FILE_DIR}/HEK293_input_Hughes_7.bam
PULLDOWN_BAM=${FILE_DIR}/GSM2026781_hg38.bam
SUMMIT_RANGE=400
GEN_FA=${FILE_DIR}/genome.fa


./JAMS --task DATA GLM \
       --experiment ${EXPERIMENT_ID} \
       --peaks ${PEAKS} \
       --wgbs_met_data ${METH} \
       --wgbs_unmet_data ${UNMETH} \
       --dna_acc_map ${DNA_ACC} \
       --genome_fa ${GEN_FA} \
       --pfm ${PFM} \
       --chr_sizes ${CHR_SIZES} \
       --mask_regions ${MASK_REGIONS} \
       --bam_control ${CONTROL_BAM} \
       --bam_pulldown ${PULLDOWN_BAM} \
       --tag_summit_range ${SUMMIT_RANGE} \
       --flanking 20 \
       --output_dir ${OUT_DIR}
