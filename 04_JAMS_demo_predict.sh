#!/bin/bash
which R
DATA=./data/CTCF_demo/

FILE_DIR=${DATA}/01_input_files
METH=${FILE_DIR}/HEK293_hg38_MET.bw
UNMETH=${FILE_DIR}/HEK293_hg38_UNMET.bw
DNA_ACC=${FILE_DIR}/DNA_ACC_HEK293_HG38_ENCFF148BGE.bw
GEN_FA=${FILE_DIR}/genome.fa
CHR_SIZES=${FILE_DIR}/hg38.chrom.sizes
PFM=${FILE_DIR}/CTCF_HUMAN.pfm.txt

REGION=data/CTCF_demo/01_input_files/GSM2026781_hg38_CTCF_input_Hughes_7_hg38_pval_1e-02_peaks_test_predict.bed
COEFFICIENTS=${DATA}/04_predict/CTCF_HEK293_GSM2026781_small_mCpG_only_model_coefficients_with_FDR_TF_binding.txt

OUT_DIR=${DATA}/04_predict
mkdir -p ${OUT_DIR}


EXPERIMENT_ID=CTCF_HEK293_GSM2026781_small
JAMS --task PREDICT \
     --motif_in_regions SEARCH \
     --model_coeficients ${COEFFICIENTS} \
     --flanking 20 \
     --regions ${REGION} \
     --output_dir ${OUT_DIR} \
     --experiment ${EXPERIMENT_ID} \
     --wgbs_met_data ${METH} --wgbs_unmet_data ${UNMETH} \
     --dna_acc_map ${DNA_ACC} --genome_fa ${GEN_FA} \
     --chr_sizes ${CHR_SIZES} --pfm ${PFM}

JAMS --task PREDICT \
     --motif_in_regions START \
     --model_coeficients ${COEFFICIENTS} \
     --flanking 20 \
     --regions ${REGION} \
     --output_dir ${OUT_DIR} \
     --experiment ${EXPERIMENT_ID} \
     --wgbs_met_data ${METH} --wgbs_unmet_data ${UNMETH} \
     --dna_acc_map ${DNA_ACC} --genome_fa ${GEN_FA} \
     --chr_sizes ${CHR_SIZES} --pfm ${PFM}


REGION=data/CTCF_demo/01_input_files/random.bed 

EXPERIMENT_ID=CTCF_HEK293_GSM2026781_small_random_regions
JAMS --task PREDICT \
     --motif_in_regions START \
     --model_coeficients ${COEFFICIENTS} \
     --flanking 20 \
     --regions ${REGION} \
     --output_dir ${OUT_DIR} \
     --experiment ${EXPERIMENT_ID} \
     --wgbs_met_data ${METH} --wgbs_unmet_data ${UNMETH} \
     --dna_acc_map ${DNA_ACC} --genome_fa ${GEN_FA} \
     --chr_sizes ${CHR_SIZES} --pfm ${PFM}

JAMS --task PREDICT \
     --motif_in_regions SEARCH \
     --model_coeficients ${COEFFICIENTS} \
     --flanking 20 \
     --regions ${REGION} \
     --output_dir ${OUT_DIR} \
     --experiment ${EXPERIMENT_ID} \
     --wgbs_met_data ${METH} --wgbs_unmet_data ${UNMETH} \
     --dna_acc_map ${DNA_ACC} --genome_fa ${GEN_FA} \
     --chr_sizes ${CHR_SIZES} --pfm ${PFM}
