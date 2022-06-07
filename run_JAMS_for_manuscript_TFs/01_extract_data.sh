#!/bin/bash

## Check software requirements
which bwtool bedtools python3

EXPERIMENT_ID=TF_cell_type_file_ID

OUT_DATA_DIR=./data/JAMS_input

## MACS peaks are found in Zenodo
### Zenodo DOI: 10.5281/zenodo.5573260
### File: macs_peaks.tar.gz
PEAKS=peak_example_hg38_pval_1e-02_peaks.xls

##  Methylated and Unmethylated count bigwig files are
### created as described here:
###  https://github.com/ahcorcha/bismark-examples
METH=${FILE_DIR}/${CELL_TYPE}_MET.bw
UNMETH=${FILE_DIR}/${CELL_TYPE}_UNMET.bw

## ENCODE DNase-seq read-depth normalized signal
### download with:
### wget https://www.encodeproject.org/files/ENCFF148BGE/@@download/ENCFF148BGE.bigWig
# HEK293 ENCFF148BGE
# GM12878 ENCFF743ULW
# H1 ENCFF131HMO
# HEK293 ENCFF148BGE
# HeLa-S3 ENCFF256QQH
# HepG2 ENCFF577SOF
# K562 ENCFF413AHU

DNA_ACC=${FILE_DIR}/DNA_ACC_${CELL_TYPE}.bw

## Seed motif from CIS-BP (i.e. M04716_2.00.txt) or RCADE2
### CIS-BP motifs can be found here:
### RCADE2 motifs can be found on Zenodo DOI: 10.5281/zenodo.5573260
### File: 
PWM=${SEED_MOTIF}.pfm.txt


CHR_SIZES=hg38.chrom.sizes


## Regions to be excludede from the analysis.
### For example, ENCODE blacklisted motifs or repeatmasker regions
MASK_REGIONS=blacklisted_and_repeats_sorted.bed


## TF ChIP-seq and input bam files
CONTROL_BAM=example_control.bam
PULLDOWN_BAM=example_input.bam


## Genome fasta dile
GEN_FA=genome.fa


JAMS --task DATA \
     --experiment ${EXPERIMENT_ID} \
     --region ${PEAKS} \
     --wgbs_met_data ${METH} \
     --wgbs_unmet_data ${UNMETH} \
     --dna_acc_map ${DNA_ACC} \
     --genome_fa ${GEN_FA} \
     --pfm ${PWM} \
     --chr_sizes ${CHR_SIZES} \
     --mask_regions ${MASK_REGIONS} \
     --bam_control ${CONTROL_BAM} \
     --bam_pulldown ${PULLDOWN_BAM} \
     --tag_summit_range 400 \
     --data_dir ${OUT_DATA_DIR}
