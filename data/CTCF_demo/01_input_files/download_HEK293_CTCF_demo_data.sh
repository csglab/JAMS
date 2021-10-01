#!/bin/bash

## Macs peak and bed files
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/GSM2026781_CTCF_rep2.macs.out.txt
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/GSM2026781_hg38_CTCF_input_Hughes_7_hg38_pval_1e-02_peaks.xls
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/GSM2026781_hg38_CTCF_input_Hughes_7_hg38_pval_1e-02_peaks_10k.xls
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/GSM2026781_hg38_CTCF_input_Hughes_7_hg38_pval_1e-02_peaks_test_predict.bed
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/blacklisted_and_repeats_sorted.bed
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/random.bed
## CTCF protein sequence and PFM
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/CTCF_HUMAN.fasta
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/CTCF_HUMAN.pfm.txt
## CTCF ChIP-seq pulldown
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/GSM2026781_hg38.bam
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/GSM2026781_hg38.bam.bai
## ChIP-seq input
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/HEK293_input_Hughes_7.bam
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/HEK293_input_Hughes_7.bam.bai
## Methylated and Unmethylated count bigwig files (WGBS) - 4.5 Gb and 8.2 Gb
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/HEK293_hg38_MET.bw
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/HEK293_hg38_UNMET.bw
## Normalized DNA accessibility from ENCODE (DNase-seq)
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/DNA_ACC_HEK293_HG38_ENCFF148BGE.bw
## UCSC hg38 fasta
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/genome.fa
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/genome.fa.fai
wget https://data.cyverse.org/dav-anon/iplant/home/ahcorcha/JAMS_input/hg38.chrom.sizes
