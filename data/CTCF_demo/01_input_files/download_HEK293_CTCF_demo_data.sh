#!/bin/bash

## Macs peak and bed files
wget https://usegalaxy.org/api/datasets/f9cad7b01a47213510ac7061c52ad932/display?to_ext=bed --output-document=blacklisted_and_repeats_sorted.bed
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135454b1d37676afb41/display?to_ext=tabular --output-document=GSM2026781_hg38_CTCF_input_Hughes_7_hg38_pval_1e-02_peaks_test_predict.bed
wget https://usegalaxy.org/api/datasets/f9cad7b01a4721358fcf84dd05153fb1/display?to_ext=tabular --output-document=GSM2026781_hg38_CTCF_input_Hughes_7_hg38_pval_1e-02_peaks_10k.xls
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135fa63bb50d2370265/display?to_ext=tabular --output-document=GSM2026781_hg38_CTCF_input_Hughes_7_hg38_pval_1e-02_peaks.xls
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135afb632b0a872f359/display?to_ext=tabular --output-document=GSM2026781_CTCF_rep2.macs.out.txt
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135b9f8b268efff7465/display?to_ext=bed --output-document=random.bed

## CTCF protein sequence and PFM
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135c26102a78d40c02d/display?to_ext=fasta --output-document=CTCF_HUMAN.fasta
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135a24fc3865f5cb2ce/display?to_ext=txt --output-document=CTCF_HUMAN.pfm.txt


## CTCF ChIP-seq pulldown
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135aef1b502acd85799/display?to_ext=bam --output-document=GSM2026781_hg38.bam
wget https://usegalaxy.org/api/datasets/f9cad7b01a4721356c8233d8f4b71550/display?to_ext=bai --output-document=GSM2026781_hg38.bam.bai

## ChIP-seq input
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135535e55c55d2c3d54/display?to_ext=bam --output-document=HEK293_input_Hughes_7.bam
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135ce84c96e2b15fa0c/display?to_ext=bai --output-document=HEK293_input_Hughes_7.bam.bai

## Normalized DNA accessibility from ENCODE (DNase-seq)
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135b6adc9de003de73d/display?to_ext=bigwig --output-document=DNA_ACC_HEK293_HG38_ENCFF148BGE.bw

## UCSC hg38 fasta
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135467c41f0d3462f2e/display?to_ext=gg --output-document=hg38.chrom.sizes
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135242ebe75ef6905b1/display?to_ext=fasta --output-document=genome.fa
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135db4e21b81f555709/display?to_ext=bed --output-document=genome.fa.fai

## Methylated and Unmethylated count bigwig files (WGBS) - 4.5 Gb and 8.2 Gb
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135b5f71ccde458fc27/display?to_ext=bigwig --output-document=HEK293_hg38_MET.bw
wget https://usegalaxy.org/api/datasets/f9cad7b01a472135490cd7c33b6e7cc2/display?to_ext=bigwig --output-document=HEK293_hg38_UNMET.bw

