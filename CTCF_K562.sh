#!/bin/bash

#SBATCH --account=def-hsn
#SBATCH --time=00:10:00
#SBATCH --mincpus=1
#SBATCH --job-name=binding_affinity
#SBATCH --mail-user=ahcorcha@gmail.com
#SBATCH --output=logs/%x_job_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mem=5000M

date
## WBGS data
wgbs_met=./data/methyl_maps/hg38/K562/ENCSR765JPC_K562_methyl_map_rep_1_methylated.bw
wgbs_unmet=./data/methyl_maps/hg38/K562/ENCSR765JPC_K562_methyl_map_rep_1_unmethylated.bw


peaks=/home/ahcorcha/projects/rrg-hsn/ahcorcha/P3_tf_methyl/prepare_data/peaks/hg19/CTCF_ENCFF000PYD_K562_rep1_new.bed

pfm=./data/demo/CTCF/in/CTCF_HUMAN.pfm.txt
uni_id=CTCF_HUMAN
prot_name=CTCF
exp_id=ENCFF000PYD_CTCF_K562_rep1
range=100 # For pre-methyl-chip
out=./out2/CTCF_K562

mkdir -p ${out}
hg_38=./data/hg38/hg38.fa
prot_seq=./data/demo/CTCF/in/CTCF_HUMAN.fasta
PFM=./data/demo/CTCF/in/CTCF_HUMAN.pfm.txt
repeats=./data/hg38/rmsk.hg38.sorted.bed
chr_sizes=./data/hg38/hg38.chrom.sizes
dna_acc=./data/dna_accessibility/K562/ENCFF936BDN_K562_dna_acc_rep2.bw

# dna_acc=./data/dna_accessibility/K562/ENCFF413AHU_K562_dna_acc_rep1.bw

# tasks:  Methyl-ChIP
./binding_affinity --peaks ${peaks} \
		   --wgbs_met_data ${wgbs_met} \
		   --wgbs_unmet_data ${wgbs_unmet} \
		   --uniprot_id ${uni_id} \
		   --prot_name ${prot_name} \
		   --experiment_id ${exp_id} \
		   --genome ${hg_38} \
		   --range ${range} \
		   --prot_fam C2H2 \
		   --c2h2_seq ${prot_seq} \
		   --out_dir ${out} \
		   --peaks_format MACS2 \
 		   --pfm ${PFM} \
		   --task Methyl-ChIP \
		   --dna_acc_map ${dna_acc} \
		   --rc_top_peaks 500
date
