#!/bin/bash
#SBATCH --account=def-hsn
#SBATCH --time=00:30:00
#SBATCH --mincpus=1
#SBATCH --job-name=binding_affinity
#SBATCH --mail-user=ahcorcha@gmail.com
#SBATCH --output=logs/%x_job_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mem=5000M
## This will be an array job, reading a tab file with the all the experiments.
date


source activate methyl_chip

# DEMO: CTCF
## WBGS data
wgbs_met=./data/methyl_maps/hg19/hek293/GSM1254259_HEK293-CT.ME.cout.bw
wgbs_unmet=./data/methyl_maps/hg19/hek293/GSM1254259_HEK293-CT.NM.cout.bw
## Chip-seq data from macs2
peaks=./data/demo/CTCF/in/GSM2026781_CTCF_rep2.macs.out.txt
pfm=./data/demo/CTCF/in/CTCF_HUMAN.pfm.txt
uni_id=CTCF_HUMAN
prot_name=CTCF
exp_id=CTCF_Hek293_HS_rep2
range=100 # For pre-methyl-chip
out=./out/CTCF_Hek293

# $PS=/home/ahcorcha/projects/rrg-hsn/ahcorcha
hg_19=$PS/resources/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
prot_seq=./data/demo/CTCF/in/CTCF_HUMAN.fasta
PFM=./data/demo/CTCF/in/CTCF_HUMAN.pfm.txt
dna_acc=./data/dna_accessibility/HEK293/DNA_ACC_HEK293_HG19_ENCFF716SFD.bw


## Defaults are included in data/hg19, you can specify your own
repeats=./data/hg19/rmsk.hg19.sorted.bed
chr_sizes=./data/hg19/hg19.chrom.sizes

# tasks: RCADE2 Methyl-ChIP # For now, it's just running pre-Methyl
# 		   --pfm ${PFM} \

./binding_affinity --peaks ${peaks} \
		   --wgbs_met_data ${wgbs_met} \
		   --wgbs_unmet_data ${wgbs_unmet} \
		   --uniprot_id ${uni_id} \
		   --prot_name ${prot_name} \
		   --experiment_id ${exp_id} \
		   --genome ${hg_19} \
		   --range ${range} \
		   --prot_fam C2H2 \
		   --c2h2_seq ${prot_seq} \
		   --out_dir ${out} \
		   --peaks_format MACS2 \
		   --dna_acc_map ${dna_acc} \
		   --task RCADE2 \
		   --rc_top_peaks 1000 \
		   --chr_sizes ${chr_sizes} \
		   --repeats ${repeats}
date

# /home/ahcorcha/projects/rrg-hsn/ahcorcha/P3_tf_methyl/binding_affinity/data/demo/CTCF/out/methyl_chip/CTCF_Hek293_HS_rep2
