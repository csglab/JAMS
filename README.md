# JAMS: Joint Accessibility-Methylation-Sequence models
#### Quantitative modeling of sequence, DNA accessibility, and methylation preferences of transcription factors using ChIP-seq data.

The JAMS command  has two main tasks (`DATA`, and `GLM`). `DATA` extracts the data required to build the JAMS model. `GLM` builds the model.  

### Requirements 
- Unix-compatible OS  
- R version 3.0.1 or later [Download](http://www.r-project.org/)  
- R libraries: data.table, ggplot2, ggseqlogo, ggforce, gridExtra, ggpubr, 
optparse, ComplexHeatmap, circlize, MASS, and patchwork  
- Python version 3.6.0 or later  
- [bwtool](https://github.com/CRG-Barcelona/bwtool)  
- [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)

No installation is require for this version.  

#### USAGE:  

```
./JAMS --task DATA,GLM \
       --experiment ${EXPERIMENT_ID} \
       --peaks ${PEAKS} \
       --wgbs_met_data ${METH} \
       --wgbs_unmet_data ${UNMETH} \
       --dna_acc_map ${DNA_ACC} \
       --genome_fa ${GEN_FA} \
       --pfm ${PFM} \
       --chr_sizes ${CHR_SIZES} \
       --mask_regions ${mask_regions} \
       --bam_control ${CONTROL_BAM} \
       --bam_pulldown ${PULLDOWN_BAM} \
       --tag_summit_range ${SUMMIT_RANGE} \
       --data_dir ${DATA_DIR}
       --flanking ${FLANKING} \
       --output_dir ${OUT_DIR}
```

To run an example of JAMS for CTCF in HEK293:  

`bash 01_JAMS_demo.sh` 

The input directory can be downloaded from Figshare with the command:

```
mkdir ./data/CTCF_demo/01_input_files
cd ./data/CTCF_demo/01_input_files
wget <Figshare link>
```

### ARGUMENTS:  
- `--peaks <PEAKS>`

Peaks from ChiP-seq experiment, output from MACS "*_peaks.xls".

- `--wgbs_met_data <METH>`

Methylated counts (BigWig).

- `--wgbs_unmet_data <UNMETH>`

Unmethylated counts (BigWig).

- `--dna_acc_map <DNA_ACC>`

Normalized DNA accessibility from ENCODE in bigwig file.

- `--bam_control <CONTROL_BAM>`

Indexed sorted bam file for the ChIP-seq control.

- `--bam_pulldown <PULLDOWN_BAM>`

Indexed and sorted bam file for ChIP-seq pulldown.

- `--genome_fa <GEN_FA>`

Reference sequence in the FASTA format (indexed).

- `--pfm <PFM>`

Position frequency matrix, Transfac format.

- `--chr_sizes <CHR_SIZES>`

Chromosome sizes.

- `--mask_regions <mask_regions>`

Genomic regions to be masked (repeats or/and blacklisted regions).

- `--experiment <EXPERIMENT_ID>`

Experiment ID. ex: CTCF_rep2_HEK293.

- `--tag_summit_range <SUMMIT_RANGE>`

Range from the best motif match where the read tags will be extracted. Default= 400.

- `--flanking <FLANKING>`

The number of base pairs around the core motif that will be included in the model.


- `--data_dir <DATA_DIR>`

The path to the folder that contains the required files to build the JAMS model.

- `--output_dir <OUT_DIR>`

The path to the folder that will contain the output files.




This should create a `./data/CTCF_demo/03_output/JAMS_CTCF_HEK293_GSM2026781_small_flanking_20bps` folder, with the output files described in the following sections.  

**Output:**

JAMS outputs the following files in the `<OUT_DIR>/<EXPERIMENT_ID>_flanking_<FLANKING>bps` folder:

- `<experimentName>_log.txt` – Text file containing messages produced by the script.

- `<experimentName>_XX_train.tab` –
- `<experimentName>_c_train.tab` –
- `<experimentName>_t_train.tab` –

- `<experimentName>_predictors_response.rda` –	
- `<experimentName>_dna_acc_by_bins.rda` –
- `<experimentName>_fit_CpG_only_model.rda` –
- `<experimentName>_fit_CpG_only_x.Met.rda` –
- `<experimentName>_methyl_chip_data.rda` –
- `<experimentName>_motif_pulldown.rda`	–

- `<experimentName>_dna_acc_vs_tags_heatmap.pdf` –
- `<experimentName>_dna_coefficients.pdf` –
- `<experimentName>_DNA_accessibility_around_motif.pdf`	–
- `<experimentName>_meth_vs_tags_heatmap_max_coverage.pdf` –
- `<experimentName>_meth_vs_tags_heatmap_sensitive_coverage.pdf` –
- `<experimentName>_methylation_correlation_with_binding.pdf` –	 PDF file summarizing the correlation of methylation at each position with the strength of TF binding. The statistics are presented for CpG's, non-CpG methylation events, and all methylation events combined. Note that since methylation events are locally correlated with each other, interpretation of these graphs can be difficult.
- `<experimentName>_methylation_stats.pdf` – PDF file summarizing the frequency of methylation at different positions of the sequences. The statistics are presented for CpG's, non-CpG methylation events, and all methylation events combined.
- `<experimentName>_sequence_vs_tags_heatmap.pdf` –
- `<experimentName>_tags_vs_lenght.pdf`	–

- `<experimentName>_mCpG_only_model_coefficients.txt` –	Text file containing the coefficients of a Poisson regression model that jointly considers DNA sequence and methylation at CpG sites to explain ChIP-seq tag counts.
- `<experimentName>_mCpG_only_model_coefficients_with_FDR.txt` –
- `<experimentName>_mCpG_only_model_total_cv_scatterplot.pdf` –	 A PDF file showing the scatterplot of predicted vs. observed logarithm of tag density.

- `<experimentName>_mCpG_only_model_control_n_pulldown_motif_logo.pdf` – PDF file containing a graphical representation of the above model. The PDF essentially contains a motif logo (after scaling each position to have an average of zero), plus arrows that represent the effect of methylation: blue arrows represent methylation on the positive strand, and light brown arrows represent methylation on the reverse strand; upward arrows indicate positive and statistically significant effect on binding, and downward arrows represent negative and significant effect on binding.
- `<experimentName>_mCpG_only_model_control.motif_scaled.txt` –	Text file containing the exact values shown in the above PDF. Non-significant methylation effects are represented by NA's.
- `<experimentName>_mCpG_only_model_control_motif_logo.pdf` –
- `<experimentName>_mCpG_only_model_pulldown.motif_scaled.txt` –
- `<experimentName>_mCpG_only_model_pulldown_motif_logo.pdf` –	


### Running JAMS' task: `DATA`

Only extracting data in and around the best motif match in ChIP-seq peaks. 

An example is included in the `02_JAMS_demo_format_data.sh` script:  

```
./JAMS --task DATA \
       --experiment ${EXPERIMENT_ID} \
       --peaks ${PEAKS} \
       --wgbs_met_data ${METH} \
       --wgbs_unmet_data ${UNMETH} \
       --dna_acc_map ${DNA_ACC} \
       --genome_fa ${GEN_FA} \
       --pfm ${PWM} \
       --chr_sizes ${CHR_SIZES} \
       --mask_regions ${mask_regions} \
       --bam_control ${CONTROL_BAM} \
       --bam_pulldown ${PULLDOWN_BAM} \
       --tag_summit_range ${SUMMIT_RANGE} \
       --data_dir ${DATA_DIR}
```

This creates the input files to build the JAMS model:

1. AffiMx output file: `<experimentName>_affimx.affinity.txt`
2. Aligned sequence file: `<experimentName>_aligned_sequences.tabulated.mx.txt`
3. Methylated read count file: `<experimentName>_aligned_sequences.fasta.methylreads`
4. Non-methylated read count file: `<experimentName>_aligned_sequences.fasta.nonmethylreads`
5. Normalized DNA accessibility file: `<experimentName>_aligned_sequences.fasta.accessibility`
6. ChIP-seq tag count file: `<experimentName>_summits.vRepeats.scores.txt`
7. PFM file for the motif that was used to align the sequences: `*_pfm.txt`


### Running JAMS' task: `GLM`

To build a JAMS model with the `GLM` task you need to specify the directory with the files created in the `DATA` task. 

For example:

```
./JAMS --task GLM \
       --experiment ${EXPERIMENT_ID} \
       --flanking {FLANKING} \
       --data_dir ${DATA_DIR} \
       --output_dir ${OUT_DIR}
```

The script `03_JAMS_demo_GLM.sh` is an example of the `GLM` task. Before, running extract the input files with:

```
gunzip ./data/CTCF_demo/02_formatted_data/small_demo/*.gz
bash 03_JAMS_demo_GLM.sh
```





