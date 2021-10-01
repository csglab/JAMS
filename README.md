# JAMS: Joint Accessibility-Methylation-Sequence models
#### Quantitative modeling of sequence, DNA accessibility, and methylation preferences of transcription factors using ChIP-seq data.


JAMS creates models for inferring the effect of CpG methylation on TF binding in vivo. It creates quantitative models that connect the strength of the binding signal observed in ChIP-seq to the DNA accessibility of the binding site, regional methylation level, DNA sequence, and base-resolution cytosine methylation. 

JAMS has three main tasks (`DATA`, `GLM`, and `PREDICT`). Briefly, `DATA` extracts the information required to build the JAMS model from BAM and BIGWIG files and creates a directory with the files required by `GLM`. `GLM` builds a generalized linear model as described in the Method section of JAMS' [preprint](https://doi.org/10.1101/2021.08.27.457995). `PREDICT` predicts TF or background binding signal specified by a  BED file, it can either predict on the first position of each genomic region or look for the best motif match. 

#### **Requirements:** 

- Unix-compatible OS.  
- [R version 3.0.1](http://www.r-project.org/) or later.  
- R libraries: [data.table](https://www.rdocumentation.org/packages/data.table/versions/1.14.2), [ggplot2](https://www.rdocumentation.org/packages/ggplot2/versions/3.3.5), [ggseqlogo](https://github.com/omarwagih/ggseqlogo), [ggforce](https://www.rdocumentation.org/packages/ggforce/versions/0.3.3), [gridExtra](https://rdrr.io/cran/gridExtra/), [ggpubr](https://www.rdocumentation.org/packages/ggpubr/versions/0.4.0), [optparse](https://www.rdocumentation.org/packages/optparse/versions/1.6.6), [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html), [circlize](https://jokergoo.github.io/circlize/), [MASS](https://rdrr.io/cran/MASS/), and [patchwork](https://www.rdocumentation.org/packages/patchwork/versions/1.1.1).  
- [Python version 3.6.0](https://www.python.org/downloads/) or later.  
- [bwtool](https://github.com/CRG-Barcelona/bwtool) (requires [libbeato](https://github.com/CRG-Barcelona/libbeato), which has a compilation issue, a workaround can be found [here](https://github.com/CRG-Barcelona/libbeato/issues/6)).  
- [bedtools](https://bedtools.readthedocs.io/en/latest/index.html).

No installation is require for this version.  

#### **Usage:**  

Input files can be downloaded with `wget`:

```bash
cd ./data/CTCF_demo/01_input_files/
bash download_HEK293_CTCF_demo_data.sh
```

`01_JAMS_demo.sh` contains an example of JAMS for CTCF in HEK293 (tasks `DATA` and `PREDICT` in the run). This should create a `./data/CTCF_demo/03_output/JAMS_CTCF_HEK293_GSM2026781_small_flanking_20bps` folder.

For example:

```bash
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

#### **ARGUMENTS:**  

##### **Input files:** 

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

#### **Input and output directories:**

- `--data_dir <DATA_DIR>`

The path to the folder that contains the required files to build the JAMS model.

- `--output_dir <OUT_DIR>`

The path to the folder that will contain the output files.

#### **Other Arguments:**

- `--experiment <EXPERIMENT_ID>`

Experiment ID. ex: CTCF_rep2_HEK293.

- `--tag_summit_range <SUMMIT_RANGE>`

Range from the best motif match where the read tags will be extracted. Default= 400.

- `--flanking <FLANKING>`

The number of base pairs around the core motif that will be included in the model.


#### **Output:**

An example of JAMS' output is provided as `./data/CTCF_demo/03_output/CTCF_HEK293_GSM2026781_hg38_from_motif_RANGE_400_neg_binomial_CpG_av_Met_methylation_flank_20.tar.gz`.Next we describe the output files.

JAMS outputs the following files in the `<OUT_DIR>/<EXPERIMENT_ID>_flanking_<FLANKING>bps` folder:

- `<experimentName>_log.txt` – Text file containing the logs of the GLM task.
- `<experimentName>_mCpG_only_model_total_cv_scatterplot.pdf` –	 A PDF file showing the scatterplot of predicted vs. observed logarithm of tag density.
- `<experimentName>_dna_acc_vs_tags_heatmap.pdf` – PDF file comparing the DNA accessibility and the observed tags of each ChIP-seq peak. 
- `<experimentName>_DNA_accessibility_around_motif.pdf`	– Observed DNA accessibility on and around the best motif match.
- `<experimentName>_meth_vs_tags_heatmap_max_coverage.pdf`, and `<experimentName>_meth_vs_tags_heatmap_sensitive_coverage.pdf` – Heatmap with the coverage and percentage of WGBS reads on the best motif match. `max` plots the coverage without the scale being cap, `sensitive` being capped at 60 is more informative.
- `<experimentName>_methylation_correlation_with_binding.pdf` –	 PDF file summarizing the correlation of methylation at each position with the strength of TF binding. The statistics are presented for CpG's, non-CpG methylation events, and all methylation events combined. Note that since methylation events are locally correlated with each other, interpretation of these graphs can be difficult.
- `<experimentName>_methylation_stats.pdf` – PDF file summarizing the frequency of methylation at different positions of the sequences. The statistics are presented for CpG's, non-CpG methylation events, and all methylation events combined.
- `<experimentName>_sequence_vs_tags_heatmap.pdf` – DNA sequence of the motif and the observed tag counts from ChIP-seq.
- `<experimentName>_tags_vs_lenght.pdf`	– Scatterplot showing the tag counts from ChIP-seq and the length of the peak. 
- `<experimentName>_dna_coefficients.pdf` – DNA coefficients for the TF-specific and background parts of the JAMS model.
- `<experimentName>_mCpG_only_model_coefficients.txt` –	Text file containing the coefficients of a Negative binomial regression model that jointly considers DNA sequence, DNA accessibility  and methylation at CpG sites to explain ChIP-seq tag counts.`<experimentName>_mCpG_only_model_coefficients_with_FDR.txt` includes the FDR.
- `<experimentName>_mCpG_only_model_control_n_pulldown_motif_logo.pdf` – PDF file containing a graphical representation of the above model. The PDF essentially contains a motif logo (after scaling each position to have an average of zero), plus arrows that represent the effect of methylation: blue arrows represent methylation on the positive strand, and light brown arrows represent methylation on the reverse strand; upward arrows indicate positive and statistically significant effect on binding, and downward arrows represent negative and significant effect on binding. For both the TF-specific and background parts of the model.
- `<experimentName>_mCpG_only_model_control_motif_logo.pdf`, and `<experimentName>_mCpG_only_model_pulldown_motif_logo.pdf` – Same as above but the TF-specific and background motifs are in a different PDF file.
- `<experimentName>_mCpG_only_model_control.motif_scaled.txt`, and `<experimentName>_mCpG_only_model_pulldown.motif_scaled.txt` – Text file containing the exact values shown in the above PDF. Non-significant methylation effects are represented by NA's.
- `<experimentName>_mCpG_only_model_TF_binding.txt`, `<experimentName>_mCpG_only_model_TF_background.txt` – Coefficient files that can be used by the `PREDICT` task. 
- `<experimentName>_XX_train.tab`, `<experimentName>_c_train.tab`, and `<experimentName>_t_train.tab` – the XX, c and t matrices and vectors described in the Methods section of the [preprint](https://doi.org/10.1101/2021.08.27.457995).

<!-- 
- `<experimentName>_predictors_response.rda` –	
- `<experimentName>_dna_acc_by_bins.rda` –
- `<experimentName>_fit_CpG_only_model.rda` –
- `<experimentName>_fit_CpG_only_x.Met.rda` –
- `<experimentName>_methyl_chip_data.rda` –
- `<experimentName>_motif_pulldown.rda`	–
--> 

### **Task:** `DATA`

Only extracting data in and around the best motif match in ChIP-seq peaks. 

An example is included in the `02_JAMS_demo_format_data.sh` script:  

```bash
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


### **Task:** `GLM`

To build a JAMS model with the `GLM` task you need to specify the directory with the files created in the `DATA` task. 

For example:

```bash
./JAMS --task GLM \
       --experiment ${EXPERIMENT_ID} \
       --flanking ${FLANKING} \
       --data_dir ${DATA_DIR} \
       --output_dir ${OUT_DIR}
```

The script `03_JAMS_demo_GLM.sh` is an example of the `GLM` task. Before, running extract the input files with:

```bash
gunzip ./data/CTCF_demo/02_formatted_data/small_demo/*.gz
bash 03_JAMS_demo_GLM.sh
```

### **Task:** `Predict`

Predict TF or background signal on user specified genomic regions (BED format). It can either predict binding of the first position of each region (--motif_in_regions START) or look for the best motif match (--motif_in_regions SEARCH, if a PFM is provided). 260 high quality pretrain TF models are provided in `./data/JAMS_models/representative_JAMS_models.tar.gz`. `*_TF_binding.txt` and `*_background.txt` files can be used with the `--model_coeficients` argument.

For example:

```bash
./JAMS --task PREDICT \
       --motif_in_regions ${SEARCH_or_START} \
       --model_coeficients ${COEFFICIENTS} \
       --flanking 20 \
       --regions ${REGION} \
       --data_dir ${DATA_DIR} \
       --predicted_binding_prefix ${OUT_PREFIX} \
       --wgbs_met_data ${METH} --wgbs_unmet_data ${UNMETH} \
       --dna_acc_map ${DNA_ACC} --genome_fa ${GEN_FA} \
       --chr_sizes ${CHR_SIZES} --pfm ${PFM}
```

The script `04_JAMS_demo_predict.sh` is an example of the `PREDICT` task. 

#### **Citation:**
Hernandez-Corchado, A., & Najafabadi, H. S. (2021). A base-resolution panorama of the in vivo impact of cytosine methylation on transcription factor binding. BioRxiv. https://doi.org/10.1101/2021.08.27.457995








