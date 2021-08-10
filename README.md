# JAMS

### Quantitative modeling of sequence, DNA accessibility, and methylation preferences of transcription factors using ChIP-seq data

## Requirements 
- Unix-compatible OS
- R version 3.0.1 or later [Download](http://www.r-project.org/) 
- R “data.table” library
- R “ggplot2” library
- R “ggseqlogo” library

## Installation 
- No installation is require for this version
- To test the script, execute this command: 


RCADE2=/home/ahcorcha/tools/RCADE2
Methyl_ChIP=/home/ahcorcha/tools/MethylChIP
## bedtools, bwtools, meme-suite, RCADE2 and MethylChIP. Python3 modules: numpy (for RCADE2),



>
    Rscript glm_poisson.methylation.R CTCF_rep2 20 ./example/affimx.100bp ./methyl_glm.out
This should create a `./methyl_glm.out/CTCF_rep2` folder, with the output files described in the following sections.

## Running MethylChIP
### Required files
The current version is designed to work with the file structure used by the CSG lab for analysis of C2H2-ZF proteins. MethylChIP requires all the following files to be present in the same folder:

1. AffiMx output file `<experimentName>.<UniProtID>.affimx.affinity.txt`
2. Aligned sequence file `<experimentName>.<UniProtID>.aligned_sequences.tabulated.mx.txt`
3. Methylated read count file `<experimentName>.<UniProtID>.aligned_sequences.fasta.methylreads`
4. Non-methylated read count file `<experimentName>.<UniProtID>.aligned_sequences.fasta.nonmethylreads`
5. ChIP-seq tag count file `<experimentName>.summits.vRepeats.scores.txt`
6. PFM file for the motif that was used to align the sequences `<UniProtID>.pfm.txt`

###Usage
Use the glm_poisson.methylation.R script to run MethylChIP:

    Rscript glm_poisson.methylation.R <experimentName> <flankingSequenceLength> <inputFolder> <outputFolder>

**Inputs:**

`<experimentName>`: The ChIP-seq experiment label, e.g. `CTCF_rep2`.

`<flankingSequenceLength>`: The number of base pairs around the core motif that will be included in the model.

`<inputFolder>`: The path to the folder that contains the required files.

`<outputFolder>`: The path to the folder that will contain the output files. If it does not exist, this folder will be created automatically.

## Output
MethylChIP outputs the following files in the `<outputFolder>/<experimentName>` folder:

- `log.txt` – Text file containing messages produced by the script

- `methylation_stats.pdf` – PDF file summarizing the frequency of methylation at different positions of the sequences. The statistics are presented for CpG's, non-CpG methylation events, and all methylation events combined.

- `methylation_correlation_with_binding` – PDF file summarizing the correlation of methylation at each position with the strength of TF binding. The statistics are presented for CpG's, non-CpG methylation events, and all methylation events combined. Note that since methylation events are locally correlated with each other, interpretation of these graphs can be difficult.

- `mCpG_only_model.coefficients.txt` – Text file containing the coefficients of a Poisson regression model that jointly considers DNA sequence and methylation at CpG sites to explain ChIP-seq tag counts.

- `mCpG_only_model.motif_logo.pdf` – PDF file containing a graphical representation of the above model. The PDF essentially contains a motif logo (after scaling each position to have an average of zero), plus arrows that represent the effect of methylation: blue arrows represent methylation on the positive strand, and light brown arrows represent methylation on the reverse strand; upward arrows indicate positive and statistically significant effect on binding, and downward arrows represent negative and significant effect on binding.

- `mCpG_only_model.motif_scaled.txt` – Text file containing the exact values shown in the above PDF. Non-significant methylation effects are represented by NA's.

- `mCpG_only_model.scatterplot.jpg` – A JPEG file showing the scatterplot of predicted vs. observed logarithm of tag density.

- `All_mC_only_model.coefficients.txt` – Text file containing the coefficients of a Poisson regression model that jointly considers DNA sequence and all methylation events (including CpG and non-CpG sites) to explain ChIP-seq tag counts.

- `All_mC_only_model.motif_logo.pdf` – PDF file containing a graphical representation of the above model.

- `All_mC_only_model.motif_scaled.txt` – Text file containing the exact values shown in the above PDF.

- `All_mC_only_model.scatterplot.jpg` – A JPEG file showing the scatterplot of predicted vs. observed logarithm of tag density.

- `Merged_model.motif_logo.pdf` – PDF file containing a graphical representation of the methylation associations that are statistically significant (and in the same direction) in both the CpG and all-C models.

- `Merged_model.motif_scaled.txt` – Text file containing the exact values shown in the above PDF.
