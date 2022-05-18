suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggseqlogo))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reticulate))
options( error = traceback, nwarnings = 10000 )

########################################################   IN and load data ####
option_list = list(
  make_option(c("-e", "--experiment"), type="character",
              default="test",
              help="Experiment ID", metavar="character"),

  make_option(c("-f", "--flanking"), type="integer", default=20,
              help="length of flanking sequence around the motif", 
              metavar="character"),
  
  make_option(c("-l", "--pfm_length"), type="integer", default=15,
              help="", metavar="character"),
  
  make_option(c("-d", "--input_dir"), type="character", metavar="character",
              default="./data/CTCF_demo/02_formatted_data/smallest_demo",
              help="Input directory with PFM, methylation counts etc ..."),

  make_option(c("-i", "--iterations"), type="character", metavar="character",
              default=10,
              help="Input directory with PFM, methylation counts etc ..."),  

  make_option(c("-p", "--path_to_JAMS"), type="character", metavar="character",
              default="/home/ahcorcha/repos/tools/JAMS/",
              help="Input directory with PFM, methylation counts etc ..."),    
    
  make_option(c("-o", "--output_dir"), type="character",
              default="./data/CTCF_demo/05_motif_discovery",
              help="", metavar="character"),
  
  make_option(c("-m", "--exclude_meth"), type="logical",
              action = "store_true", default = "FALSE",
              help="", metavar="character")
    
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)


source( paste0( opt$path_to_JAMS, "/src/Methyl_ChIP_ftns.R") )
## Source_python is required for source de_novo_discovery_ftns.R
source_python(paste0( opt$path_to_JAMS, "/src/motif_discovery.py" ) )
source( paste0( opt$path_to_JAMS, "/src/de_novo_discovery_ftns.R" ) )

outdir <- opt$output_dir
pfm_length <- opt$pfm_length
flanking <- opt$flanking
input_root <- opt$input_dir
iterations <- opt$iterations
experiment <- paste0( opt$experiment, "_motif_length_", pfm_length, 
                      "_flanking_", flanking )

# Used for testing with Rstudio, normally commented
# setwd(opt$path_to_JAMS) # ahcorcha

outdir <- paste0( outdir, "/de_novo_motif_", experiment )
dir.create( outdir )
prefix <- paste0( outdir, "/", experiment )

sink( paste0( prefix, "_log.txt" ) )
cat(paste0( experiment, "\n"))
cat(paste0( "Start wall time: ", Sys.time(), "\n"))

if ( opt$exclude_meth ) {
  cat( paste0( "Training TF models without intra-motif methylation, ",
               "these coefficients will be added to *coefficients_with_FDR.txt ", 
               "for conformity with other functions\n"))
}
  
cat("Loading data ...\n")
dat_all <- load_dat( input_root, pfm_length = pfm_length )


###########################################################   Pre iteration ####
## The starting position here makes the center of the motif is at the center of the peaks

## Start at the peak's center, intuitive
# start_pos <- rep_len(x = 101 , # Start at peak's middle
#                      length.out = nrow(dat_all$x.Met.all)) # Number of ChIP-seq peaks

## Start at random positions, for testing
set.seed(5)
start_pos <- floor( runif( nrow( dat_all$x.A.all ),
                           min=flanking+1,
                           max= ncol(dat_all$x.A.all) - pfm_length - flanking  ) )

start_pos_list <-list( start_pos )

## The furthest right we can go, position at the start of the left flank
possible_position <- ncol( dat_all$x.C.all ) - 2 * flanking - pfm_length
## Pre compute position specific predictors 
### Create list of length upper_limit_pos, each entry is the dat_all for the 
### correspondent position.
cat("Pre-calc. predictors per position ...\n")
predictors_list <- pre_calc_by_pos_dat( this_dat_all = dat_all, 
                                        possible_position = possible_position, 
                                        flanking = flanking, 
                                        pfm_length = pfm_length )
## Random (but constant from iteration to iteration) peaks for motif heatmap 
rnd_num <- sort( sample.int( nrow( dat_all$x.A.all ), 7500 ) )

###############################################################   Iteration ####
cat(paste0( "Start iterations wall time: ", Sys.time(), "\n"))
# iterations <- 20
for (i in 1:as.integer(iterations)) {
  ### Starting iteration
  # i <- 1
  cat( paste0( "Iteration: ", i, "\n" ) )
  prefix_iteration <- paste0( prefix, "_iteration_", format_iteration(i) )

  #############################################################   Train GLM ####
  this_glm <- train_GLM_at_shifted_pos( flanking = flanking, 
                                        pfm_length = pfm_length, 
                                        dat_all = dat_all,
                                        start_pos = start_pos,
                                        exclude_meth = opt$exclude_meth
                                        )
  
  ############### Evaluate every position within +/- 200 bps of peak center ####
  pdwn_coeffs <- as.data.frame( coefficients( summary( this_glm ) ) )
  
  ## Get the pulldown coefficients
  pdwn_coeffs <- pdwn_coeffs[ grepl(pattern = ":t$",
                              x = rownames(pdwn_coeffs)),]
  
  ## Get the name of the variables included in th model
  X_var_names <- gsub( ":t", "", rownames( pdwn_coeffs ) )
  
  ## Make sure those variables are in the correct order for the matrix*vector
  ## Mult. the PULLDOWN coeffs vector and the variable matrix for each position
  c_pdwn_predicted <- sapply( X = predictors_list, 
                              FUN = eval_coeffs, 
                              pdn_coeff = pdwn_coeffs$Estimate, 
                              X_names = X_var_names )
    
  rownames(c_pdwn_predicted) <- rownames(predictors_list[[1]])
  c_pdwn_predicted <- as.data.frame( c_pdwn_predicted )
  
  ## Per row, get column name of max value (max TF signal).
  new_start_pos <- colnames(c_pdwn_predicted)[max.col(c_pdwn_predicted, 
                                                      ties.method = "first")]

  new_start_pos <- as.numeric( gsub( "V", "", new_start_pos ) )
  new_start_pos <- new_start_pos + flanking ## ahcorcha
  
  ## Change start_pos to the ones with max TF signal
  len <- length(start_pos_list)
  start_pos_list[[len+1]] <- new_start_pos
  
  # compare new_start_pos with start_pos
  mean_abs_pos_change <- mean( abs( start_pos - new_start_pos ) )
  cat( paste0("   Mean absolute position change: ", round(mean_abs_pos_change, 4), "\n") )
  
  median_abs_pos_change <- median( abs( start_pos - new_start_pos ) )
  cat( paste0(" Median absolute position change: ", round(median_abs_pos_change, 4), "\n") )
  
  abs_change_pos <- ( start_pos - new_start_pos )
  abs_change_pos_df <- as.data.frame( abs_change_pos )
  
  phist <- ggplot( abs_change_pos_df, aes( x = abs_change_pos ) ) + 
                   geom_histogram( aes(y=..count..), binwidth = 2, 
                                   colour="black", fill="white") +
                   xlim( -( ncol(dat_all$x.C.all)-2*flanking), 
                          ( ncol(dat_all$x.C.all)-2*flanking) ) +
                   labs(x = "Position change", y = "Density") +
                   theme_light()
  
  start_pos <- new_start_pos 
  
  #############################   Visualize motif start pos over iterations ####
  ###### Save run's information: write coefficients / draw logo and DNA coeffs
  dna_acc_plot_name <- paste0( prefix_iteration, "_dna_coefficients.pdf" )

  p_dna_coeffs <- plot_dna_acc_coefficients( this_glm )
  motif_coefs <- write.sequence.model.av.met( seq_fit = this_glm, opt$exclude_meth )

  write.table( x =  motif_coefs[[1]], quote = FALSE, sep = "\t",
               col.names = TRUE, row.names = TRUE,
               file = paste0( prefix_iteration, "_coefficients_with_FDR.txt") )

  p_motif_coefs <- ( motif_coefs[[2]] / p_dna_coeffs )
  
  p_motif_coefs <- p_motif_coefs +
             plot_annotation( title = paste0(experiment, ", iteration: ", i ) ) +
             plot_layout( heights = c(2, 0.75) )

  # ggsave( filename = paste0( prefix_iteration, "_motifs_and_dna_coeffs.pdf" ),
  #         plot = p_motif_coefs, height = 9, width = 9 )
  
  sample_start_pos <- start_pos[rnd_num]
  
  motif_ht <- motif_pos_heatmap( sample_start_pos, n_cols = ncol(dat_all$x.C.all), pfm_length, i )
  
  p_motif_ht <- grid.grabExpr( draw(motif_ht, heatmap_legend_side = "bottom") )
  
  layout <- "AADD
             AADD
             BBDD
             CCDD"

  
  this.tittle <- paste0( experiment,
                         ", flanking = ", flanking,
                         ", motif length = ", pfm_length,
                         ", iteration = ", i,
                         ", n peaks = ", nrow(dat_all$x.C.all))
  
  if(opt$exclude_meth) { 
    this.tittle <- paste0(this.tittle, ", TF model without intra-motif methylation") }
  
  p1 <- p_motif_coefs + phist +  p_motif_ht  +
        plot_layout( design = layout ) +
        plot_annotation(title = this.tittle )

  ggsave( filename = paste0( prefix_iteration, "_logo_acc_coeffs_motif_ht.pdf" ),
          p1, height = 10, width = 14 )
  
  ggsave( filename = paste0( prefix_iteration, "_only_ht.pdf" ),
          p_motif_ht, height = 10, width = 7 )
}

######################## Format and write start positions across iterations ####
start_pos_list_df <- as.data.frame( start_pos_list )
colnames(start_pos_list_df) <- paste0( "iteration_", 1:ncol(start_pos_list_df) )
rownames(start_pos_list_df) <- rownames(predictors_list[[1]])
start_pos_list_df$peak_id <-rownames(start_pos_list_df)
start_pos_list_df$peak_chr <- gsub( ":.*", "", start_pos_list_df$peak_id )

last_cols <- c( ( ncol(start_pos_list_df) - 1 ), ncol(start_pos_list_df) )

start_pos_list_df <- cbind( start_pos_list_df[,last_cols], start_pos_list_df[,-last_cols] )

write.csv( x = start_pos_list_df, row.names = FALSE, 
           file = paste0( prefix, "_start_position_across_iterations.csv" ) )


####################################################################### End ####
warning()
cat(paste0( "End wall time: ", Sys.time(), "\n"))
sink()





