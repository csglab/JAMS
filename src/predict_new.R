require(stringr)
require(data.table)
require(optparse)
suppressPackageStartupMessages(require(ComplexHeatmap))
suppressPackageStartupMessages(require(circlize))

options( error = traceback, nwarnings = 10000 )
source(  "./src/Methyl_ChIP_ftns.R" )

# setwd( "/home/ahcorcha/repos/tools/JAMS" )
# source(  "./src/Methyl_ChIP_ftns.R" )

option_list = list(
  make_option(c("-c", "--coeffs"), type="character",
              default="./data/CTCF_demo/04_predict/CTCF_HEK293_GSM2026781_small_mCpG_only_model_coefficients_with_FDR_TF_binding.txt",
              help="original coefficients with FDR", metavar="character"),

  make_option(c("-i", "--input_dir"), type="character",
              default="./data/CTCF_demo/04_predict/data_CTCF_search",
              help="", metavar="character"),
  
    make_option(c("-f", "--flanking"), type="integer",
              default=20,
              help="", metavar="character"),
  
    make_option(c("-e", "--experiment"), type="character",
              default="CTCF_HEK293_GSM2026781_small",
              help="", metavar="character"), 
  
    make_option(c("-o", "--output_dir"), type="character",
              default="./data/CTCF_demo/04_predict",
              help="", metavar="character")
  
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)
flanking <- opt$flanking

# opt$input_dir <- "./data/CTCF_demo/04_predict"
# input_root <- "./data/CTCF_demo/04_predict"

coeffs <- read.table( file = opt$coeffs, header = TRUE, sep = "\t" )
in_data <- load_data(opt$input_dir)


prepare_X <- function( in_data, flanking ){
  
  # seq <- in_data$seq
  # CpGs <- in_data$CpGs
  # methyl <- in_data$methyl
  # nonmethyl <- in_data$nonmethyl
  # dna_acc <- in_data$acc
  # target <- in_data$target
  # pfm_length <- in_data$pfm_length

  
  ## Preparation
  # extract the sequence matrix and convert it to one-hot encoding
  x.A <- ( in_data$seq[,-1] == "A" ) * 1
  x.C <- ( in_data$seq[,-1] == "C" ) * 1
  x.G <- ( in_data$seq[,-1] == "G" ) * 1
  x.T <- ( in_data$seq[,-1] == "T" ) * 1
  rownames(in_data$CpGs) <- in_data$CpGs$Name
  x.CpG <- in_data$CpGs[,-1] * 1
  
  ncolx <- ncol(x.A)
  x.CG <- x.C[,-ncolx] * x.G[,-1]
  
  ## Calculate the fraction of methylated C's in the reverse strand (in front of every G)
  x.W <- x.G * ( in_data$methyl[,-1] + 1 ) / ( in_data$methyl[,-1] + in_data$nonmethyl[,-1] + 2 ) 
  ## Calculate the fraction of methylated C's in the forward strand
  x.M <- x.C * ( in_data$methyl[,-1] + 1 ) / ( in_data$methyl[,-1] + in_data$nonmethyl[,-1] + 2 ) 
  
  ## Methylation in reverse strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.W[ is.na(x.W) ] <- 0
  # only keep the W's that are within a CpG
  x.W <- x.W * x.CpG
  ## Methylation in forward strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.M[ is.na(x.M) ] <- 0
  # only keep the W's that are within a CpG
  x.M <- x.M * x.CpG
  
  
  ## Methylation average M (forward) W (reverse).
  x.Met <- ( ( x.M[,-ncolx] + x.W[,-1] ) / 2 )
  
  
  motif <- paste0( "pos(", 0:(in_data$pfm_length - 1) , ")" )
  upsteam_flank <- paste0( "pos(", (-flanking):-1 , ")" )
  downsteam_flank <- paste0( "pos(", (in_data$pfm_length ):(in_data$pfm_length + flanking - 1), ")" )

  
  x.M_up <- rowMeans( x.M[, upsteam_flank] )
  x.W_up <- rowMeans( x.W[, upsteam_flank] )
  x.T_up <- rowSums( x.T[, upsteam_flank] )
  x.C_up <- rowSums( x.C[, upsteam_flank] )
  x.G_up <- rowSums( x.G[, upsteam_flank] )
  
  x.M_down <- rowMeans( x.M[, downsteam_flank] )
  x.W_down <- rowMeans( x.W[, downsteam_flank] )
  x.T_down <- rowSums(  x.T[, downsteam_flank] )
  x.C_down <- rowSums(  x.C[, downsteam_flank] )
  x.G_down <- rowSums(  x.G[, downsteam_flank] )
  
  x.W <- x.W[, motif]; x.M <- x.M[, motif]
  x.C <- x.C[, motif]; x.A <- x.A[, motif]
  x.G <- x.G[, motif]; x.T <- x.T[, motif]
  x.CpG <- x.CpG[, motif]
  
  
  x.Met <- x.Met[, motif]
  x.CG <- x.CG[, motif]

  
  ## DNA accessibility
  acc <- calculate_mean_by_bins(acc = in_data$acc, pfm_length = in_data$pfm_length)
  
  X <- data.frame( acc = acc[,-1], 
                   x.C = x.C, x.G = x.G, x.T = x.T,
                   x.CG = x.CG, x.Met = x.Met,
                   x.T_up = x.T_up, x.C_up = x.C_up, x.G_up = x.G_up,
                   x.M_up = x.M_up, x.W_up = x.W_up,
                   x.T_down = x.T_down, x.C_down = x.C_down,x.G_down = x.G_down,
                   x.M_down = x.M_down, x.W_down = x.W_down)
  
  return(X)
}


X <- prepare_X( in_data, flanking )

## Make sure those variables are in the correct order for the matrix*vector
X <- as.matrix( X[, as.vector(coeffs$Coefficient_name) ] )


log_c_predicted <- as.vector( X %*% coeffs$Estimate )

log_c_predicted <- as.data.frame( log_c_predicted )

rownames(log_c_predicted) <- rownames(X)
log_c_predicted$region <- rownames(X)

# colnames(log_c_predicted) <- c( "region_name", "log_c" )

write.table( x = log_c_predicted[,c( "region", "log_c_predicted")], 
             file = paste0( opt$output_dir, "/", opt$experiment, "_predicted_binding.tab" ), 
             sep = "\t", quote = FALSE,
             row.names = FALSE, col.names = TRUE)


all_data_ht <- function( all_data = all_data, flanking = flanking, htmp_name = htmp_name ){
  
  all_data <- all_data[ order(all_data$log_c_predicted, decreasing = TRUE), ]
  
  acc_col <- c( "acc.bin_up_5", "acc.bin_up_4", 
                "acc.bin_up_3", "acc.bin_up_2", 
                "acc.bin_up_1", "acc.motif_bin", "acc.bin_down_1", 
                "acc.bin_down_2", "acc.bin_down_3", 
                "acc.bin_down_4", "acc.bin_down_5" )
  
  pfm_length <- sum( str_count(colnames(all_data), pattern = "CG") )
  
  meth_col <- paste0( "x.Met.pos.", 0:(pfm_length-1), "." )
  
  C_col <- paste0( "x.C.pos.", 0:(pfm_length-1), "." )
  G_col <- paste0( "x.G.pos.", 0:(pfm_length-1), "." ) # x.G.pos.0.
  T_col <- paste0( "x.T.pos.", 0:(pfm_length-1), "." )
  
  Gs <- all_data[, G_col]
  Cs <- all_data[, C_col]
  Ts <- all_data[, T_col]
  
  Gs[Gs == 1] <- 2
  Cs[Cs == 1] <- 3
  
  seqs <- Gs + Cs + Ts
  
  seqs[seqs == 0] <- "A"
  seqs[seqs == 1] <- "T"
  seqs[seqs == 2] <- "G"
  seqs[seqs == 3] <- "C"
  
  colnames( seqs ) <- paste0("pos_", 1:(pfm_length))
  
  meths <- all_data[, meth_col]
  colnames( meths ) <- paste0("pos_", 1:(pfm_length))
  
  acc <- all_data[, acc_col]
  colnames(acc) <- gsub("acc.bin_", "", colnames(acc) )
  colnames(acc) <- gsub("acc.", "", colnames(acc) )
  
  colors_seq <- structure( c("green", "orange", "blue", "red"), 
                           names = c( "A", "G", "C", "T" ) )
  
  colors_meth <- colorRamp2( c( 0, 0.5, 1 ), c( "white", "red", "black" ) )
  
  
  mid <- min( acc ) + abs( max( acc ) -  min( acc ) ) / 2
  
  colors_acc <- colorRamp2( c( min( acc ), mid, max( acc ) ), 
                            c( "white", "red", "black" ) )
    
  mid <- min( all_data$log_c_predicted ) + 
          abs( max( all_data$log_c_predicted ) -  
                 min( all_data$log_c_predicted ) ) / 2
  
  colors_pred <- colorRamp2( c( min( all_data$log_c_predicted ), 
                                mid,
                                max( all_data$log_c_predicted ) ),
                             c( "white", "red", "black" ) )
  
  pdf( htmp_name, width = 15, height =  20 )
  ## sequence
  ht1 <- ComplexHeatmap::Heatmap( matrix = as.matrix( seqs ),
                                  column_title_rot = 0,
                                  column_title = "Motif sequence",
                                  column_title_side = "bottom",
                                  cluster_columns = F, cluster_rows = F,
                                  col = colors_seq,
                                  show_row_names = F, show_column_names = T,
                                  name = "Sequence", raster_device = "png",
                                  use_raster = T, raster_quality = 2 )
  ## DNA acc.
  ht2 <- ComplexHeatmap::Heatmap( matrix = as.matrix( acc ),
                                  column_title_rot = 0,
                                  column_title = "DNA accessibility",
                                  column_title_side = "bottom",
                                  cluster_columns = F, cluster_rows = F,
                                  col = colors_acc,
                                  show_row_names = F, show_column_names = T,
                                  name = "DNA acc.", raster_device = "png",
                                  use_raster = T, raster_quality = 2 )
  ## Methylation
  ht3 <- ComplexHeatmap::Heatmap( matrix = as.matrix( meths ),
                                  column_title_rot = 0,
                                  column_title = "Methylation",
                                  column_title_side = "bottom",
                                  cluster_columns = F, cluster_rows = F,
                                  col = colors_meth,
                                  show_row_names = F, show_column_names = T,
                                  name = "Methylation", raster_device = "png",
                                  use_raster = T, raster_quality = 2 )
  ## log(c)
  ht4 <- ComplexHeatmap::Heatmap( matrix = as.matrix( all_data$log_c_predicted ),
                                  column_title_rot = 0,
                                  column_title = "",
                                  column_title_side = "bottom",
                                  cluster_columns = F, cluster_rows = F,
                                  col = colors_pred,
                                  show_row_names = F, show_column_names = T,
                                  name = "Predicted\nbinding", raster_device = "png",
                                  use_raster = T, raster_quality = 2 )
   
  ht_list <- ht2 + ht1 + ht3 + ht4
  draw( ht_list, column_title = paste0( "Predicted binding and predictors, n = ", 
                                        nrow(all_dat) ) )
  dev.off() }

X <- as.data.frame(X)
X$name <- rownames(X)

all_dat <- merge(x = X, y = log_c_predicted, by.x = "name", by.y = "region")
rownames(all_dat) <- all_dat$name

write.table( x = all_dat, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
             file = paste0( opt$output_dir, "/", opt$experiment, "_predicted_binding_n_predictors.tab" )
             )

if ( nrow(all_dat) >= 1000){
  all_dat <- all_dat[ sample(1:nrow(all_dat), 1000), ]
  }

htmp_name <- paste0( opt$output_dir, "/", opt$experiment, "_predicted_binding_heatmap.pdf" )
all_data_ht( all_data = all_dat, flanking = flanking, htmp_name = htmp_name )








