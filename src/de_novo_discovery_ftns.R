#### Notes
###############################################################################
# For the de novo motif finding, 
# I think the main idea is to “dynamically” align the peaks based on the
# motif/accessibility/methylation score of the JAMS model, retrain the model, 
# and repeat. So, something like this:
# 
# 1. Start with peak sequence/accessibility/methylation matrices that 
#    are initially aligned simply based on the peak center.
# 
# 2. Train the JAMS TF and background parameters.
# 
# 3. For each peak, identify the position that would result in the 
#    maximum JAMS TF score if it was used as the “center point”.
# 
# 4. Re-align the peaks based on the newly identified center points.
#    By re-aligning, I mean simply sliding each row of the input 
#    matrices to the left or right, so that the center 
#    points are now aligned.
# 
# 5. Repeat 3-4, until convergence.
################################################################################


###################################################################### #
load_dat <- function( input_root, pfm_length = 10 ){
  
  # pfm_length <- 8
  # input_root <- input_root
  
  ###################################################################### #
  fileName <- list.files( path=input_root,full.names = T, 
                          pattern="*_aligned_sequences_tabulated_mx.txt")
  seq <- fread(fileName,sep="\t", data.table = F, header = F)
  
  rownames(seq) <- seq$V1
  seq <- seq[, -c(1,2,3,4,5)]
  
  Cs <- (seq=="C") * 1 # find all the Cs
  Gs <- (seq=="G") * 1 # find all the Gs
  CpGs <- cbind( Cs[,-ncol(Cs)] * Gs[,-1], rep(0,nrow(Cs)) ) 
  ## find all Cs that are followed by a G (the C in a CpG)
  #   also find any entry that follows a mark from the previous step 
  #   (i.e. also mark the Gs in a CpG)
  CpGs[,-1] <- CpGs[,-1] + CpGs[,-ncol(CpGs)] 
  
  rownames(CpGs) <- rownames(seq)
  
  # remove sequences with Ns. These sequences will later be removed 
  #  from the CpG matrix based on sequence name
  #  seq <- seq[ apply( seq[,-1], 1, function(x) sum(x=="N") ) == 0, ]
    
  # Set the column names
  colnames(seq) <- c(paste0("pos_",1:ncol(seq) ) ) # ah_flank
  colnames(CpGs) <- c(paste0("pos_",1:ncol(CpGs) ) ) # ah_flank
  
  
  
  ### Read acc
  ##############################################################################
  fileName <- list.files( path=input_root,
                          pattern=paste0( 
                          "*_aligned_sequences.fasta.accessibility" ), 
                          full.names = T )
  acc <- fread(fileName,sep="\t", data.table = F, header = F)
  # remove peaks with all NAs in the range
  acc <- acc[complete.cases(acc), ]                         
  
  acc_names <- acc$V1
  acc <- acc[,6:ncol(acc)]
  
  ## Replace 0 for pseudocount. As we will use log of DNA_acessibility 
  #   for model and visualization.
  acc[acc == 0] <- NA # All the columns except for the names.
  minimum <- min( acc, na.rm = TRUE )
  pseudocount <- median( sort(as.numeric(as.vector(as.matrix( acc))) ), na.rm = TRUE ) * 0.01
  
  cat( paste0( "\nMinimum value for DNA accessibility (before average or log): ", 
               minimum, "\n", "Pseudocount: ", pseudocount, "\n" ) )
  acc[ is.na( acc ) ] <- pseudocount
  
  ## We are assuming that the input files are centered 
  # take the part of the sequence that's in range
  center <- as.integer(floor(ncol(acc)/2))
  
  rownames(acc) <- acc_names
  # acc$Name <- row.names(acc)
  
  ## 1000 covers the bins in each side, 100 (200) on total covers the searching 
  ## region, centered at the summit (actual PWM hit on the test)
  left_lim <- center - 1100
  right_lim <- center + pfm_length + 1100 # the center position is the PWM start
  
  acc <- acc[, left_lim:right_lim ]
  
  ##############################################################################  
  
  
  
  
  
  ###################################################################### #
  fileName <- list.files(path=input_root, full.names = T, 
                         pattern=paste0("*_sequences.fasta.methylreads"))
  methyl <- fread(fileName,sep="\t", data.table = F, header = F)
  
  rownames(methyl) <- methyl$V1
  methyl <- methyl[, -c(1,2,3,4,5) ]
  methyl <- methyl[ apply( methyl, 1, function(x) sum(!is.na(x)) ) > 0, ]
  # Set the column names
  colnames(methyl) <- c(paste0("pos_",1:ncol(methyl) ) ) # ah_flank
  
    
  ###################################################################### #
  fileName <- list.files(path=input_root, 
                         pattern=paste0("*_aligned_sequences.fasta.nonmethylreads"),
                         full.names = T)
  nonmethyl <- fread(fileName,sep="\t", data.table = F, header = F)
  
  rownames(nonmethyl) <- nonmethyl$V1
  nonmethyl <- nonmethyl[, -c(1,2,3,4,5) ]
  
  # remove peaks with all NAs in the range  
  nonmethyl <- nonmethyl[ apply( nonmethyl, 1, function(x) sum(!is.na(x)) ) > 0, ]
  # Set the column names
  colnames(nonmethyl) <- c(paste0("pos_",1:ncol(nonmethyl) ) ) # ah_flank
  
  
  fileName <- list.files( path = input_root, full.names = TRUE,
                          pattern = "*_summits_vRepeats_scores.txt" )
  target <- fread( fileName, sep = "\t", data.table = FALSE, header = TRUE )
  rownames(target) <- target$Name
  
  include <- Reduce( intersect, list( rownames(seq), 
                                      rownames(methyl), 
                                      rownames(nonmethyl), 
                                      rownames(acc),
                                      rownames(target) ) )
  
  # include <- include[1:10000] # test with smaller number of peaks.
  acc <- acc[ include , ]
  seq <- seq[ include , ]
  CpGs <- CpGs[ include , ]
  methyl <- methyl[ include , ]
  nonmethyl <- nonmethyl[ include , ]
  target <- target[ include, ]
  
  ## Methyl and non methyl reads are not over A or T, only over G and C.
  #####
  if( sum(is.na(methyl[seq=="C"])) > 0 | sum(is.na(methyl[seq=="G"])) > 0 |
      sum(is.na(nonmethyl[seq=="C"])) > 0 | sum(is.na(nonmethyl[seq=="G"])) > 0 ) {
    stop("NA's found where methylation read counts were expected.")
  }
  if( sum(!is.na(methyl[seq=="A"])) > 0 | sum(!is.na(methyl[seq=="T"])) > 0 |
      sum(!is.na(nonmethyl[seq=="A"])) > 0 | sum(!is.na(nonmethyl[seq=="T"])) > 0 ) {
    stop("Methylation read counts found where NA's were expected.")
  }
  #####
  
  
  ###### Prepare data
  ## Preparation
  # extract the sequence matrix and convert it to one-hot encoding
  x.A.all <- ( seq == "A" ) * 1
  x.C.all <- ( seq == "C" ) * 1
  x.G.all <- ( seq == "G" ) * 1
  x.T.all <- ( seq == "T" ) * 1
  # rownames(CpGs) <- CpGs$Name
  x.CpG.all <- CpGs * 1
  
  ncolx <- ncol(x.A.all)
  x.CG.all <- x.C.all[,-ncolx] * x.G.all[,-1]
  
  ## Calculate the fraction of methylated C's in the reverse strand 
  ##  (in front of every G)
  x.W.all <- x.G.all * ( methyl + 1 ) / ( methyl + nonmethyl + 2 ) 
  ## Calculate the fraction of methylated C's in the forward strand
  x.M.all <- x.C.all * ( methyl + 1 ) / ( methyl + nonmethyl + 2 ) 
  
  ## Methylation in reverse strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.W.all[ is.na(x.W.all) ] <- 0
  # only keep the W's that are within a CpG
  x.W.all <- x.W.all * x.CpG.all
  ## Methylation in forward strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.M.all[ is.na(x.M.all) ] <- 0
  # only keep the W's that are within a CpG
  x.M.all <- x.M.all * x.CpG.all
  
  ## Methylation average M (forward) W (reverse).
  x.Met.all <- ( ( x.M.all[,-ncolx] + x.W.all[,-1] ) / 2 )
  
  
  rm(CpGs, ncolx, fileName, include, pseudocount, center, acc_names, seq, Gs, Cs, 
     methyl, nonmethyl, input_root, minimum )
  
  target$log_c_observed <- log(target$pulldown.tag.old / target$ctrl.tag.old)
  
  dat_all <- list( acc = as.data.frame(acc),
                   x.Met.all = as.data.frame(x.Met.all),
                   x.A.all = as.data.frame(x.A.all),
                   x.C.all = as.data.frame(x.C.all), 
                   x.G.all = as.data.frame(x.G.all),
                   x.T.all = as.data.frame(x.T.all), 
                   x.CpG.all = as.data.frame(x.CpG.all),
                   x.CG.all = as.data.frame(x.CG.all),
                   target = as.data.frame(target),
                   x.M.all = as.data.frame(x.M.all), 
                   x.W.all = as.data.frame(x.W.all ) )
  
  rm( x.Met.all, x.A.all, x.C.all, x.G.all, x.T.all, x.CpG.all, x.CG.all, 
      acc, x.M.all, x.W.all )
  
  return( dat_all )
}

############################################################################## #
shift_per_row <- function( shift_pos, df, region_len ){
  
  new_df <- shift_per_row_py( as.vector( shift_pos ), 
                              as.matrix( df ), 
                              as.integer( region_len ) )
  
  rownames( new_df ) <- rownames( df ) 
  
  return( new_df ) }

############################################################################## #
# This function receives an acc matrix that is already shifted 
# by different values per row (by shift per row function).
get_bin_acc <- function( acc = acc, pfm_length = pfm_length ){
  
    acc <- data.frame( bin_down_5 = rowMeans( acc[,1:200] ),
                     bin_down_4 = rowMeans( acc[,200:400] ),
                     bin_down_3 = rowMeans( acc[,400:600] ),
                     bin_down_2 = rowMeans( acc[,600:800] ),
                     bin_down_1 = rowMeans( acc[,800:1000] ),
                     bin_motif =  rowMeans( acc[,1000:(1000+pfm_length)] ),
                     bin_up_1 =   rowMeans( acc[,(1000+pfm_length):(1200+pfm_length)] ),
                     bin_up_2 =   rowMeans( acc[,(1200+pfm_length):(1400+pfm_length)] ),
                     bin_up_3 =   rowMeans( acc[,(1400+pfm_length):(1600+pfm_length)] ),
                     bin_up_4 =   rowMeans( acc[,(1600+pfm_length):(1800+pfm_length)] ),
                     bin_up_5 =   rowMeans( acc[,(1800+pfm_length):(2000+pfm_length)] ))
  return(acc) }


############################################################################## #
# 
# 
train_GLM_at_shifted_pos <- function( flanking, pfm_length, dat_all, start_pos, exclude_meth ){
  
  # flanking <- flanking
  # pfm_length <- pfm_length
  # dat_all <- dat_all
  # exclude_meth <- FALSE
  
  region_len <- 2*flanking + pfm_length
  col_names <- paste0( "pos.", 1:(pfm_length), "." )
  
  start_extract_pos <- start_pos - flanking 
  
  x.A <- shift_per_row( start_extract_pos, dat_all$x.A.all, region_len)
  x.C <- shift_per_row( start_extract_pos, dat_all$x.C.all, region_len)
  x.G <- shift_per_row( start_extract_pos, dat_all$x.G.all, region_len)
  x.T <- shift_per_row( start_extract_pos, dat_all$x.T.all, region_len)
  x.CpG <- shift_per_row( start_extract_pos, dat_all$x.CpG.all, region_len)
  x.CG <- shift_per_row( start_extract_pos, dat_all$x.CG.all, region_len)
  x.Met <- shift_per_row( start_extract_pos, dat_all$x.Met.all, region_len)
  x.M <- shift_per_row( start_extract_pos, dat_all$x.M.all, region_len)
  x.W <- shift_per_row( start_extract_pos, dat_all$x.W.all, region_len)
 
  
  region_len <- 2000 + pfm_length
  acc_start_extract_pos <- start_pos
  acc <- shift_per_row( acc_start_extract_pos, dat_all$acc, region_len )
  acc <- get_bin_acc( acc, pfm_length = pfm_length )
  
  
  upstream_flank <- 1:(flanking)
  motif <- (flanking+1):(flanking+pfm_length)
  downstream_flank <- (flanking+pfm_length+1):(2*flanking+pfm_length)

  #### 
  x.T_up <- rowSums( x.T[, upstream_flank] )
  x.C_up <- rowSums( x.C[, upstream_flank] )
  x.G_up <- rowSums( x.G[, upstream_flank] )
  x.M_up <- rowMeans( x.M[, upstream_flank] )
  x.W_up <- rowMeans( x.W[, upstream_flank] )
  
  x.T_down <- rowSums(x.T[, downstream_flank])
  x.C_down <- rowSums(x.C[, downstream_flank])
  x.G_down <- rowSums(x.G[, downstream_flank])
  x.M_down <- rowSums(x.M[, downstream_flank])
  x.W_down <- rowSums(x.W[, downstream_flank])
  
  x.C <- x.C[, motif ]
  x.A <- x.A[, motif ]
  x.G <- x.G[, motif ]
  x.T <- x.T[, motif ]
  x.CpG <- x.CpG[, motif ]
  x.Met <- x.Met[, motif ]
  x.CG <- x.CG[, motif ]
  
  colnames( x.A ) <- col_names
  colnames( x.C ) <- col_names
  colnames( x.G ) <- col_names
  colnames( x.T ) <- col_names
  colnames( x.CpG ) <- col_names
  colnames( x.CG ) <- col_names
  colnames( x.Met ) <- col_names
  
  c_tags <- unlist( c( dat_all$target$ctrl.tag.old, 
                       dat_all$target$pulldown.tag.old ) )

  t <- unlist( c( rep( 0, nrow( dat_all$target ) ), 
                  rep( 1, nrow( dat_all$target ) ) ) )
  
  XX <- data.frame( acc = acc, x.C = x.C, x.G = x.G, x.T = x.T, 
                    x.CG = x.CG, x.Met = x.Met, x.A = x.A,
                    x.T_up = x.T_up, x.C_up = x.C_up, 
                    x.G_up = x.G_up,x.M_up = x.M_up, 
                    x.W_up = x.W_up, x.T_down = x.T_down, 
                    x.C_down = x.C_down, x.G_down = x.G_down,
                    x.M_down = x.M_down, x.W_down = x.W_down, 
                    stringsAsFactors = TRUE )
  
  rownames( XX ) <- NULL
  XX <- rbind( XX, XX )
  rownames( XX ) <- c( paste0( "control.", dat_all$target$Name ), 
                       paste0( "pulldown.", dat_all$target$Name ) )
  
  
  predictor.names.ctrl <- colnames( XX )
  
  # ahcorcha
  predictor.names.ctrl <- predictor.names.ctrl[
    !grepl("x\\.A\\.pos\\..*\\.$", predictor.names.ctrl)] 
  
  predictor.names.pdwn <- paste0( predictor.names.ctrl, ":t" )
  
  
  if ( exclude_meth ) {
    predictor.names.pdwn <- predictor.names.pdwn[
      !grepl("x\\.Met\\.pos\\..*\\.:t$", predictor.names.pdwn)] 
    }
  
  predictor.names.pdwn <- paste(predictor.names.pdwn, collapse = "+")
  predictor.names.ctrl <- paste(predictor.names.ctrl, collapse = "+")
  
  # c_tags ~ . + t + .:t
  this.formula <- as.formula( 
                   paste0( "c_tags ~ ", predictor.names.ctrl, "+t+",
                            predictor.names.pdwn) )
  
  return( MASS::glm.nb( this.formula, data = XX ) )
}


############################################################################## #
# 
# Create the motif from the regression coefficients
# Store the coefficients and their associated statistics
write.sequence.model.av.met <- function( seq_fit, exclude_meth ) {

  # label <- paste0( experiment, "_", iteration_name ) 
  # seq_fit <- this_glm
  # exclude_meth <- FALSE
  # outdir <- "./data/de_novo_discovery_test"
  
  coefs <- coefficients( summary( seq_fit ) )
  
  ## In case there are missing values in the coefs matrix, add them back
  coefs <- fix.coefs( coefs, seq_fit$coefficients )
  
  ## Calculate the FDR for each coefficient
  coefs <- cbind( coefs, p.adjust( coefs[,4], method = "fdr" ) ) 
  colnames(coefs)[5] <- "FDR" # add the FDR to the matrix
  if(exclude_meth){ coefs <- add.meth.coeffs( coefs, pfm_length )  }
  
  # write.table( x =  coefs, quote = FALSE, sep = "\t", row.names = TRUE, 
  #              col.names = TRUE, 
  #              file = paste0( prefix, "_coefficients_with_FDR.txt") )
  
  
  ## Get coef. for bases and methylation.
  coefs_bases_control <- coefs[ grepl( pattern = ".\\.pos\\..*\\.$", x = rownames(coefs)), ]
  coefs_bases_control <- data.frame( coefs_bases_control )
  
  
  ## Get coef. for bases and methylation.
  coefs_bases_pulldown <- coefs[grepl( pattern = ".\\.pos\\..*\\.:t$", x = rownames(coefs)), ]
  coefs_bases_pulldown <- data.frame( coefs_bases_pulldown )
  
  
  # Get base/bases specific per position
  base <- gsub( pattern = "x.", replacement = "", x = rownames( coefs_bases_control ) )
  coefs_bases_control$base <- gsub( pattern = "\\.pos\\..*\\.", replacement = "", x = base )
  
  base <- gsub( pattern = "x.", replacement = "", x = rownames( coefs_bases_pulldown ) )
  coefs_bases_pulldown$base <- gsub( pattern = "\\.pos\\..*\\.:t", replacement = "", x = base )
  
  
  coefs_control_meth <- coefs_bases_control[ ( coefs_bases_control$base == "Met" ), ]
  coefs_pulldown_meth <- coefs_bases_pulldown[ ( coefs_bases_pulldown$base == "Met" ), ]
  
  
  ## 
  coefs_control_meth$FDR[ is.na( coefs_control_meth$FDR ) ] <- 1
  coefs_control_meth[ coefs_control_meth$FDR >= 0.00001, ]$Estimate <- NA

  coefs_pulldown_meth$FDR[ is.na( coefs_pulldown_meth$FDR ) ] <- 1
  coefs_pulldown_meth[ coefs_pulldown_meth$FDR >= 0.00001, ]$Estimate <- NA
  
  
  # Take out CpG position coefficients
  coefs_bases_control <- coefs_bases_control[ ( coefs_bases_control$base != "Met" ) & 
                                              ( coefs_bases_control$base != "CG" ), ]
  
  coefs_bases_pulldown <- coefs_bases_pulldown[ ( coefs_bases_pulldown$base != "Met" ) & 
                                                ( coefs_bases_pulldown$base != "CG" ), ]
  
  
  #### Control
  motif_control <- matrix( coefs_bases_control[, "Estimate"], nrow = 3, byrow = TRUE )
  ## Add a row in the beginning, corresponding to A
  motif_control <- rbind( rep.int( 0, ncol( motif_control ) ), motif_control )  
  motif_control <- rbind( motif_control, coefs_control_meth$Estimate )
  
  rownames( motif_control ) <- c( "A", "C", "G", "T", "CtoM" )
  colnames( motif_control ) <- 1:ncol( motif_control )
  motif_control[1:4, ] <- scale( motif_control[1:4,], center = TRUE, scale = FALSE )
  
  
  #### Pulldown
  motif_pulldown <- matrix( coefs_bases_pulldown[, "Estimate"], nrow = 3, byrow = TRUE )
  ## Add a row in the beginning, corresponding to A
  motif_pulldown <- rbind( rep.int( 0, ncol( motif_pulldown ) ), motif_pulldown )  
  motif_pulldown <- rbind( motif_pulldown, coefs_pulldown_meth$Estimate )
  
  rownames( motif_pulldown ) <- c( "A", "C", "G", "T", "CtoM" )
  colnames( motif_pulldown ) <- 1:ncol( motif_pulldown )
  motif_pulldown[1:4, ] <- scale( motif_pulldown[1:4,], center = TRUE, scale = FALSE )
  
  
  ## Scale the motif so that each column has a mean of zero, and also correct the
  #  intercept so that the motif scores would still exactly match the regression results
  ## Show control and pulldown motif in the same figure
  pp_motif <- draw.ctrl.pdwn.motif.av.meth( pulldown_coefs = coefs_pulldown_meth,
                                            control_coefs = coefs_control_meth,
                                            pulldown_motif = motif_pulldown,
                                            control_motif = motif_control )
  return( list(coefs = coefs, pp_motif =pp_motif) )
}


draw.ctrl.pdwn.motif.av.meth <- function( pulldown_coefs = coefs_pulldown_meth,
                                          control_coefs = coefs_control_meth,
                                          pulldown_motif = motif_pulldown,
                                          control_motif = motif_control ) {

  # pulldown_coefs = coefs_pulldown_meth
  # control_coefs = coefs_control_meth
  # pulldown_motif = motif_pulldown
  # control_motif = motif_control
  
  area_scale <- max( max( abs( pulldown_motif[1:4, ] ) ), max( abs( control_motif[1:4, ] ) ) )
  area_scale_ctrl <- max( abs( control_motif[1:4, ] ) )
  
  
  this_width <- max( 2 +  0.85*ncol( pulldown_motif ), 8 )

  pdwn_motif <- motif.logo.av.met( motif = pulldown_motif, coefs_bases = pulldown_coefs )
  pdwn_htmp <- motif.heatmap.av.met( motif = pulldown_motif, coefs_bases = pulldown_coefs,
                                     area_scale = area_scale )

  # ggsave( filename = paste0(outdir,"/",label,"_pulldown_motif_logo.pdf"), 
  #         plot = pdwn_motif / pdwn_htmp, 
  #         units = "cm", device = "pdf", width = this_width, height = 23 )
  

  ctrl_motif <- motif.logo.av.met( motif = control_motif, coefs_bases = control_coefs )
  ctrl_htmp <- motif.heatmap.av.met( motif = control_motif, coefs_bases = control_coefs,
                                     area_scale = area_scale_ctrl )
  
  # ggsave( filename = paste0(outdir,"/",label,"_control_motif_logo.pdf"), plot = ctrl_motif / ctrl_htmp,
  #         units = "cm", device = "pdf", width = this_width, height = 23 )
  
  ctrl_htmp <- motif.heatmap.av.met( motif = control_motif, coefs_bases = control_coefs,
                                     area_scale = area_scale )
  
  
  ## Y positive limit
  pos_y <- cbind( pulldown_motif, control_motif )
  pos_y[ pos_y < 0 ] <- 0
  pos_y <- max( colSums( pos_y, na.rm = TRUE ) )*1

  ## Y negative limit
  neg_y <- cbind( pulldown_motif, control_motif )
  neg_y[ neg_y > 0 ] <- 0
  neg_y <- min( colSums( neg_y, na.rm = TRUE ) )*1
  

  
  pdwn_motif <- pdwn_motif + ylim( neg_y, pos_y) + ggtitle( "Pulldown") + 
                 theme( plot.title = element_text( hjust = 0.5 ) )
  
  ctrl_motif <- ctrl_motif + ylim( neg_y, pos_y)
  
  ctrl_motif <- ctrl_motif + ylab("") + ggtitle( "Control" ) +
                 theme( plot.title = element_text( hjust = 0.5 ) )
  
  ctrl_htmp <- ctrl_htmp + ylab("")
  
  
  pp <- ( ( pdwn_motif / pdwn_htmp ) | ( ctrl_motif  / ctrl_htmp ) ) 
  #       + plot_annotation( theme( title = "Sequence & mCpG coefficients"  ) )
  
  # ggsave( filename = paste0( prefix, "_control_n_pulldown_motif_logo.pdf"), 
  #         plot = pp, units = "cm", device = "pdf", 
  #         width = this_width*2, height = 23 )
  return(pp)
  }


pre_calc_by_pos_dat <- function( this_dat_all, possible_position, flanking = 20, pfm_length = pfm_length ){
  
  ## This should run once
  # this_dat_all <- dat_all
  # possible_position <- possible_position
  # flanking <- 20
  # pfm_length <- pfm_length
  
  by_pos_predictor_lists <- list()
  col_names <- paste0("pos.", 1:pfm_length, "." )
  
  for ( off.pos in 1:possible_position ){
    
    # off.pos <- 1
    if ( off.pos == 1 ){
      acc <- get_bin_acc(this_dat_all$acc, pfm_length=pfm_length)
    } else if ( off.pos == 2 ){
      acc <- get_bin_acc(this_dat_all$acc[,-1], pfm_length=pfm_length)
    } else{
      acc <- get_bin_acc(this_dat_all$acc[,-1:-off.pos], pfm_length=pfm_length)
    }
    
    #### 
    upstream_flank <- (off.pos):(off.pos+flanking-1)
    motif <- (off.pos+flanking):(off.pos+flanking+pfm_length-1)
    downstream_flank <- (off.pos+flanking+pfm_length):(off.pos+2*flanking+pfm_length-1)
    
    
    x.T_up <- rowSums( this_dat_all$x.T.all[, upstream_flank] )
    x.C_up <- rowSums( this_dat_all$x.C.all[, upstream_flank] )
    x.G_up <- rowSums( this_dat_all$x.G.all[, upstream_flank] )
    x.M_up <- rowMeans( this_dat_all$x.M.all[, upstream_flank] )
    x.W_up <- rowMeans( this_dat_all$x.W.all[, upstream_flank] )
    
    x.T_down <- rowSums(this_dat_all$x.T.all[, downstream_flank])
    x.C_down <- rowSums(this_dat_all$x.C.all[, downstream_flank])
    x.G_down <- rowSums(this_dat_all$x.G.all[, downstream_flank])
    x.M_down <- rowSums(this_dat_all$x.M.all[, downstream_flank])
    x.W_down <- rowSums(this_dat_all$x.W.all[, downstream_flank])
    
    x.C <- this_dat_all$x.C.all[, motif ]
    x.A <- this_dat_all$x.A.all[, motif ]
    x.G <- this_dat_all$x.G.all[, motif ]
    x.T <- this_dat_all$x.T.all[, motif ]
    x.CpG <- this_dat_all$x.CpG.all[, motif ]
    x.Met <- this_dat_all$x.Met.all[, motif ]
    x.CG <- this_dat_all$x.CG.all[, motif ]
    
    colnames( x.A ) <- col_names
    colnames( x.C ) <- col_names
    colnames( x.G ) <- col_names
    colnames( x.T ) <- col_names
    colnames( x.CpG ) <- col_names
    colnames( x.CG ) <- col_names
    colnames( x.Met ) <- col_names
    
    X <- data.frame( acc = acc, x.A = x.A,
                     x.C = x.C, x.G = x.G, x.T = x.T,
                     x.CG = x.CG, x.Met = x.Met,
                     x.T_up = x.T_up, x.C_up = x.C_up, 
                     x.G_up = x.G_up,x.M_up = x.M_up, 
                     x.W_up = x.W_up, x.T_down = x.T_down, 
                     x.C_down = x.C_down, x.G_down = x.G_down,
                     x.M_down = x.M_down, x.W_down = x.W_down, 
                     stringsAsFactors = TRUE )
    
    len <- length(by_pos_predictor_lists)
    by_pos_predictor_lists[[len+1]] <- X
    
    
    
    }
  return( by_pos_predictor_lists )
}


plot_dna_acc_coefficients <- function(fit_model){
  
  # fit_model <- fit_CpG_only$fit
  # plot_name <- dna_acc_plot_name
  
  ## Subset coefficients
  coefs <- as.data.frame( coef( summary( fit_model ) ) )
  coefs$FDR <- as.double( p.adjust( coefs[,4], method = "fdr" ) )
  coefs$FDR[ coefs$FDR >= 0.1 ] <- NA
  
  coefs$FDR_color[ is.na( coefs$FDR ) ] <- "< 0.1"
  coefs$FDR_color[ ! is.na( coefs$FDR ) ] <- "> 0.1"
  
  dna_coefs <- coefs[ grepl( pattern = "acc", x = rownames( coefs ) ), ]
  
  dna_coefs_ctrl <- dna_coefs[ ! grepl( pattern = ":t", x = rownames( dna_coefs ) ), ]
  dna_coefs_ctrl$Name <- rownames( dna_coefs_ctrl )
  dna_coefs_ctrl$Name <- gsub( "_", "\\.", dna_coefs_ctrl$Name )
  
  dna_coefs_pdwn <- dna_coefs[ grepl( pattern = ":t", x = rownames( dna_coefs ) ), ]
  dna_coefs_pdwn$Name <- rownames( dna_coefs_pdwn )
  dna_coefs_pdwn$Name <- gsub( "_", "\\.", dna_coefs_pdwn$Name )
  
  
  dna_coefs_pdwn$Name <- factor(dna_coefs_pdwn$Name, levels = dna_coefs_pdwn$Name)
  dna_coefs_ctrl$Name <- factor(dna_coefs_ctrl$Name, levels = dna_coefs_ctrl$Name)
  
    
  
  #################################################################### ##
  y_min <- min( c( ( dna_coefs_ctrl$Estimate - dna_coefs_ctrl$`Std. Error` ), 
                   ( dna_coefs_pdwn$Estimate - dna_coefs_pdwn$`Std. Error` ) ) )
  
  y_max <- max( c( ( dna_coefs_ctrl$Estimate + dna_coefs_ctrl$`Std. Error` ), 
                   ( dna_coefs_pdwn$Estimate + dna_coefs_pdwn$`Std. Error` ) ) )
  
  mag <- abs( y_max - y_min )
  center <- y_min + ( mag / 2 )
  
  y_min <- center - ( mag / 2 )*1.2
  y_max <- center + ( mag / 2 )*1.2
  #################################################################### ## 
  
  
  dna_pdwn_plot <- ggplot( data = dna_coefs_pdwn, aes( x = Name, y = Estimate, 
                                                       color = FDR_color ) ) +
                           geom_point( size = 1 ) +
                           geom_errorbar( aes( x = Name, width=0.5,
                                               ymin = Estimate - `Std. Error`, 
                                               ymax = Estimate + `Std. Error` ) ) +
                           theme_minimal() + labs( x = "Pulldown", y = "DNA accessibility\n coefficients" ) +
                           theme( axis.text.x = element_text( angle = 290, hjust = 0 ),
                                  legend.position = "none" ) + 
                           ylim( y_min, y_max ) + 
                           scale_color_manual( values = c( "grey", "black" ) )
  
  
  dna_ctrl_plot <- ggplot( data = dna_coefs_ctrl, aes( x = Name, y = Estimate, 
                                                       color = FDR_color ) ) +
                           geom_point( size = 1 ) +
                           geom_errorbar( aes( x = Name, width=0.5,
                                               ymin = Estimate - `Std. Error`, 
                                               ymax = Estimate + `Std. Error` ) ) +
                           theme_minimal() + labs( x = "Control", y = "" ) + 
                           theme( axis.text.x = element_text( angle = 290, hjust = 0 ) ) +
                           ylim( y_min, y_max ) + 
                           scale_color_manual( values = c( "grey", "black" ),
                                               labels = c( "> 0.1", "< 0.1" ) ) + 
                           labs(color = "FDR")
  
  coefficient_plot <- ( dna_pdwn_plot + dna_ctrl_plot ) + 
                        plot_annotation( title = "GLM coefficients for DNA accessibility" )
  
  # ggsave( filename = plot_name, plot = coefficient_plot, 
  #         units = "cm", width = 20, height = 10 ) 
  return(coefficient_plot)
  }


format_iteration <- function(i){
  i <- as.integer(abs(i))
  if ( ( 1 <= i) & ( i <= 9 ) ){ new_i <- paste0("00", i) }
  if ( ( 10 <= i) & ( i <= 99 ) ){ new_i <- paste0("0", i) }
  if ( ( 100 <= i) & ( i <= 99 ) ){ new_i <- paste0("", i) }
  return( new_i ) }


create_pos_vector <- function(this_start_pos, n_cols, pfm_length){
  return(c( rep( NA, this_start_pos-1 ), 
            rep( "motif", pfm_length ), 
            rep( NA, n_cols - ( this_start_pos+pfm_length ) ) ) ) }


motif_pos_heatmap <- function(this_start_pos, n_cols = 201, pfm_length, iteration  ){
  
  # this_start_pos <- start_pos
  # n_cols <- 201
  
  viz_highest_TF <- as.data.frame( t( sapply( X = this_start_pos, 
                                              FUN = create_pos_vector, 
                                              n_cols = n_cols,
                                              pfm_length = pfm_length ) ) )
  
  centered_line <- list( c( rep(NA, floor(n_cols/2)), 
                            rep("center", pfm_length), 
                            rep(NA, floor(n_cols/2)-pfm_length) ) )
  
  centered_lines <- rep(list(centered_line), round(0.01*length(this_start_pos) ) )
  centered_lines <- t( as.data.frame( centered_lines ) )
  rownames(centered_lines) <- NULL
  viz_highest_TF <- rbind( centered_lines, viz_highest_TF )
  
  col_fun_tf <- structure( c("black", "red" ), names = c( "motif", "center" ) )
  
  this_title <- paste0( "Motif over peak center (+/-200bp)", 
                        ", n = ", length( this_start_pos ) )
  
  ht_tf <- ComplexHeatmap::Heatmap( matrix = as.matrix( viz_highest_TF ),
                                    cluster_columns = F, cluster_rows = F,
                                    show_row_names = F, show_column_names = F,
                                    na_col = "white",
                                    col = col_fun_tf,
                                    # border_gp = gpar(col = "black"),
                                    border = TRUE,
                                    name = " ",
                                    column_title = this_title,
                                    raster_device = "png",
                                    use_raster = T, raster_quality = 2 )
  return(ht_tf) }


eval_coeffs <- function( pos_predictor, pdn_coeff, X_names ){
  return( as.vector( as.matrix( pos_predictor[, X_names]) %*% pdn_coeff ) ) }


add.meth.coeffs <- function( pdwn_coeffs, pfm_length ){
  
  dummy_line1 <- rep_len( x = 0, length.out = pfm_length )
  dummy_line2 <- rep_len( x = 1, length.out = pfm_length )
  
  meth_coeffs <-  data.frame( "Estimate" = dummy_line1,
                              "Std. Error" = dummy_line1,
                               "z value" = dummy_line1,
                               "Pr(>|z|)" = dummy_line2, 
                               "FDR" = dummy_line2 
                              )
  # x.Met.pos.6.
  rownames(meth_coeffs) <- paste0("x.Met.pos.", 1:pfm_length, ".:t" )
    
  colnames(meth_coeffs) <- colnames(pdwn_coeffs)
    
  pdwn_coeffs <- rbind( pdwn_coeffs, meth_coeffs )
    
  return( as.data.frame( pdwn_coeffs ) ) }


rev_complement_predictor <- function( pos_specific_predictor, pfm_length ){
  
  rename_cols <- function( predictor_cols ){
    
    predictor_cols <- gsub( "^acc\\.bin_down_5$", "acc\\.bin_up_5_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^acc\\.bin_down_4$", "acc\\.bin_up_4_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^acc\\.bin_down_3$", "acc\\.bin_up_3_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^acc\\.bin_down_2$", "acc\\.bin_up_2_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^acc\\.bin_down_1$", "acc\\.bin_up_1_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^acc\\.bin_up_1$", "acc\\.bin_down_1_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^acc\\.bin_up_2$", "acc\\.bin_down_2_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^acc\\.bin_up_3$", "acc\\.bin_down_3_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^acc\\.bin_up_4$", "acc\\.bin_down_4_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^acc\\.bin_up_5$", "acc\\.bin_down_5_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^x\\.T_up$", "x\\.T_down_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^x\\.C_up$", "x\\.C_down_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^x\\.G_up$", "x\\.G_down_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^x\\.M_up$", "x\\.M_down_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^x\\.W_up$", "x\\.W_down_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^x\\.T_down$", "x\\.T_up_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^x\\.C_down$", "x\\.C_up_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^x\\.G_down$", "x\\.G_up_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^x\\.M_down$", "x\\.M_up_asdsd_tmp", predictor_cols )
    predictor_cols <- gsub( "^x\\.W_down$", "x\\.W_up_asdsd_tmp", predictor_cols )

    predictor_cols <- gsub( "_asdsd_tmp$", "", predictor_cols )
    
  }
    
  # pos_specific_predictor <- predictors_list[[1]]
  
  colnames(pos_specific_predictor) <- rename_cols(colnames(pos_specific_predictor))
    
  
  ############################################## Name seq/CG/meth col names #### 
  A_cols <- paste0("x.A.pos.", 1:pfm_length, ".")
  T_cols <- paste0("x.T.pos.", 1:pfm_length, ".")
  C_cols <- paste0("x.C.pos.", 1:pfm_length, ".")
  G_cols <- paste0("x.G.pos.", 1:pfm_length, ".")
  CG_cols <- paste0("x.CG.pos.", 1:pfm_length, ".")
  Met_cols <- paste0("x.Met.pos.", 1:pfm_length, ".")

  CG_cols_shift <- paste0("x.CG.pos.", 2:(pfm_length+1), ".")
  Met_cols_shift <- paste0("x.Met.pos.", 2:(pfm_length+1), ".")
    
  rev_A_cols <- paste0("x.A.pos.", pfm_length:1, ".")
  rev_T_cols <- paste0("x.T.pos.", pfm_length:1, ".")
  rev_C_cols <- paste0("x.C.pos.", pfm_length:1, ".")
  rev_G_cols <- paste0("x.G.pos.", pfm_length:1, ".")
  rev_CG_cols <- paste0("x.CG.pos.", pfm_length:1, ".")
  rev_Met_cols <- paste0("x.Met.pos.", pfm_length:1, ".")      
  
  ############################## Inplace reverse complement of meth and CpG #### 
  ## Reverses and renames
  pos_specific_predictor[, CG_cols] <- pos_specific_predictor[,rev_CG_cols]
  pos_specific_predictor[, Met_cols] <- pos_specific_predictor[,rev_Met_cols]
  
  pos_specific_predictor <- pos_specific_predictor[, !colnames(pos_specific_predictor) %in% 
                                                      c( CG_cols[1], Met_cols[1]) ]
  
  ## add one column at the end for CpG and Met
  pos_specific_predictor[, CG_cols_shift[length(CG_cols_shift)] ] <- 0
  pos_specific_predictor[, Met_cols_shift[length(Met_cols_shift)] ] <- 0
  
  
  ## Rename columns 2-16 to 1-15  
  colnames(pos_specific_predictor)[colnames(pos_specific_predictor)
            %in% CG_cols_shift] <- CG_cols

  colnames(pos_specific_predictor)[colnames(pos_specific_predictor)
            %in% Met_cols_shift] <- Met_cols    
  
  
  ################################### Get the reverse complement of T/A/C/G ####
  rev_compl_T <- pos_specific_predictor[, rev_A_cols]
  colnames(rev_compl_T) <- T_cols
  
  rev_compl_A <- pos_specific_predictor[, rev_T_cols]
  colnames(rev_compl_A) <- A_cols
  
  rev_compl_C <- pos_specific_predictor[, rev_G_cols]
  colnames(rev_compl_C) <- C_cols
  
  rev_compl_G <- pos_specific_predictor[, rev_C_cols]
  colnames(rev_compl_G) <- G_cols
  
  rev_compl <- cbind(rev_compl_T,rev_compl_A, rev_compl_C, rev_compl_G)
  rm(rev_compl_T,rev_compl_A, rev_compl_C, rev_compl_G); gc(verbose=FALSE)
  
  pos_specific_predictor <- pos_specific_predictor[,
                     !colnames(pos_specific_predictor) %in% colnames(rev_compl)]
  
  return( list( cbind( pos_specific_predictor, rev_compl ) ) )
  
  
  
}


vis_sequence <- function( predictors, pfm_length, ht_path ){
  
  # predictors <- predictors_list[[101]]
  predictors <- head(predictors, n = 5)
  C_bases <- predictors[, paste0( "x.C.pos.", 1:(pfm_length), "." ) ]
  G_bases <- predictors[, paste0( "x.G.pos.", 1:(pfm_length), "." ) ]
  G_bases[] <- base::lapply( G_bases, function(x) as.numeric( gsub("1", "2", x )) )
  
  T_bases <- predictors[, paste0( "x.T.pos.", 1:(pfm_length), "." ) ]
  T_bases[] <- base::lapply( T_bases, function(x) as.numeric( gsub("1", "3", x )) )
  
  bases <- C_bases + G_bases + T_bases
  
  bases[] <- base::lapply( bases, function(x) as.character( gsub("0", "A" , x )) )
  bases[] <- base::lapply( bases, function(x) as.character( gsub("1", "C" , x )) )
  bases[] <- base::lapply( bases, function(x) as.character( gsub("2", "G" , x )) )
  bases[] <- base::lapply( bases, function(x) as.character( gsub("3", "T" , x )) )
  col_seq <- paste0( "seq.pos.", 1:pfm_length )
  colnames(bases) <- col_seq
  
  acc_cols <- c( "acc.bin_down_5", "acc.bin_down_4", "acc.bin_down_3", 
                 "acc.bin_down_2", "acc.bin_down_1", "acc.bin_motif",  
                 "acc.bin_up_1", "acc.bin_up_2", "acc.bin_up_3", 
                 "acc.bin_up_4", "acc.bin_up_5")
  
  acc <- predictors[, acc_cols ]
  
  met <- predictors[, paste0( "x.Met.pos.", 1:pfm_length, "." ) ]
  met_col <- paste0( "met.pos.", 1:pfm_length )
  colnames( met ) <- met_col
  
  res_dat <- cbind( met, bases, acc )
  
  col_fun <- colorRamp2( c( 0, 0.5, 1 ),
                         c( "white", "red" ,"black" ) )

  my_use_raster <- TRUE
  
  htm_met <- ComplexHeatmap::Heatmap( as.matrix( res_dat[, met_col] ), 
                                      use_raster = my_use_raster,
                                      raster_quality = 2,
                                      cluster_columns = FALSE,
                                      cluster_rows = FALSE,
                                      show_column_names =  TRUE,
                                      show_row_names  =  FALSE,
                                      col = col_fun,
                                      name = "CpG Methylation",
                                      column_title = "meCpG" )
  
  #### Sequence
  colors_bases <- structure( c("green", "orange", "blue", "red"), 
                             names = c( "A", "G", "C", "T" ) )
      
  htm_seq <- ComplexHeatmap::Heatmap( as.matrix( res_dat[, col_seq] ), 
                                      use_raster = my_use_raster,
                                      raster_quality = 2,
                                      cluster_columns = FALSE, 
                                      cluster_rows = FALSE,
                                      show_column_names =  TRUE,
                                      show_row_names  =  FALSE,
                                      col = colors_bases,
                                      name = "Sequence",
                                      column_title = "Sequence" )
  
  #### Accessibility
  # diff_tmp <-  abs( max( res_dat[, acc_cols] ) -  min( res_dat[, acc_cols] ) )
  # mid <- min( res_dat[, acc_cols] ) + diff_tmp/2
  col_fun <- colorRamp2( c( 0, 3, 7),
                         c( "white", "red", "black" ) )
  
  htm_acc <- ComplexHeatmap::Heatmap( as.matrix( res_dat[, acc_cols] ), 
                                      use_raster = my_use_raster,
                                      raster_quality = 2,
                                      cluster_columns = FALSE, 
                                      cluster_rows = FALSE,
                                      show_column_names =  TRUE,
                                      show_row_names  =  FALSE,
                                      col = col_fun,
                                      name = "Accessibility",
                                      column_title = "Accessibility" )
  
  all_htm <- htm_seq + htm_acc + htm_met
  
  pdf( file = ht_path, width = 8, height = 5 )
  
  draw( all_htm, ht_gap = unit( 0.25, "cm" ), 
        column_title = paste0( "\nn = ", nrow(res_dat)),
        merge_legends = TRUE,
        heatmap_legend_side = "bottom" )
  
  dev.off()
}


pos_to_wiggle <- function( pdwn_coeffs, shifting_pos = 5, inf_pct = 0.3 ){
  
  # shifting_pos <- 5
  # inf_pct <- 0.3
  # pdwn_coeffs <- pdwn_coeffs
  
  # motif <- get_motif(pdwn_coeffs, magnitud = "Estimate" )
  motif <- get_motif(pdwn_coeffs, magnitud = "z value" )
  
  motif <- as.data.frame(abs( t( motif ) ))
  
  motif$sum <- rowSums( motif[,1:4] )
  
  top <- mean( motif[ order(motif$sum, decreasing = TRUE), "sum" ][1:4] )
  
  motif$pct <- ( motif$sum / top )
  
  mean_start <- mean( motif$pct[1:shifting_pos] )
  mean_stop <- mean( motif$pct[ ( nrow(motif)-shifting_pos ):nrow(motif)] )
  
  if((mean_start > inf_pct)&(mean_stop < inf_pct)){shift_pos <- -shifting_pos}
  if((mean_start < inf_pct)&(mean_stop > inf_pct)){shift_pos <- +shifting_pos}
  if((mean_start > inf_pct)&(mean_stop > inf_pct)){shift_pos <- 0}
  if((mean_start < inf_pct)&(mean_stop < inf_pct)){shift_pos <- 0}
  
  return( shift_pos )
}


get_motif <- function(pdwn_coeffs, magnitud = "Estimate" ){
  
  # magnitud <- "Estimate"
  # magnitud <- "z value"
  C_coefs <- pdwn_coeffs[ grepl( pattern = "^x\\.C\\.pos\\..*\\.:t$", 
                                 x = rownames(pdwn_coeffs)),]
  T_coefs <- pdwn_coeffs[ grepl( pattern = "^x\\.T\\.pos\\..*\\.:t$", 
                                 x = rownames(pdwn_coeffs)),]
  G_coefs <- pdwn_coeffs[ grepl( pattern = "^x\\.G\\.pos\\..*\\.:t$", 
                                 x = rownames(pdwn_coeffs)),]
  
  motif <- data.frame( C_seq = C_coefs[, magnitud ],
                       T_seq = T_coefs[, magnitud ],
                       G_seq = G_coefs[, magnitud ],
                       A_seq = rep.int( 0, nrow( C_coefs ) ) )
  
  motif <- t( motif )
  rownames( motif ) <- c( "C", "T", "G", "A" )
  motif <- scale( motif, center = TRUE, scale = FALSE )
  return(as.data.frame(motif))
}


write.pfm.cisbp <- function(filename, motif, motif_name = "motif_name" ){
  
  # filename = paste0( prefix, "_motif.pfm.txt" )
  # motif = motif
  # motif_name = opt$experiment
  
  header <- paste0( "Gene\tUnknown\n", 
                    "Motif\t", motif_name, "\n",
                    "Family\tUnknown\n", 
                    "Species\tUnknown" )
  
  motif <- rbind( motif, Pos = 1:ncol(motif) )
  
  motif_header <- c("Pos", "A", "C", "G", "T")
  out_motif <- t(motif)[, motif_header ]
  out_motif <- rbind(motif_header, out_motif)
  motif_string <- paste0(out_motif[,1], "\t", 
                         out_motif[,2], "\t",
                         out_motif[,3], "\t",
                         out_motif[,4], "\t",
                         out_motif[,5] )
    
  fileConn <- file( filename )
  writeLines( c( header, motif_string, "\n" ), fileConn)
  close(fileConn)
  
}

trim_motif <- function( motif, base_th = 0.2 ){
  
  motif <- as.data.frame(t(motif))
  motif$sum <- rowSums(abs(motif))
  top <- mean( motif[ order(motif$sum, decreasing = TRUE), "sum" ][1:2] )
  motif$pct <- ( motif$sum / top )
  change_flag <- TRUE
  
  while ( change_flag == TRUE ) {
    
    change_flag <- FALSE
    
    if ( motif[1,"pct"] < base_th ) {
      motif <- motif[-1,]
      change_flag <- TRUE }
    
    if ( motif[ nrow(motif), "pct" ] < base_th ) {
      motif <- motif[-nrow(motif),]
      change_flag <- TRUE }
  }
  
  return(  t( motif[, c( "C", "T", "G", "A" ) ] )  )
}











