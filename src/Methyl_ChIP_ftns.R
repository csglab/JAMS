############################################################################### #
###################### Functions for fitting sequence model ################### #
#################################################################################
# function that creates per-position statistics regarding methylation
# seq: the sequence matrix
# CpGs: the binary matrix denoting CpG locations
# methyl: the read count matrix for methylated C/Gs
# nonmethyl: the read count matrix for non-methylated C/Gs
# cutoff: read count cutoff for considering a positin
##
# fit a model that predicts TF binding probability 
#    (and, therefore, ChIP tag counts) based on sequence and methylation
# seq: the sequence matrix
# CpGs: the binary matrix denoting CpG locations
# methyl: the read count matrix for methylated C/Gs
# nonmethyl: the read count matrix for non-methylated C/Gs
# tags: the tag count per peak
# length: the length of each peak
# CpG_only: whether only DNA methylation at CpG sites should be considered in the model
############################################################################### #
plot_methylation_stats <- function( methyl_chip_data = methyl_chip_data, 
                                    pseudocount = 0.0001, cutoff = 10, 
                                    outdir = outdir, experiment = opt$experiment ){
  # methyl_chip_data <- methyl_chip_data
  # seq <- methyl_chip_data$seq
  # CpGs <- methyl_chip_data$CpGs
  # methyl <- methyl_chip_data$methyl
  # nonmethyl <- methyl_chip_data$nonmethyl
  # pseudocount = 0.0001
  # cutoff <- 10
  # experiment <- opt$experiment  
  
  tags_no_zero <- replace(methyl_chip_data$target$tags, methyl_chip_data$target$tags == 0, pseudocount)
  log_density <- log(tags_no_zero/methyl_chip_data$target$length)
  
  # extract the sequence matrix and convert it to one-hot encoding
  x.C <- ( methyl_chip_data$seq[,-1] == "C" ) * 1
  x.C[x.C == 0] <- NA # replace zeros with NAs
  x.G <- ( methyl_chip_data$seq[,-1] == "G" ) * 1
  
  x.G[x.G == 0] <- NA # replace zeros with NAs
  x.CpG <- methyl_chip_data$CpGs[,-1] * 1 # the matrix of CpGs
  x.nonCpG <- 1 - x.CpG # the things that are not CpG
  x.CpG[x.CpG == 0] <- NA # replace zeros with NAs
  x.nonCpG[x.nonCpG == 0] <- NA # replace zeros with NAs
  
  # calculate the beta matrix, replacing ambiguous values (i.e. places with total read count below cutoff) with NA
  sum.counts <- methyl_chip_data$methyl[,-1] + methyl_chip_data$nonmethyl[,-1]
  beta <- ( methyl_chip_data$methyl[,-1] + 1 ) / ( sum.counts + 2 )
  beta[ sum.counts < cutoff ] <- NA
  
  # convert the beta matrix to data frames that can be used in the boxplot function
  M.all <- data.frame( pos = as.factor( rep(1:ncol(beta),nrow(beta)) ), beta = as.vector( t( as.matrix( beta * x.C ) ) ) )
  MpG <-   data.frame( pos = as.factor( rep(1:ncol(beta),nrow(beta)) ), beta = as.vector( t( as.matrix( beta * x.C * x.CpG ) ) ) )
  MpH <-   data.frame( pos = as.factor( rep(1:ncol(beta),nrow(beta)) ), beta = as.vector( t( as.matrix( beta * x.C * x.nonCpG ) ) ) )
  W.all <- data.frame( pos = as.factor( rep(1:ncol(beta),nrow(beta)) ), beta = as.vector( t( as.matrix( beta * x.G ) ) ) )
  CpW <-   data.frame( pos = as.factor( rep(1:ncol(beta),nrow(beta)) ), beta = as.vector( t( as.matrix( beta * x.G * x.CpG ) ) ) )
  HpW <-   data.frame( pos = as.factor( rep(1:ncol(beta),nrow(beta)) ), beta = as.vector( t( as.matrix( beta * x.G * x.nonCpG ) ) ) )
  
  # write the boxplots of beta values in a PDF file
  # pdf(file=paste0(outdir,"/methylation_stats.pdf"),width=1500/100,height=800/100)
  y_lim <- c(0, 1)
  
  
  pdf(file=paste0(outdir,"/", experiment, "_methylation_stats.pdf"),width=15,height=8)
  
  par( mfrow=c(2,3) )
  boxplot( beta~pos, data = M.all, ylab="Methylation (forward strand)", main= "All", 
           outcol=rgb(0.5,0.5,0.5), ylim = y_lim )
  
  boxplot( beta~pos, data = MpG, main= "CpG", outcol=rgb(0.5,0.5,0.5), ylim = y_lim )
  
  boxplot( beta~pos, data = MpH, main= "non-CpG", outcol=rgb(0.5,0.5,0.5), ylim = y_lim ) 
  boxplot( beta~pos, data = W.all, xlab = "Position", ylab="Methylation (reverse strand)", 
           outcol=rgb(0.5,0.5,0.5), ylim = y_lim )
  
  boxplot( beta~pos, data = CpW, xlab = "Position", outcol=rgb(0.5,0.5,0.5), ylim = y_lim)
  boxplot( beta~pos, data = HpW, xlab = "Position", outcol=rgb(0.5,0.5,0.5), ylim = y_lim)
  
  null_message <- dev.off()
  # close the PDF
  
  # calculate the correlation between beta at each position and binding strength
  correls.M.all <- apply(
    ( beta * x.C ), 2,
    function(x) { y <- cor.test(c(x,rnorm(5,0,0.0001)),c(log_density,rnorm(5,0,0.0001)), method="pearson", conf.level=0.95); return(c(y$estimate,y$conf.int[1],y$conf.int[2])) } )
  correls.MpG <- apply(
    ( beta * x.C * x.CpG ), 2,
    function(x) { y <- cor.test(c(x,rnorm(5,0,0.0001)),c(log_density,rnorm(5,0,0.0001)), method="pearson", conf.level=0.95); return(c(y$estimate,y$conf.int[1],y$conf.int[2])) } )
  correls.MpH <- apply(
    ( beta * x.C * x.nonCpG ), 2,
    function(x) { y <- cor.test(c(x,rnorm(5,0,0.0001)),c(log_density,rnorm(5,0,0.0001)), method="pearson", conf.level=0.95); return(c(y$estimate,y$conf.int[1],y$conf.int[2])) } )
  correls.W.all <- apply(
    ( beta * x.G ), 2,
    function(x) { y <- cor.test(c(x,rnorm(5,0,0.0001)),c(log_density,rnorm(5,0,0.0001)), method="pearson", conf.level=0.95); return(c(y$estimate,y$conf.int[1],y$conf.int[2])) } )
  correls.CpW <- apply(
    ( beta * x.G * x.CpG ), 2,
    function(x) { y <- cor.test(c(x,rnorm(5,0,0.0001)),c(log_density,rnorm(5,0,0.0001)), method="pearson", conf.level=0.95); return(c(y$estimate,y$conf.int[1],y$conf.int[2])) } )
  correls.HpW <- apply(
    ( beta * x.G * x.nonCpG ), 2,
    function(x) { y <- cor.test(c(x,rnorm(5,0,0.0001)),c(log_density,rnorm(5,0,0.0001)), method="pearson", conf.level=0.95); return(c(y$estimate,y$conf.int[1],y$conf.int[2])) } )
  
  # define a function for plotting the correlations and their confidence intervals
  plot.correlation <- function( correls, xlab, ylab, main ) {
    
    # correls <- correls.CpW 
    plot( correls[1,] ,type="l" , ylim=range(c(correls,0)), ylab=ylab, xlab=xlab, main=main )
    
    suppressWarnings( arrows(1:ncol(correls),correls[2,],1:ncol(correls),correls[3,], 
                             length=0.05, angle=90, code=3) )
    suppressWarnings( lines(c(1,ncol(correls)),c(0,0),lty=2) )
    
  }
  
  # write the correlation plots in a PDF file
  # pdf(file=paste0(outdir,"/methylation_correlation_with_binding.pdf"),width=1500/100,height=800/100)
  pdf(file=paste0(outdir,"/", experiment, "_methylation_correlation_with_binding.pdf"),
      width = 15, height = 8)
  
  par( mfrow=c(2,3) )
  plot.correlation(correls.M.all,"","Correlation (forward strand)","All methylation")
  plot.correlation(correls.MpG,"","","CpG methylation")
  plot.correlation(correls.MpH,"","","Non-CpG methylation")
  plot.correlation(correls.W.all,"Position","Correlation (reverse strand)","")
  plot.correlation(correls.CpW,"Position","","")
  plot.correlation(correls.HpW,"Position","","")
  null_message <- dev.off()
  # close the PDF
}

############################################################################### #
motif.logo <- function( motif ) {
  
  # draw the motif logo in a PDF file
  # calculate the position of the top (or bottom) of Cs and Gs in the motif logo
  C.end <- motif["C",] + 
    (sign(motif["C",]) == sign(motif["A",]) & abs(motif["C",]) > abs(motif["A",] )) * motif["A",] +
    (sign(motif["C",]) == sign(motif["G",]) & abs(motif["C",]) > abs(motif["G",] )) * motif["G",] +
    (sign(motif["C",]) == sign(motif["T",]) & abs(motif["C",]) > abs(motif["T",] )) * motif["T",]
  
  
  G.end <- motif[3,] +
    (sign(motif["G",]) == sign(motif["A",]) & abs(motif["G",]) > abs(motif["A",] )) * motif["A",] +
    (sign(motif["G",]) == sign(motif["C",]) & abs(motif["G",]) > abs(motif["C",] )) * motif["C",] +
    (sign(motif["G",]) == sign(motif["T",]) & abs(motif["G",]) > abs(motif["T",] )) * motif["T",]
  
  
  p <- ggseqlogo(motif[1:4,], method='custom', seq_type='dna') + 
    ylab('Regression coefficient') +
    geom_segment( aes( x = ( 1:ncol( motif ) - 0.1)[ !is.na( motif[5,]) ], 
                       y = C.end[ !is.na( motif[5, ] ) ],
                       xend = ( 1:ncol( motif ) - 0.1 )[ !is.na( motif[5,] ) ], 
                       yend = ( C.end + motif[5,])[ !is.na( motif[5,] ) ] ), 
                  arrow = arrow( angle = 10, length = unit( 0.15, "inches" ) ),
                  color = rgb( 0.5, 0.5, 1) ) +
    
    geom_segment( aes( x = ( 1:ncol( motif ) + 0.1 )[ !is.na( motif[6,] ) ], 
                       y = G.end[ !is.na( motif[6,] ) ],
                       xend = ( 1:ncol( motif ) + 0.1)[ !is.na(motif[ 6, ] ) ], 
                       yend = ( G.end + motif[6, ])[ !is.na( motif[6,] ) ] ), 
                  arrow = arrow( angle = 10, length = unit( 0.15, "inches" ) ), 
                  color = rgb( 0.7, 0.7, 0) )

  return(p) 
}

############################################################################### #
motif.heatmap <- function( motif = motif, coefs_bases = coefs_bases ) {
  
  # motif <- matrix( coefs[,1], nrow = 5, byrow = T) # create the motif matrix based on regression coefficients
  # motif <- rbind( rep.int(0,ncol(motif)), motif ) # add a row in the begining, corresponding to A
  # motif[1:4,] <- scale(motif[1:4,],center = T, scale = F)
  # rownames(motif) <- c("A","C","G","T","CtoM","GtoW")
  # colnames(motif) <- 1:ncol(motif)
  motif[is.na(motif)] <- 0
  
  # Create a matrix of FDR values for filtering out insignificant coefficients
  p.mx <- matrix( coefs_bases[,4], nrow = 5, byrow = T) 
  p.mx[is.na(p.mx)] <- 1
  
  tbl <- data.frame( base = rep( c( 1, 2, 3, 4), each = ncol( motif ) ),
                     pos = rep( 1:ncol( motif ), 4),
                     base_effect_size = c( motif[ 1, ], motif[ 2, ], motif[ 3, ], motif[ 4, ] ),
                     methyl_effectSize_forward = c( rep( 0, times = ncol( motif ) ),
                                                    motif[5, ],
                                                    rep( 0, times = ncol( motif ) * 2 ) ),
                     methyl_logP_forward = c( rep( 0, times = ncol( motif ) ),
                                              -sign( motif[5, ] ) * log10( p.mx[4, ] + 1e-20 ),
                                              rep( 0, times = ncol( motif ) * 2 ) ),
                     methyl_effectSize_reverse = c( rep( 0, times = ncol( motif ) * 2 ),
                                                    motif[6, ],
                                                    rep( 0, times = ncol( motif ) ) ),
                     methyl_logP_reverse = c( rep( 0, times = ncol( motif ) * 2 ),
                                              - sign( motif[6, ]) * log10( p.mx[ 5, ] + 1e-20),
                                              rep( 0, times = ncol( motif ) ) )
                     )
  
  tbl.f <- tbl
  tbl.f$base <- 5 + tbl$base
  tbl.f$methyl_logP <- tbl$methyl_logP_forward
  
  tbl.r <- tbl
  tbl.r$base <- 5 - tbl$base
  tbl.r$methyl_logP <- tbl$methyl_logP_reverse
  
  tbl <- rbind( tbl.f, tbl.r )
  area_scale <- max( abs( tbl$base_effect_size ) )
  tbl$r <- sqrt( abs( tbl$base_effect_size ) / area_scale ) / 2.5
  
  p <- ggplot(tbl) + ## global aes
        geom_tile(aes(x=pos,y=base,fill = methyl_logP),width=1,height=1) + ## to get the rect filled
        geom_circle(aes(x0=pos,y0=base,r=r,
                        fill = sign(base_effect_size)*20 ),
                    color = "white", size=0.5 )  +    ## geom_point for circle illusion
    scale_fill_gradient2( limits=c(-20,20), low = "blue", mid = "white", high = "red" ) +
    geom_hline(yintercept=5)+
    # scale_color_gradient(low = "white",  
    #                       high = "white")+       ## color of the corresponding aes
    scale_y_continuous(breaks=c(1:4,6:9),labels=c("A","C","G","T","A","C","G","T"))+
    scale_x_continuous(breaks=1:ncol(motif))+
    labs(x="Position",fill="Signed log10(P-value) for methylation")+
    coord_fixed() +
    theme_test() +
    theme(legend.position="bottom")
  
  return(p)
}

############################################################################### #
motif.logo2 <- function( coefs ) {
  
  # motif <- matrix( coefs[,1], nrow = 5, byrow = T) # create the motif matrix based on regression coefficients
  # motif <- rbind( rep.int(0,ncol(motif)), motif ) # add a row in the begining, corresponding to A
  # motif[1:4,] <- scale(motif[1:4,],center = T, scale = F)
  # rownames(motif) <- c("A","C","G","T","CtoM","GtoW")
  # colnames(motif) <- 1:ncol(motif)
  
  motif[is.na(motif)] <- 0
  # create a matrix of FDR values for filtering out insignificant coefficients
  # fdr.mx <- matrix( coefs[-1,5], nrow = 5, byrow = T) 
  fdr.mx[is.na(fdr.mx)] <- 1
  
  fdr_max <- 0.005
  fdr_min <- 1e-10
  col_range_max <- 0.05
  col_range_min <- 1e-10
  
  motif[5,fdr.mx[4,]>=fdr_max] <- 0
  motif[6,fdr.mx[5,]>=fdr_max] <- 0
  
  fdr.mx[fdr.mx>fdr_max] <- col_range_max
  fdr.mx[fdr.mx<fdr_min] <- col_range_min
  
  # draw the motif logo in a PDF file
  # calculate the position of the top (or bottom) of Cs and Gs in the motif logo
  C.end <- motif[2,] +
    ( sign(motif[2,])==sign(motif[1,]) & abs(motif[2,])>abs(motif[1,]) ) * motif[1,] +
    ( sign(motif[2,])==sign(motif[3,]) & abs(motif[2,])>abs(motif[3,]) ) * motif[3,] +
    ( sign(motif[2,])==sign(motif[4,]) & abs(motif[2,])>abs(motif[4,]) ) * motif[4,]
  G.end <- motif[3,] +
    ( sign(motif[3,])==sign(motif[1,]) & abs(motif[3,])>abs(motif[1,]) ) * motif[1,] +
    ( sign(motif[3,])==sign(motif[2,]) & abs(motif[3,])>abs(motif[2,]) ) * motif[2,] +
    ( sign(motif[3,])==sign(motif[4,]) & abs(motif[3,])>abs(motif[4,]) ) * motif[4,]
  
  p <- ggseqlogo(motif[1:4,], method='custom', seq_type='dna') +
    ylab('Regression coefficient') +
    geom_segment( aes(
      x=(1:ncol(motif)-0.1)[!is.na(motif[5,])], y=C.end[!is.na(motif[5,])],
      xend=(1:ncol(motif)-0.1)[!is.na(motif[5,])], yend=(C.end+motif[5,])[!is.na(motif[5,])],
      alpha=-log10(fdr.mx[4,]) ),
      arrow=arrow(angle=10,length=unit(0.15,"inches") ), color = rgb(0.5,0.5,1) ) +
    geom_segment( aes(
      x=(1:ncol(motif)+0.1)[!is.na(motif[6,])], y=G.end[!is.na(motif[6,])],
      xend=(1:ncol(motif)+0.1)[!is.na(motif[6,])], yend=(G.end+motif[6,])[!is.na(motif[6,])],
      alpha=-log10(fdr.mx[5,]) ),
      arrow=arrow(angle=10,length=unit(0.15,"inches") ), color = rgb(0.7,0.7,0) ) +
    scale_alpha(limits=c(-log10(col_range_max),-log10(col_range_min)),range=c(0,1)) +
    labs( alpha = "Log10(FDR)" )
  
  return(p)
}

############################################################################### #
draw.motif <- function( coefs, motif, outdir, label ) {

  # write the motif in a file
  write.table( motif, file=paste0(outdir,"/",label,".motif_scaled.txt"),
               sep="\t", row.names = T, col.names = NA, quote = F )
  
  p1 <- motif.logo( motif )
  # p2 <- motif.heatmap( motif = motif, coefs_bases = coefs ) 
  
  ggsave( filename = paste0(outdir,"/",label,"_motif_logo.pdf"),
          plot = grid.arrange(p1, nrow = 1),
          width=25*ncol(motif)/100, height=700/100 )
  
  # ggsave( filename = paste0(outdir,"/",label,"_motif_logo.pdf"),
  #         plot = grid.arrange(p1, p2, nrow = 2),
  #         width=25*ncol(motif)/100, height=700/100 )
  
}

draw.motif.original <- function( coefs, motif, outdir, label ) {
  
  # write the motif in a file
  write.table( motif, file=paste0(outdir,"/",label,".motif_scaled.txt"),
               sep="\t", row.names = T, col.names = NA, quote = F )
  
  p1 <- motif.logo( motif )
  p2 <- motif.heatmap( motif = motif, coefs_bases = coefs ) 
  
  ggsave( filename = paste0(outdir,"/",label,"_motif_logo.pdf"),
          plot = grid.arrange(p1, p2, nrow = 2),
          width=25*ncol(motif)/100, height=700/100 )
  
  # ggsave( filename = paste0(outdir,"/",label,"_motif_logo.pdf"),
  #         plot = grid.arrange(p1, p2, nrow = 2),
  #         width=25*ncol(motif)/100, height=700/100 )
  
}



############################################################################### #

write.sequence.model.av.met <- function( seq_fit, outdir, label ) {
  ## Create the motif from the regression coefficients
  #  Store the coefficients and their associated statistics

  # seq_fit <- fit_CpG_only
  # coefs <- seq_fit$coefficients ## ahcorcha
  coefs <- seq_fit$coefs_summary
  
  ## In case there are missing values in the coefs matrix, add them back
  coefs <- fix.coefs( coefs, seq_fit$fit$coefficients )
  
  ## Calculate the FDR for each coefficient
  coefs <- cbind( coefs, p.adjust( coefs[,4], method = "fdr" ) ) 
  colnames(coefs)[5] <- "FDR" # add the FDR to the matrix
  
  
  write.table( x =  coefs, file = paste0( outdir, "/", label, "_coefficients_with_FDR.txt"),
               quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE )
  
  
  
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
  # intercept <- coefs[1,1] + mean(motif[1:4,])*ncol(motif)
  
  
  # label_control <- paste0( label, "_control" )
  # draw.motif.original.av.meth( coefs_control_meth, motif_control, outdir, label_control )
  # 
  # label_pulldown <- paste0( label, "_pulldown" )
  # draw.motif.original.av.meth( coefs_pulldown_meth, motif_pulldown, outdir, label_pulldown )

  ## Show control and pulldown motif in the same figure
  draw.ctrl.pdwn.motif.av.meth( pulldown_coefs = coefs_pulldown_meth,
                                control_coefs = coefs_control_meth,
                                pulldown_motif = motif_pulldown,
                                control_motif = motif_control, 
                                outdir = outdir,
                                label = label )
  
    
  
  return( list( motif_control = motif_control, 
                motif_pulldown = motif_pulldown ) )
  
}

draw.motif.original.av.meth <- function( coefs, motif, outdir, label ) {
  
  # write the motif in a file
  write.table( motif, file=paste0(outdir,"/",label,".motif_scaled.txt"),
               sep="\t", row.names = T, col.names = NA, quote = F )
  
  p1 <- motif.logo.av.met( motif )
  p2 <- motif.heatmap.av.met( motif = motif, coefs_bases = coefs ) 
  
  this_width <- max( 2 +  0.85*ncol( motif ), 1 )
  
  ggsave( filename = paste0(outdir,"/",label,"_motif_logo.pdf"),
          plot = p1 / p2, 
          units = "cm", device = "pdf",
          width = this_width, height = 23 ) 
          #  width = 25 * ncol( motif ) / 100, height = 700 / 100 ) 

}


draw.ctrl.pdwn.motif.av.meth <- function( pulldown_coefs = coefs_pulldown_meth,
                                          control_coefs = coefs_control_meth,
                                          pulldown_motif = motif_pulldown,
                                          control_motif = motif_control, 
                                          outdir = outdir,
                                          label = "label" ) {

  # pulldown_coefs = coefs_pulldown_meth
  # control_coefs = coefs_control_meth
  # pulldown_motif = motif_pulldown
  # control_motif = motif_control
  
  # write the control motif in a file
  write.table( control_motif, file=paste0(outdir,"/",label,"_control.motif_scaled.txt"),
               sep="\t", row.names = T, col.names = NA, quote = F )
  # write the pulldown motif in a file
  write.table( pulldown_motif, file=paste0(outdir,"/",label,"_pulldown.motif_scaled.txt"),
               sep="\t", row.names = T, col.names = NA, quote = F )
                      

  
  # area_scale <- max( abs( motif_pulldown[1:4, ] ) )
  # area_scale <- max( abs( control_motif[1:4, ] ) )
  area_scale <- max( max( abs( pulldown_motif[1:4, ] ) ), max( abs( control_motif[1:4, ] ) ) )
  area_scale_ctrl <- max( abs( control_motif[1:4, ] ) )
  
  
  this_width <- max( 2 +  0.85*ncol( pulldown_motif ), 8 )

  
      
  pdwn_motif <- motif.logo.av.met( motif = pulldown_motif, coefs_bases = pulldown_coefs )
  pdwn_htmp <- motif.heatmap.av.met( motif = pulldown_motif, coefs_bases = pulldown_coefs,
                                     area_scale = area_scale )
  
  
  ggsave( filename = paste0(outdir,"/",label,"_pulldown_motif_logo.pdf"), 
          plot = pdwn_motif / pdwn_htmp, 
          units = "cm", device = "pdf", width = this_width, height = 23 )
  
  
  
  ctrl_motif <- motif.logo.av.met( motif = control_motif, coefs_bases = control_coefs )
  ctrl_htmp <- motif.heatmap.av.met( motif = control_motif, coefs_bases = control_coefs,
                                     area_scale = area_scale_ctrl )
  
  ggsave( filename = paste0(outdir,"/",label,"_control_motif_logo.pdf"), plot = ctrl_motif / ctrl_htmp,
          units = "cm", device = "pdf", width = this_width, height = 23 )
  
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
  
  
  pp <- ( ( pdwn_motif / pdwn_htmp ) | ( ctrl_motif  / ctrl_htmp ) ) +
            plot_annotation( theme( title = ""  ) )
  
  
  
  
  ggsave( filename = paste0(outdir,"/",label,"_control_n_pulldown_motif_logo.pdf"), plot = pp,
          units = "cm", device = "pdf", width = this_width*2, height = 23 )

  }

motif.logo.av.met <- function( motif, coefs_bases ) {
  
  # motif <- motif_pulldown
  # coefs_bases <- pulldown_coefs
  
   
  # draw the motif logo in a PDF file
  # calculate the position of the top (or bottom) of Cs and Gs in the motif logo
  C.end <- motif["C",] + 
           ( sign(motif["C",]) == sign(motif["A",]) & abs(motif["C",]) > abs(motif["A",] )) * motif["A",] +
           ( sign(motif["C",]) == sign(motif["G",]) & abs(motif["C",]) > abs(motif["G",] )) * motif["G",] +
           ( sign(motif["C",]) == sign(motif["T",]) & abs(motif["C",]) > abs(motif["T",] )) * motif["T",]

  G.end <- motif[3,] +
            ( sign(motif["G",]) == sign(motif["A",]) & abs(motif["G",]) > abs(motif["A",] )) * motif["A",] +
            ( sign(motif["G",]) == sign(motif["C",]) & abs(motif["G",]) > abs(motif["C",] )) * motif["C",] +
            ( sign(motif["G",]) == sign(motif["T",]) & abs(motif["G",]) > abs(motif["T",] )) * motif["T",]
  
  
  arrow_pval <- - log10( coefs_bases[ !is.na( coefs_bases$Estimate ), "Pr...z.." ] + 1e-20 )
  arrow_pval <- arrow_pval / 20
  
  
  ############################################################################################################## ##
  
  
  p <- ggseqlogo( motif[1:4,], method='custom', seq_type='dna') + 
                  ylab('Regression coefficient') +
                  geom_segment( aes( x = ( 1:ncol( motif ) - 0.1 )[ !is.na( motif[5,] ) ],
                                     y = C.end[ !is.na( motif[5, ] ) ],
                                     xend = ( 1:ncol( motif ) - 0.1 )[ !is.na( motif[5,] ) ], 
                                     yend = ( C.end + motif[5,] )[ !is.na( motif[5,] ) ],
                                     ), alpha = arrow_pval,
                                arrow = arrow( angle = 10, length = unit( 0.15, "inches" ) ),
                                color = rgb( 100, 149, 237, maxColorValue = 255 ) ) 
  
  return(p) }

motif.heatmap.av.met <- function( motif = motif, coefs_bases = coefs_bases, area_scale = area_scale ) {
  
  # coefs_bases <- coefs_pulldown_meth
  # motif <- motif_pulldown
  # area_scale <- max( abs( motif_pulldown[1:4, ] ) )
  
  
  motif[is.na(motif)] <- 0
  
  # Create a matrix of FDR values for filtering out insignificant coefficients
  coefs_bases[, "Pr...z.."][is.na( coefs_bases[, "Pr...z.."] )] <- 1
  
  tbl <- data.frame( base = rep( c( 1, 2, 3, 4), each = ncol( motif ) ),
                     pos = rep( 1:ncol( motif ), 4),
                     base_effect_size = c( motif[ 1, ], motif[ 2, ], motif[ 3, ], motif[ 4, ] ),
                     methyl_effectSize_forward = c( rep( 0, times = ncol( motif ) ),
                                                    motif[5, ],
                                                    rep( 0, times = ncol( motif ) * 2 ) ),
                     methyl_logP_forward = c( rep( 0, times = ncol( motif ) ),
                                              -sign( motif[5, ] ) * log10( coefs_bases[, "Pr...z.."] + 1e-20 ),
                                              rep( 0, times = ncol( motif ) * 2 ) )
                     )
  
  # area_scale <- max( abs( tbl$base_effect_size ) )
  
  tbl$r <- sqrt( abs( tbl$base_effect_size ) / area_scale ) / 2.5
  ########################################### ###
  #####   ahcorcha
  
  
  
  ########################################### ###
  
  p <- ggplot( tbl ) +
               # To get the rect filled
               geom_tile( aes( x = pos, y = base, fill = methyl_logP_forward ), width = 1, height = 1 ) + 
               geom_circle( aes( x0 = pos, y0 = base, r = r, fill = sign( base_effect_size ) * 20 ),
                            color = "white", size = 0.5 )  + # geom_point for circle illusion
               scale_fill_gradient2( limits = c(-20, 20), low = "blue", mid = "white", high = "red" ) +
               # scale_color_gradient(low = "white", high = "white") + ## color of the corresponding aes
               scale_y_continuous( breaks = c(1:4), labels = c( "A","C","G","T" ) ) +
               scale_x_continuous( breaks = 1:ncol( motif ) ) +
               labs( x = "Position\n\nSigned log10(P-value) for methylation", 
                     fill = "" ) +
               coord_fixed() + theme_test() + theme( legend.position = "bottom") 
  
  return(p)
}

############################################################################### #



############################################################################### #
merge.models <- function( model1, model2, outdir, label ) {
  
  motif <- ( model1$motif + model2$motif ) / 2 # take the average of the two motifs
  motif[ which(sign(model1$motif) != sign(model2$motif)) ] <- NA # anywhere that the sign of the two motifs are different, replace with NA in the new motif
  motif[1:4,] <- ( model1$motif[1:4,] + model2$motif[1:4,] ) / 2 # return the values of the first four rows
  
  #ggsave( filename = paste0(outdir,"/",label,".motif_logo.pdf"),
  #        plot = motif.logo( motif ),
  #        width=25*ncol(motif)/100, height=400/100 )
  
  ggsave( filename = paste0(outdir,"/",label,".motif_logo.pdf"),
          plot = motif.logo( motif ),
          width=25*ncol(motif)/100, height=400/100 )
}

############################################################################### #
calculate_mean_by_bins <- function(acc = NULL, pfm_length = NULL){
  
  # bin_sizes <- c(988, 269, 73, 20)
  bin_sizes <- head(rev(seq(0,1000, by=200)), -1)
  
  init <- 1
  n <- (length(bin_sizes) * 2) + 1 
  bin_reg <- array(NA, dim=c(n,3))
  
  while(init <= length(bin_sizes) ){
    
    bin_reg[init, 1] <- paste0("bin_up_", length(bin_sizes) - init + 1)
    bin_reg[init, 2] <- - bin_sizes[init]
    bin_reg[init, 3] <- - bin_sizes[init + 1] - 1
    
    bin_reg[(length(bin_reg[,2]) - init), 1] <- paste0("bin_down_", length(bin_sizes) - init + 1)
    bin_reg[(length(bin_reg[,2]) - init), 2] <- pfm_length + bin_sizes[init + 1]
    bin_reg[(length(bin_reg[,2]) - init), 3] <- pfm_length + bin_sizes[init] - 1
    
    if (is.na(bin_sizes[init + 1]) ){
      ## This will only run once.
      bin_reg[init, 3] <- 0
      bin_reg[length(bin_sizes) + 1, 1] <- "motif_bin"
      bin_reg[length(bin_sizes) + 1, 2] <- 0
      bin_reg[length(bin_sizes) + 1, 3] <- pfm_length - 1            
      bin_reg[length(bin_reg[,1]), 1] <- "bin_down_1"
      bin_reg[length(bin_reg[,1]), 2] <- pfm_length
      bin_reg[length(bin_reg[,1]), 3] <- pfm_length + bin_sizes[length(bin_sizes)] - 1
    }
    init <- init + 1
  }
  
  old_col <- colnames(acc[, -1])
  ## Add column with means
  init2 <- 1
  
  while(init2 <= length(bin_reg[,1]) ){
    
    acc[bin_reg[init2, 1]] <- rowMeans( acc[, c(paste0("pos","(",( bin_reg[init2, 2] ):( bin_reg[init2, 3] ),")") ) ]) 
    init2 <- init2 + 1
  }
  
  acc <- acc[ , !(names(acc) %in% old_col)]
  
  # Reorder columns
  df_order <- c(1,2,3,4,5,6,7,12,8,9,10,11)
  acc <- acc[,df_order]
  acc[,-1] <- log(acc[,-1])
  
  return(acc)
}

############################################################################### #
se <- function(x) sqrt(var(x)/length(x))

############################################################################### #
sequence_fit <- function( seq, CpGs, methyl, nonmethyl, dna_acc, 
                          target, pfm_length, flanking = 0, CpG_only = T) {

  # seq <- methyl_chip_data$seq
  # CpGs <- methyl_chip_data$CpGs
  # methyl <- methyl_chip_data$methyl
  # nonmethyl <- methyl_chip_data$nonmethyl
  # dna_acc <- methyl_chip_data$acc
  # target <- methyl_chip_data$target
  # pfm_length <- methyl_chip_data$pfm_length
  # flanking <- opt$flanking
  # CpG_only <- T
  
  ## Preparation
  # extract the sequence matrix and convert it to one-hot encoding
  x.A <- ( seq[,-1] == "A" ) * 1
  x.C <- ( seq[,-1] == "C" ) * 1
  x.G <- ( seq[,-1] == "G" ) * 1
  x.T <- ( seq[,-1] == "T" ) * 1
  rownames(CpGs) <- CpGs$Name
  x.CpG <- CpGs[,-1] * 1
  
  ncolx <- ncol(x.A)
  x.CG <- x.C[,-ncolx] * x.G[,-1]
  
  ## Calculate the fraction of methylated C's in the reverse strand (in front of every G)
  x.W <- x.G * ( methyl[,-1] + 1 ) / ( methyl[,-1] + nonmethyl[,-1] + 2 ) 
  ## Calculate the fraction of methylated C's in the forward strand
  x.M <- x.C * ( methyl[,-1] + 1 ) / ( methyl[,-1] + nonmethyl[,-1] + 2 ) 
  
  ## Methylation in reverse strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.W[ is.na(x.W) ] <- 0
  # only keep the W's that are within a CpG
  if( CpG_only ) { x.W <- x.W * x.CpG }
  ## Methylation in forward strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.M[ is.na(x.M) ] <- 0
  # only keep the W's that are within a CpG
  if( CpG_only ) { x.M <- x.M * x.CpG }
  
  
  ## Methylation average M (forward) W (reverse).
  x.Met <- ( ( x.M[,-ncolx] + x.W[,-1] ) / 2 )
  
  
  motif <- paste0( "pos(", 0:(pfm_length - 1) , ")" )
  upsteam_flank <- paste0( "pos(", (-flanking):-1 , ")" )
  downsteam_flank <- paste0( "pos(", (pfm_length ):(pfm_length + flanking - 1), ")" )

  
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
  acc <- calculate_mean_by_bins(acc = dna_acc, pfm_length = pfm_length)

  ## Tags
  tags <- as.factor(target$pulldown.tag.old)
  tags_ctrl <- as.factor(target$ctrl.tag.old)
  
  
  c <- c( target$ctrl.tag.old, target$pulldown.tag.old )
  c <- as.data.frame( c )
  
  rownames( c ) <- c( paste0("control.", target$Name), paste0( "pulldown.", target$Name ) )
  
  
  t_all <- c( rep( 0, nrow( target ) ), rep( 1, nrow( target ) ) )
  t_all <- as.data.frame( t_all )
  rownames( t_all ) <- c( paste0("control.", target$Name), paste0( "pulldown.", target$Name ) )
  
  
  
  XX <- data.frame( acc = acc[,-1], 
                         x.C = x.C, x.G = x.G, x.T = x.T,
                         x.CG = x.CG, x.Met = x.Met,
                         x.T_up = x.T_up, x.C_up = x.C_up, x.G_up = x.G_up,
                         x.M_up = x.M_up, x.W_up = x.W_up,
                         x.T_down = x.T_down, x.C_down = x.C_down,x.G_down = x.G_down,
                         x.M_down = x.M_down, x.W_down = x.W_down)
  
  rownames( XX ) <- NULL
  XX <- rbind( XX, XX )
  rownames( XX ) <- c( paste0("control.", target$Name), paste0( "pulldown.", target$Name ) )
  
  
  
  ## Cross-validation
  cv <- 10; counter <- 1; set.seed(001) # just to make it reproducible
  
  rows <- rownames( target )
  rows <- sample(rows); models <- list()
  correlations <- vector(mode="double")
  
  # size <- floor( length(rows) / cv )
  parts <- split(rows, cut(seq_along(rows), cv, labels = FALSE))
  
  all_c_predicted <- vector(mode="double")
  all_c_real <- vector(mode="double")
  
  
  for(i in parts){
    
    ################################## ###
    ### Train negbin GLM 
    
    train_names <- unlist(parts[-counter])
    train_names <- c( paste0( "control.", train_names ), paste0( "pulldown.", train_names ) )
    
    XX_train <- as.data.frame( XX[ train_names, ] )
    t <- t_all[ train_names, ]
    t_train <- t_all[ train_names, ]
    c_train <- c[ train_names, ]

    fit <- glm.nb( c_train ~ . + t + .:t, data = XX_train )
    
    
    ################################## ###
    ### Get Observed tags  for test set
    ### log( c_pulldown_observed / c_control_observed )
    
    test_names <- unlist(parts[counter])
    test_names <- c( paste0( "control.", test_names ), paste0( "pulldown.", test_names ) )
    c_observed <- c[ test_names, ]
    c_ctrl_observed <- c_observed[ 1 : ( length( test_names ) / 2 ) ]
    c_pdwn_observed <- c_observed[ ( (length( test_names ) / 2) + 1 ) : ( length( test_names ) ) ]

    ## log( pulldown_observed / control_observed )
    log_c_ratio_observed <- log( c_pdwn_observed / c_ctrl_observed )
    
    
    
    ################################## ###
    ### Define predicted tags
    
    ## Get all coefficients from the GLM (trained set)
    fit_coefficients <- as.data.frame( coefficients( summary( fit ) ) )

    ## Get the pulldown coefficients
    fit_pdwn_coefficients <- fit_coefficients[ grepl( pattern = ":t$", 
                                                      x = rownames( fit_coefficients ) ), ]

    ## Get the name of the variables included in th model    
    X_var_names <- gsub( ":t", "", rownames( fit_pdwn_coefficients ) )
    ## Make sure those variables are in the correct order for the matrix*vector
    XX_test <- as.matrix( XX[ test_names, X_var_names ] )
    
    ## Mult. the PULLDOWN coeffs vector and the variable matrix
    c_pdwn_predicted <- as.vector( XX_test %*% fit_pdwn_coefficients$Estimate ) # 1st must be 1.2476048 in test
    # Only the first half of XX
    c_pdwn_predicted <- c_pdwn_predicted[ 1 : ( length( test_names ) / 2 ) ]
    
    
    log_c_ratio_predicted <- c_pdwn_predicted
    
    
    all_c_predicted <- append(all_c_predicted, log_c_ratio_predicted)
    all_c_real <- append(all_c_real, log_c_ratio_observed )
    
    correlations[counter] <- cor( log_c_ratio_observed, log_c_ratio_predicted )
    counter <- counter + 1    
    coefs_summary <- coefficients( summary( fit ) )
  }
  
  return( list(fit = fit, 
               correlation = correlations, 
               all_c_predicted = all_c_predicted, 
               all_c_real = all_c_real,
               coefs_summary = coefs_summary,
               XX_train = XX_train,
               t_train = t_train,
               c_train = c_train,
               x.Met = x.Met
               ) )
}

############################################################################### #
fix.coefs <- function( coefs, coef.full ) {
  
  # create two data frames, one based on coefs and the other based on the coef.names
  d1 <- data.frame(Name=rownames(coefs),coefs)
  d2 <- data.frame(Name=names(coef.full),Order=1:length(coef.full))
  d3 <- merge(d2,d1,by="Name",all=T)
  d3 <- d3[order(d3$Order),]
  
  new.coefs <- as.matrix(d3[,3:6])
  rownames(new.coefs) <- d3$Name
  
  return( new.coefs )
}

############################################################################### #
write.sequence.model.for.dinucleotide <- function( seq_fit, outdir, label ) {

  ## Create the motif from the regression coefficients
  #  Store the coefficients and their associated statistics
  # seq_fit <- fit_CpG_only
  
  coefs <- seq_fit$coefs_summary
  
  ## In case there are missing values in the coefs matrix, add them back
  coefs <- fix.coefs( coefs, seq_fit$fit$coefficients )
  
  ## Calculate the FDR for each coefficient
  coefs <- cbind( coefs, p.adjust( coefs[,4], method = "fdr" ) ) 
  colnames(coefs)[5] <- "FDR" # add the FDR to the matrix
  
  
  #### Format coefs for dinucleotides ####
  
  # Get base specific coefficients.
  coefs_NN <- coefs[grepl( pattern = "*_trainpos\\(", x = rownames(coefs)), ]
  coefs_NN <- as.data.frame(coefs_NN)
  
  # Get position per base
  positions <- gsub( pattern = ".*_trainpos\\(", replacement = "", x = rownames(coefs_NN) )
  coefs_NN$position <- gsub( pattern = "\\)", replacement = "", x = positions )

  # Get base/bases specific per position
  base <- gsub( pattern = "x.", replacement = "", x = rownames(coefs_NN) )
  coefs_NN$base <- gsub( pattern = "_trainpos\\(.*", replacement = "", x = base )
  
  # Number of nucleotides per position
  coefs_NN$nucleotides <- nchar(coefs_NN$base)
  
  # Take out CpG position coefficients
  coefs_NN <- coefs_NN[ ( coefs_NN$base != "CG" ) | ( coefs_NN$base != "Met" ), ]
  
  # Just get the first base, for example: A out of AG
  coefs_NN$base <- substring( coefs_NN$base, 1, 1 )
  
  # Zero for coefficients with 0  
  coefs_NN$Estimate[ is.na( coefs_NN$Estimate ) ] <- 0
  
  # Divide coefficients/4 only for dinucleotides
  coefs_NN[ coefs_NN$nucleotides == 2, "Estimate"] <-
    coefs_NN[ coefs_NN$nucleotides == 2, "Estimate"] / 4

  #####
  
  
  # filter non-significant W and M coefficients
  coefs_NN[(coefs_NN$base == "W")&(coefs_NN$FDR >= 0.01),]$Estimate <- NA
  coefs_NN[(coefs_NN$base == "M")&(coefs_NN$FDR >= 0.01),]$Estimate <- NA
  
  
  
  #################################################### #
  # View(coefs_NN[ (coefs_NN$position == 6)& (coefs_NN$base == "A") , ])
  #################################################### #
  
  
  
  coefs_motif <- aggregate( x = coefs_NN$Estimate, 
                            by = list( coefs_NN$position, coefs_NN$base ), 
                            FUN = sum )
  
  coefs_motif$Group.1 <- as.numeric( coefs_motif$Group.1 )
  coefs_motif <- coefs_motif[ with( coefs_motif, order( Group.2, Group.1) ),]
  
  
  motif <- matrix( coefs_motif$x, nrow = 6, byrow = T )  
  rownames(motif) <- unique(coefs_motif$Group.2)
    
  
  rownames(motif) <- gsub( "W", "GtoW", x = rownames(motif) )
  rownames(motif) <- gsub( "M", "CtoM", x = rownames(motif) )
  
  motif <- motif[ c(1, 2, 3, 5, 4, 6),  ]
  
  ## Scale the motif so that each column has a mean of zero, and also correct the 
  #  intercept so that the motif scores would still exactly match the regression results
  motif[1:4,] <- scale(motif[1:4,],center = T, scale = F)

    
  # write the coefficient matrix in a file
  write.table( coefs, file = paste0( outdir, "/", label, "_coefficients.txt" ),
               sep = "\t", row.names = T, col.names = NA, quote = F )
  
  draw.motif( coefs_bases, motif, outdir, label )
  
  return( list( motif = motif, coefs = coefs_NN ) )
}

############################################################################### #
write.sequence.model <- function( seq_fit, outdir, label ) {
  ## Create the motif from the regression coefficients
  #  Store the coefficients and their associated statistics

  # seq_fit <- fit_CpG_only
  # coefs <- seq_fit$coefficients
  
  coefs <- seq_fit$coefs_summary
  
  ## In case there are missing values in the coefs matrix, add them back
  coefs <- fix.coefs( coefs, seq_fit$fit$coefficients )
  
  ## Calculate the FDR for each coefficient
  coefs <- cbind( coefs, p.adjust( coefs[,4], method = "fdr" ) ) 
  colnames(coefs)[5] <- "FDR" # add the FDR to the matrix
  
  ## Get coef. for bases and methylation.
  coefs_bases <- coefs[grepl( pattern = "pos\\(", x = rownames(coefs)), ]
  coefs_bases <- data.frame(coefs_bases)
  
  ## Remove CpG
  # Get base/bases specific per position
  base <- gsub( pattern = "x.", replacement = "", x = rownames(coefs_bases) )
  coefs_bases$base <- gsub( pattern = "_trainpos\\(.*", replacement = "", x = base )

  # Take out CpG position coefficients
  coefs_bases <- coefs_bases[ coefs_bases$base != "CG", ]
  
  ## Create the motif matrix based on regression coefficients
  motif <- matrix( coefs_bases[, 1], nrow = 5, byrow = T )   
  
  ## Add a row in the begining, corresponding to A
  motif <- rbind( rep.int(0,ncol(motif)), motif )
  rownames(motif) <- c( "A", "C", "G", "T", "CtoM", "GtoW" )
  colnames(motif) <- 1:ncol(motif)
  
  ## Process the CtoM and GtoW coefficients
  # Create a matrix of FDR values for filtering out insignificant coefficients
  fdr.mx <- matrix( coefs_bases[,5], nrow = 5, byrow = T) 
  
  motif[5,fdr.mx[4,]>=0.01] <- NA # filter non-significant CtoM coefficients
  motif[6,fdr.mx[5,]>=0.01] <- NA # filter non-significant GtoW coefficients
  
  ## Scale the motif so that each column has a mean of zero, and also correct the 
  #  intercept so that the motif scores would still exactly match the regression results
  intercept <- coefs[1,1] + mean(motif[1:4,])*ncol(motif)
  motif[1:4,] <- scale(motif[1:4,],center = T, scale = F)
  
  # write the coefficient matrix in a file
  write.table( coefs, file = paste0( outdir, "/", label, "_coefficients.txt" ),
               sep = "\t", row.names = T, col.names = NA, quote = F )
  
  draw.motif.original( coefs_bases, motif, outdir, label )
  
  return( list( motif = motif, coefs = coefs_bases ) )
}

############################################################################### #
write.sequence.model.from.text <- function( coefs, outdir, label ) {
  
  ## Create the motif matrix based on regression coefficients
  motif <- matrix( coefs[-1, 1], nrow = 5, byrow = T )   
  
  ## Add a row in the begining, corresponding to A
  motif <- rbind( rep.int(0,ncol(motif)), motif )
  rownames(motif) <- c( "A", "C", "G", "T", "CtoM", "GtoW" )
  colnames(motif) <- 1:ncol(motif)
  
  # Process the CtoM and GtoW coefficients
  coefs <- cbind( coefs, p.adjust( coefs[,4], method = "fdr" ) ) # calculate the FDR for each coefficient
  colnames(coefs)[5] <- "FDR" # add the FDR to the matrix
  
  # Create a matrix of FDR values for filtering out insignificant coefficients
  fdr.mx <- matrix( coefs[-1,5], nrow = 5, byrow = T) 
  
  motif[5,fdr.mx[4,]>=0.01] <- NA # filter non-significant CtoM coefficients
  motif[6,fdr.mx[5,]>=0.01] <- NA # filter non-significant GtoW coefficients
  
  ## Scale the motif so that each column has a mean of zero, and also correct the 
  #  intercept so that the motif scores would still exactly match the regression results
  intercept <- coefs[1,1] + mean(motif[1:4,])*ncol(motif)
  motif[1:4,] <- scale(motif[1:4,],center = T, scale = F)
  
  # write the coefficient matrix in a file
  write.table( coefs, file=paste0(outdir,"/",label,".coefficients.txt"),
               sep="\t", row.names = T, col.names = NA, quote = F )
  
  draw.motif( coefs, motif, outdir, label )
  
  return( list( motif=motif,coefs=coefs ) )
}

############################################################################### #
plot_dna_acc_coefficients <- function(fit_model, plot_name){
  
  # fit_model <- fit_CpG_only$fit
  # plot_name <- dna_acc_plot_name
  
  ## Subset coefficients
  coefs <- as.data.frame( coef( summary( fit_model ) ) )
  coefs$FDR <- as.double( p.adjust( coefs[,4], method = "fdr" ) )
  coefs$FDR[ coefs$FDR >= 0.1 ] <- NA
  
  coefs$FDR_color[ is.na( coefs$FDR ) ] <- "< 0.1"
  coefs$FDR_color[ ! is.na( coefs$FDR ) ] <- "> 0.1"
  
  dna_coefs <- coefs[ grepl( pattern = "acc", x = rownames( coefs ) ), ]
  
  dna_coefs_ctrl <- dna_coefs[ grepl( pattern = ":t", x = rownames( dna_coefs ) ), ]
  dna_coefs_ctrl$Name <- rownames( dna_coefs_ctrl )
  dna_coefs_ctrl$Name <- gsub( "_", "\\.", dna_coefs_ctrl$Name )
  
  dna_coefs_pdwn <- dna_coefs[ ! grepl( pattern = ":t", x = rownames( dna_coefs ) ), ]
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
  
  
  dna_pdwn_plot <- ggplot( data = dna_coefs_pdwn, aes( x = Name, y = Estimate, color = FDR_color ) ) +
                           geom_point( size = 1 ) +
                           geom_errorbar( aes( x = Name, width=0.5,
                                               ymin = Estimate - `Std. Error`, 
                                               ymax = Estimate + `Std. Error` ) ) +
                           theme_minimal() + labs( x = "Pulldown", y = "Estimate" ) +
                           theme( axis.text.x = element_text( angle = 290, hjust = 0 ),
                                  legend.position = "none" ) + 
                           ylim( y_min, y_max ) + 
                           scale_color_manual( values = c( "grey", "black" ) )
  
  
  dna_ctrl_plot <- ggplot( data = dna_coefs_ctrl, aes( x = Name, y = Estimate, color = FDR_color ) ) +
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
  
  ggsave( filename = plot_name, plot = coefficient_plot, 
          units = "cm", width = 20, height = 10 ) 
  
  
  }

############################################################################### #
predict_new_data <- function(fit, cell_line_data, exp_name, flanking, 
                             CpG_only = T) {
  
  # fit <- fit_hek293
  # cell_line_data <- Cell_line_1
  
  seq <- cell_line_data$seq
  CpGs <- cell_line_data$CpGs
  methyl <- cell_line_data$methyl
  nonmethyl <- cell_line_data$nonmethyl
  dna_acc <- cell_line_data$acc
  target <- cell_line_data$target
  pfm_length <- cell_line_data$pfm_length
  CpG_only <- T
  flanking <- 20
  
  
  ## Preparation
  # extract the sequence matrix and convert it to one-hot encoding
  x.A <- ( seq[,-1] == "A" ) * 1
  x.C <- ( seq[,-1] == "C" ) * 1
  x.G <- ( seq[,-1] == "G" ) * 1
  x.T <- ( seq[,-1] == "T" ) * 1
  rownames(CpGs) <- CpGs$Name
  x.CpG <- CpGs[,-1] * 1
  
  ncolx <- ncol(x.A)
  x.CG <- x.C[,-ncolx] * x.G[,-1]
  
  ## Calculate the fraction of methylated C's in the reverse strand (in front of every G)
  x.W <- x.G * ( methyl[,-1] + 1 ) / ( methyl[,-1] + nonmethyl[,-1] + 2 ) 
  ## Calculate the fraction of methylated C's in the forward strand
  x.M <- x.C * ( methyl[,-1] + 1 ) / ( methyl[,-1] + nonmethyl[,-1] + 2 ) 
  
  ## Methylation in reverse strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.W[ is.na(x.W) ] <- 0
  # only keep the W's that are within a CpG
  if( CpG_only ) { x.W <- x.W * x.CpG }
  ## Methylation in forward strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.M[ is.na(x.M) ] <- 0
  # only keep the W's that are within a CpG
  if( CpG_only ) { x.M <- x.M * x.CpG }
  
  
  ## Methylation average M (forward) W (reverse).
  x.Met <- ( ( x.M[,-ncolx] + x.W[,-1] ) / 2 )
  
  
  motif <- paste0( "pos(", 0:(pfm_length - 1) , ")" )
  upsteam_flank <- paste0( "pos(", (-flanking):-1 , ")" )
  downsteam_flank <- paste0( "pos(", (pfm_length ):(pfm_length + flanking - 1), ")" )

  
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
  acc <- calculate_mean_by_bins(acc = dna_acc, pfm_length = pfm_length)

  ## Tags
  tags <- as.factor(target$pulldown.tag.old)
  tags_ctrl <- as.factor(target$ctrl.tag.old)
  
  
  # c <- c( target$ctrl.tag.old, target$pulldown.tag.old )
  c <- c( target$pulldown.tag.old )
  c <- as.data.frame( c )
  # rownames( c ) <- c( paste0("control.", target$Name), paste0( "pulldown.", target$Name ) )
  # rownames( c ) <- c( paste0( "pulldown.", target$Name ) )
  rownames( c ) <- target$Name
  
  
  # t <- c( rep( 0, nrow( target ) ), rep( 1, nrow( target ) ) )
  t <- c( rep( 1, nrow( target ) ) )
  
  t <- as.data.frame( t )
  # rownames( t ) <- c( paste0("control.", target$Name), paste0( "pulldown.", target$Name ) )
  # rownames( t ) <- c( paste0( "pulldown.", target$Name ) )
  rownames( t ) <- target$Name
  
  
  
  XX <- data.frame( acc = acc[,-1], 
                    x.C = x.C, x.G = x.G, x.T = x.T,
                    x.CG = x.CG, x.Met = x.Met,
                    x.T_up = x.T_up, x.C_up = x.C_up, x.G_up = x.G_up,
                    x.M_up = x.M_up, x.W_up = x.W_up,
                    x.T_down = x.T_down, x.C_down = x.C_down,x.G_down = x.G_down,
                    x.M_down = x.M_down, x.W_down = x.W_down )
  
  rownames( XX ) <- NULL
  # XX <- rbind( XX, XX )
  #rownames( XX ) <- c( paste0("control.", target$Name), paste0( "pulldown.", target$Name ) )
  # rownames( XX ) <- c( paste0( "pulldown.", target$Name ) )
  rownames( XX ) <- target$Name
  
  a <- cbind(XX, c, t)
  
  
  tags_predicted <- predict( fit, a, type = "response" )
  
  return(tags_predicted)
}


predict_new_cell_line <- function( seq = seq, CpGs = CpGs, methyl = methyl, nonmethyl = nonmethyl,
                                   dna_acc = dna_acc, target = target, pfm_length = pfm_length, 
                                   flanking = 0, fit_object = fit_ctcf_hek293,
                                   CpG_only = T) {

  # seq <- methyl_chip_data$seq
  # CpGs <- methyl_chip_data$CpGs
  # methyl <- methyl_chip_data$methyl
  # nonmethyl <- methyl_chip_data$nonmethyl
  # dna_acc <- methyl_chip_data$acc
  # target <- methyl_chip_data$target
  # pfm_length <- methyl_chip_data$pfm_length
  # flanking <- opt$flanking
  # CpG_only <- T
  # fit_object = fit_ctcf_hek293
  
  ## Preparation
  # extract the sequence matrix and convert it to one-hot encoding
  x.A <- ( seq[,-1] == "A" ) * 1
  x.C <- ( seq[,-1] == "C" ) * 1
  x.G <- ( seq[,-1] == "G" ) * 1
  x.T <- ( seq[,-1] == "T" ) * 1
  rownames(CpGs) <- CpGs$Name
  x.CpG <- CpGs[,-1] * 1
  
  ncolx <- ncol(x.A)
  x.CG <- x.C[,-ncolx] * x.G[,-1]
  
  ## Calculate the fraction of methylated C's in the reverse strand (in front of every G)
  x.W <- x.G * ( methyl[,-1] + 1 ) / ( methyl[,-1] + nonmethyl[,-1] + 2 ) 
  ## Calculate the fraction of methylated C's in the forward strand
  x.M <- x.C * ( methyl[,-1] + 1 ) / ( methyl[,-1] + nonmethyl[,-1] + 2 ) 
  
  ## Methylation in reverse strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.W[ is.na(x.W) ] <- 0
  # only keep the W's that are within a CpG
  if( CpG_only ) { x.W <- x.W * x.CpG }
  ## Methylation in forward strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.M[ is.na(x.M) ] <- 0
  # only keep the W's that are within a CpG
  if( CpG_only ) { x.M <- x.M * x.CpG }
  
  
  ## Methylation average M (forward) W (reverse).
  x.Met <- ( ( x.M[,-ncolx] + x.W[,-1] ) / 2 )
  
  
  motif <- paste0( "pos(", 0:(pfm_length - 1) , ")" )
  upsteam_flank <- paste0( "pos(", (-flanking):-1 , ")" )
  downsteam_flank <- paste0( "pos(", (pfm_length ):(pfm_length + flanking - 1), ")" )

  
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
  acc <- calculate_mean_by_bins(acc = dna_acc, pfm_length = pfm_length)

  
  X <- data.frame( acc = acc[,-1], 
                    x.C = x.C, x.G = x.G, x.T = x.T,
                    x.CG = x.CG, x.Met = x.Met,
                    x.T_up = x.T_up, x.C_up = x.C_up, x.G_up = x.G_up,
                    x.M_up = x.M_up, x.W_up = x.W_up,
                    x.T_down = x.T_down, x.C_down = x.C_down,x.G_down = x.G_down,
                    x.M_down = x.M_down, x.W_down = x.W_down )
  
  
  ################################## ###
  ### c Pulldown tag density Observed
  tags <- as.double(target$pulldown.tag.old)
  tags_ctrl <- as.double(target$ctrl.tag.old)
  
  tags <- replace( tags, tags == 0, 0.1 )
  tags_ctrl <- replace( tags_ctrl, tags_ctrl == 0, 0.1 )
  
  
  # a <-  (tags / tags_ctrl)
  log_c_ratio_observed <- log( tags / tags_ctrl )

  
  
  ################################## ###
  ### Define predicted tags
  
  ## Get all coefficients from the GLM (trained set)
  fit_coefficients <- as.data.frame( coefficients( summary( fit_object ) ) )
  
  ## Get the pulldown coefficients
  fit_pdwn_coefficients <- fit_coefficients[ grepl( pattern = ":t$", 
                                                    x = rownames( fit_coefficients ) ), ]

  ## Get the name of the variables included in th model    
  X_var_names <- gsub( ":t", "", rownames( fit_pdwn_coefficients ) )
  ## Make sure those variables are in the correct order for the matrix*vector
  X_test <- as.matrix( X[ , X_var_names ] )
  
  ## Mult. the PULLDOWN coeffs vector and the variable matrix
  c_pdwn_predicted <- as.vector( X_test %*% fit_pdwn_coefficients$Estimate ) # 1st must be 1.2476048 in test
  
  log_c_ratio_predicted <- c_pdwn_predicted
    
    
  # cor( log_c_ratio_predicted, log_c_ratio_observed )
    
  
  return( list( log_c_ratio_predicted = log_c_ratio_predicted, 
                log_c_ratio_observed = log_c_ratio_observed ) )
}



#################################################################################


############################################################################### #
################################### LOAD DATA ################################# #
#################################################################################

load_data <- function(input_root){

  ###################################################################### #
  ## Read PWM.
  fileName <- list.files(path=input_root, pattern="*_pfm.txt", full.names = T)
  rcade_pfm <- fread(fileName, sep="\t", data.table = F, header = T, skip="Pos")
  
  pfm_length <- max(rcade_pfm$Pos)
  rcade_pwm <- log( unlist( rcade_pfm[,2:5] ) + 0.001 )

  rm(rcade_pfm, fileName)
  
  ###################################################################### #
  ## Read Affimx and create big dataframe
  fileName <- list.files(path=input_root, pattern="*.affimx.affinity.txt",full.names = T)
  affimx <- read.csv(fileName,sep="\t", header = T, col.names = c("Name","Affinity"),
                     stringsAsFactors = FALSE)
  affimx$logAffinity <- log10(affimx$Affinity) # create an additional column with log of affinity (i.e. energy)
  rm(fileName)
  
  ###################################################################### #
  fileName <- list.files(path=input_root, pattern="*_aligned_sequences_tabulated_mx.txt",full.names = T)
  seq <- fread(fileName,sep="\t", data.table = F, header = F)

  # take the part of the sequence that's in range
  center <- as.integer( ( ncol(seq) - 5 ) / 2 + 5 ) # take the column that represents the begining of the motif hit, and use it to define the sequence range
  range <- ( center - flanking ):( center + pfm_length + flanking ) # ah_flank

        
  Cs <- (seq[,-1:-5]=="C") * 1 # find all the Cs
  Gs <- (seq[,-1:-5]=="G") * 1 # find all the Gs
  CpGs <- cbind( Cs[,-ncol(Cs)] * Gs[,-1], rep(0,nrow(Cs)) ) 
  ## find all Cs that are followed by a G (the C in a CpG)
  # also find any entry that follows a mark from the previous step (i.e. also mark the Gs in a CpG)
  CpGs[,-1] <- CpGs[,-1] + CpGs[,-ncol(CpGs)] 

  # convert to dataframe and add the sequence names
  CpGs <- cbind( seq[,1:5], as.data.frame( CpGs ) )
  
  seq <- seq[ , c(1,range) ]
  CpGs <- CpGs[ , c(1,range) ]
  
  # remove sequences with Ns. These sequences will later be removed from the CpG matrix based on sequence name
  seq <- seq[ apply( seq[,-1], 1, function(x) sum(x=="N") ) == 0, ]
  
  # Set the column names
  colnames(seq) <- c("Name",paste0("pos(",(-flanking-1):(pfm_length+flanking-1),")" ) ) # ah_flank
  colnames(CpGs) <- c("Name",paste0("pos(",(-flanking-1):(pfm_length+flanking-1),")" ) ) # ah_flank

  
  ###################################################################### #
  fileName <- list.files(path=input_root, pattern=paste0("*_sequences.fasta.methylreads"), full.names = T)
  methyl <- fread(fileName,sep="\t", data.table = F, header = F)

  # take the part of the sequence that's in range
  center <- as.integer( ( ncol(methyl) - 5 ) / 2 + 5 ) # ah_flank
  
  range <- ( center - flanking ):( center + pfm_length + flanking ) # ah_flank
  methyl <- methyl[ , c(1,range) ]
  # remove peaks with all NAs in the range
  methyl <- methyl[ apply( methyl[,-1], 1, function(x) sum(!is.na(x)) ) > 0, ]
  # Set the column names
  colnames(methyl) <- c("Name",paste0("pos","(",(-flanking-1):(pfm_length+flanking-1),")" ) ) # ah_flank

  
  ###################################################################### #
  fileName <- list.files(path=input_root, pattern=paste0("*_aligned_sequences.fasta.nonmethylreads"),full.names = T)
  nonmethyl <- fread(fileName,sep="\t", data.table = F, header = F)

  # take the part of the sequence that's in range
  center <- as.integer( ( ncol(nonmethyl)-5 ) / 2 + 5 ) # ah_flank
  range <- ( center - flanking ):( center + pfm_length + flanking )  # ah_flank
  
  nonmethyl <- nonmethyl[ , c(1,range) ]
  # remove peaks with all NAs in the range  
  nonmethyl <- nonmethyl[ apply( nonmethyl[,-1], 1, function(x) sum(!is.na(x)) ) > 0, ]
  # Set the column names
  colnames(nonmethyl) <- c("Name",paste0("pos","(",(-flanking-1):(pfm_length+flanking-1),")" ) ) # ah_flank

  
  ###################################################################### #
  fileName <- list.files(path=input_root, pattern=paste0("*_aligned_sequences.fasta.accessibility"), 
                         full.names = T)
  acc <- fread(fileName,sep="\t", data.table = F, header = F)                          
  # remove peaks with all NAs in the range
  acc <- acc[complete.cases(acc), ]                         
  acc_flanking <- 1200

  # take the part of the sequence that's in range
  center <- as.integer( ( ncol(acc)-5 ) / 2 + 5 ) # ah_flank
  range <- ( center - acc_flanking ):( center + pfm_length + acc_flanking ) # ah_flank
  acc <- acc[ , c(1,range) ] # 67139
  
  ## Replace 0 for pseudocount. As we will use log of DNA_acessibility for model and visualization.
  acc[,-1][acc[,-1] == 0] <- NA # All the columns except for the names.
  minimum <- min( acc[,-1], na.rm = TRUE )
  pseudocount <- median( sort(as.numeric(as.vector(as.matrix( acc[,-1]))) ), na.rm = TRUE ) * 0.01
  
  cat( paste0( "\nMinimum value for DNA accessibility (before average or log): ", 
               minimum, "\n", "Pseudocount: ", pseudocount, "\n" ) )
  
  acc[,-1][ is.na( acc[,-1] ) ] <- pseudocount

  # Set the column names
  colnames(acc) <- c("Name",paste0("pos","(",(-acc_flanking -1 ):(pfm_length + acc_flanking -1 ),")" ) )

  
  ###################################################################### #
  fileName <- list.files(path=input_root, pattern="*_summits_vRepeats_scores.txt",full.names = T)
  target <- fread(fileName,sep="\t", data.table = F, header = T)
  ## Makes sure that there are no entries with cero tags.
  # target <- target[target$tag != 0, ] # How many are we dropping here?
  # target <- target[target$MACS_score >= 1 , ]
  
  row.names(seq) <- seq$Name
  row.names(affimx) <- affimx$Name
  row.names(methyl) <- methyl$Name
  row.names(nonmethyl) <- nonmethyl$Name
  row.names(target) <- target$Name
  row.names(acc) <- acc$Name
  
  include <- Reduce( intersect, list( seq$Name, methyl$Name, nonmethyl$Name, 
                                      target$Name, acc$Name, affimx$Name ) )
  
  # include <- include[1:10000] # test with smaller number of peaks.
  acc <- acc[ acc$Name %in% include , ]
  acc <- acc[ order( acc$Name ) , ]

  affimx <- affimx[ affimx$Name %in% include , ]
  affimx <- affimx[ order( affimx$Name ) , ]

  seq <- seq[ seq$Name %in% include , ]
  seq <- seq[ order( seq$Name ) , ]

  CpGs <- CpGs[ CpGs$Name %in% include , ]
  CpGs <- CpGs[ order( CpGs$Name ) , ]

  methyl <- methyl[ methyl$Name %in% include , ]
  methyl <- methyl[ order( methyl$Name ) , ]
  
  nonmethyl <- nonmethyl[ nonmethyl$Name %in% include , ]
  nonmethyl <- nonmethyl[ order( nonmethyl$Name ) , ]
  
  target <- target[ target$Name %in% include , ]
  target <- target[ order( target$Name ) , ]
  
  
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
  
  return( list(acc = acc, target = target, affimx = affimx, 
               seq = seq, CpGs = CpGs, pfm_length = pfm_length,
               methyl = methyl, nonmethyl = nonmethyl, rcade_pwm = rcade_pwm) )
}
#################################################################################


############################################################################### #
############################### Descriptive plots ############################# #
#################################################################################

plot_heatmap_methyl_vs_target <- function( methyl_chip_data = methyl_chip_data, 
                                           ht_methyl_plot = ht_methyl_plot, 
                                           show_max_depth = T ) {
  # show_max_depth <- T
  pos_col <- paste0( "pos(", 0:( methyl_chip_data$pfm_length ), ")" )
  pos_col_cov <- paste0( pos_col, "_cov" ); pos_col_meth_pct <- paste0( pos_col, "_meth_pct" )
  tags_list <- c( "tags", "pulldown.tag.old", "ctrl.tag.old" )
  
  ## Define dataframes
  motif_methyl <- methyl_chip_data$methyl[ , c( "Name", pos_col ) ]
  motif_nonmethyl <- methyl_chip_data$nonmethyl[ , c( "Name", pos_col ) ]
  
  meth_cov <- motif_methyl[, pos_col] + motif_nonmethyl[, pos_col] ### Total coverage.
  meth_percent <- motif_methyl[, pos_col] * 100 / ( meth_cov ) ### Methylation percent.
  
  
  meth_cov$Name <- rownames( meth_cov ) ; meth_percent$Name <- rownames( meth_percent )
  
  meths_vs_tags <- merge( x = meth_cov, y = meth_percent, by = "Name",
                          suffixes = c( "_cov", "_meth_pct" ) )
  
  meths_vs_tags <- merge( x = meths_vs_tags, y = methyl_chip_data$target, by = "Name" )
  
  meths_vs_tags$pulldown.tag.old <- as.double( as.character(meths_vs_tags$pulldown.tag.old) )
  meths_vs_tags$ctrl.tag.old <- as.double( as.character(meths_vs_tags$ctrl.tag.old) )
  meths_vs_tags$tags <- as.double( as.character(meths_vs_tags$tags) )
  
  
  meths_vs_tags$tags_ht <- meths_vs_tags$pulldown.tag.old

  meths_vs_tags$tags_ht[ meths_vs_tags$pulldown.tag.old == 0 ] <- 0.1
  meths_vs_tags$tags_ht <- log2( meths_vs_tags$tags_ht )
  
  meths_vs_tags <- meths_vs_tags[order(meths_vs_tags$pulldown.tag.old, decreasing = T), ]
  
  
  ## Define scale for tag list.
  min_val <- min(meths_vs_tags[, "tags" ]); max_val <- max(meths_vs_tags[, tags_list] )
  med_val <- min_val + ( ( max_val - min_val ) / 2 )
  col_fun_tag_list <- circlize::colorRamp2( c( min_val, med_val, max_val ), c("white", "violet", "blue") )
  
  
  ### tags ht
  min_val <- min(meths_vs_tags[, "tags_ht" ]); max_val <- max(meths_vs_tags[, "tags_ht" ] )
  med_val <- min_val + ( ( max_val - min_val ) / 2 )
  col_fun_tags_ht <- circlize::colorRamp2( c( min_val, med_val, max_val ), c("white", "yellow", "green") )
  
  ## Define scale for meth. percentage.
  col_fun_meth_pct <- circlize::colorRamp2( c( 0, 50, 100 ), c("white", "orange", "red") )
  
  
  ## Define scale for meth. coverage.
  min_val <- min(meths_vs_tags[, pos_col_cov ], na.rm = T ) # Probably 0.
  max_val <- max(meths_vs_tags[, pos_col_cov ], na.rm = T  )
  med_val <- min_val + ( ( max_val - min_val ) / 2 )
  
  if( !show_max_depth ){ min_val <- 0; med_val <- 25; max_val <- 50 }
  
  col_fun_cov <- circlize::colorRamp2( c( min_val, med_val, max_val ), c("white", "orange", "red") )
  
  ################################################################################## #
  
  this_width <- 3 + 0.6 * ( methyl_chip_data$pfm_length )
  this_height <- 5 + ( dim( methyl_chip_data$target )[1] / 3000 )

  this_width <- 15
  this_height <- 8
    
  pdf( file = ht_methyl_plot, width = this_width, height = this_height )
  
  ht_meth_pct <- ComplexHeatmap::Heatmap( matrix = as.matrix( meths_vs_tags[, pos_col_meth_pct ] ),
                                          cluster_columns = F, cluster_rows = F, show_row_names = F,
                                          na_col = "black",col = col_fun_meth_pct,
                                          name = "Meth. percentage", column_title = "Meth. percentage",
                                          raster_device = "png",
                                          use_raster = T, raster_quality = 2 )

  ht_cov <- ComplexHeatmap::Heatmap( matrix = as.matrix( meths_vs_tags[, pos_col_cov ] ),
                                     cluster_columns = F, cluster_rows = F, show_row_names = F, 
                                     na_col = "black",col = col_fun_cov,
                                     name = "WGBS coverage", column_title = "WGBS coverage",
                                     raster_device = "png",
                                     use_raster = T, raster_quality = 2 )

  
  ht_tag_list <- ComplexHeatmap::Heatmap( matrix = as.matrix( meths_vs_tags[, tags_list ] ),
                                          cluster_columns = F, cluster_rows = F,
                                          show_row_names = F, col = col_fun_tag_list,
                                          name = "Tags ", column_title = "",
                                          raster_device = "png",
                                          use_raster = T, raster_quality = 2 )
  
  ht_tags <- ComplexHeatmap::Heatmap( matrix = as.matrix( meths_vs_tags[, "tags_ht" ] ),
                                      cluster_columns = F, cluster_rows = F,
                                      show_row_names = F, col = col_fun_tags_ht,
                                      name = "log2( Pulldown tags )", column_title = "",
                                      raster_device = "png",
                                      use_raster = T, raster_quality = 2 )
  
  ht_list <- ht_cov + ht_meth_pct + ht_tag_list + ht_tags
  
  draw( ht_list, column_title = "Methylation vs ChIP-seq tags" )
  dev.off()
  
  rm( ht_list, ht_meth_pct, ht_cov, ht_tag_list, ht_tags, col_fun_meth_pct, 
      col_fun_tags_ht, meths_vs_tags, meth_cov, meth_percent, pos_col_cov, pos_col,
      tags_list, motif_nonmethyl, motif_methyl, col_fun_cov, min_val, med_val, max_val)
  
  
}

plot_heatmap_sequence_vs_target <- function( methyl_chip_data = methyl_chip_data, 
                                             ht_seq_plot = ht_seq_plot ){
  
  tags_list <- c( "tags", "pulldown.tag.old", "ctrl.tag.old" )
  seq_col <- paste0( "pos(", -5:(methyl_chip_data$pfm_length+5), ")" )
  
  motif_sequence <- methyl_chip_data$seq[ , c( "Name", seq_col ) ]
  seq_vs_tags <- merge( x = motif_sequence, y = methyl_chip_data$target, by = "Name")
  seq_vs_tags <- seq_vs_tags[order(seq_vs_tags$pulldown.tag.old, decreasing = T), ]
  
  
  seq_vs_tags$tags_ht <- seq_vs_tags$pulldown.tag.old
  seq_vs_tags$tags_ht[ seq_vs_tags$tags_ht <= 0 ] <- 0.1
  seq_vs_tags$tags_ht <- log2( seq_vs_tags$tags_ht )  
  
  
  min_val <- min(seq_vs_tags[, tags_list ] ); max_val <- max(seq_vs_tags[, tags_list ] )
  med_val <- min_val + ( ( max_val - min_val ) / 2 )
  col_fun_tags <- circlize::colorRamp2( c( min_val, med_val, max_val ), c("white", "violet", "blue") )
  
  
  min_val <- min(seq_vs_tags[, "tags_ht" ] ); max_val <- max(seq_vs_tags[, "tags_ht" ] )
  med_val <- min_val + ( ( max_val - min_val ) / 2 )
  col_fun_tags2 <- circlize::colorRamp2( c( min_val, med_val, max_val ), c("white", "yellow", "green") )
  
  
  min_val <- min(seq_vs_tags[, "pulldown.tag.old" ]); max_val <- max(seq_vs_tags[, "pulldown.tag.old"])
  med_val <- min_val + ( ( max_val - min_val ) / 2 )
  col_fun <- circlize::colorRamp2( c( min_val, max_val ), c("white", "blue") )
  
  colors_bases <- structure( c("green", "orange", "blue", "red"), names = c( "A", "G", "C", "T" ) )
  
  this_width <- 4 + 0.6 * ( methyl_chip_data$pfm_length )
  this_height <- 5 + ( dim( methyl_chip_data$target )[1] / 3000 )
  this_width <- 15
  this_height <- 8

  pdf( file = ht_seq_plot, width = this_width, height = this_height )

  ht_main <- ComplexHeatmap::Heatmap( matrix = as.matrix( seq_vs_tags[, seq_col ] ),
                                      cluster_columns = F, cluster_rows = F,
                                      show_row_names = F, col = colors_bases,
                                      name = "Sequence",
                                      raster_device = "png",
                                      use_raster = T, raster_quality = 2 )
  
  ht_tags <- ComplexHeatmap::Heatmap( matrix = as.matrix( seq_vs_tags[, tags_list ] ),
                                      cluster_columns = F, cluster_rows = F,
                                      show_row_names = F, col = col_fun_tags,
                                      name = "Tags ",
                                      raster_device = "png",
                                      use_raster = T, raster_quality = 2 )
  
  ht_tags2 <- ComplexHeatmap::Heatmap( matrix = as.matrix( seq_vs_tags[, "tags_ht" ] ),
                                       cluster_columns = F, cluster_rows = F,
                                       show_row_names = F, col = col_fun_tags2,
                                       name = "log2( Pulldown tags) ",
                                       raster_device = "png",
                                       use_raster = T, raster_quality = 2 )

    
  ht_list <- ht_main + ht_tags + ht_tags2
  draw( ht_list, column_title = "Motif sequence vs ChIP-seq tags" )
  dev.off()
  
  rm(ht_main, ht_list, seq_vs_tags, seq_col, motif_sequence, 
     min_val, max_val, med_val, colors_bases )
}

plot_heatmap_dna_acc_vs_target <- function( methyl_chip_data = methyl_chip_data, 
                                            dna_acc_by_bins = dna_acc_by_bins,
                                            ht_acc_plot = ht_acc_plot ){
  
  keep_dna_bins <- c( "bin_up_5", "bin_up_4", "bin_up_3", "bin_up_2", "bin_up_1", "motif_bin", 
                      "bin_down_1", "bin_down_2", "bin_down_3", "bin_down_4", "bin_down_5") 
  
  tags_list <- c( "tags", "pulldown.tag.old", "ctrl.tag.old" ) 
  
  dna_acc_vs_tags <- merge( x = dna_acc_by_bins, y = methyl_chip_data$target, by = "Name")
  dna_acc_vs_tags <- subset(dna_acc_vs_tags, select = - Name )
  
  dna_acc_vs_tags$tags_ht <- as.numeric( as.character( dna_acc_vs_tags$pulldown.tag.old ) )
  dna_acc_vs_tags$tags_ht[ dna_acc_vs_tags$tags_ht == 0 ] <- 0.1
  dna_acc_vs_tags$tags_ht <- log2( dna_acc_vs_tags$tags_ht )  
  
  dna_acc_vs_tags <- dna_acc_vs_tags[order(dna_acc_vs_tags$pulldown.tag.old, decreasing = T), ]
  
  min_val <- min(dna_acc_vs_tags[, keep_dna_bins ]); max_val <- max(dna_acc_vs_tags[, keep_dna_bins])
  med_val <- min_val + ( ( max_val - min_val ) / 2 )
  col_fun <- circlize::colorRamp2( c( min_val, med_val, max_val ), c("white", "orange", "red") )
  
  
  min_val <- min(dna_acc_vs_tags[, tags_list ] ); max_val <- max(dna_acc_vs_tags[, tags_list ] )
  med_val <- min_val + ( ( max_val - min_val ) / 2 )
  col_fun_tags <- circlize::colorRamp2( c( min_val, med_val, max_val ), c("white", "violet", "blue") )
  
  
  min_val <- min(dna_acc_vs_tags[, "tags_ht" ] ); max_val <- max(dna_acc_vs_tags[, "tags_ht" ] )
  med_val <- min_val + ( ( max_val - min_val ) / 2 )
  col_fun_tags2 <- circlize::colorRamp2( c( min_val, med_val, max_val ), c("white", "yellow", "green") )
  
  
  this_width <- 10
  this_height <- 5 + ( dim( methyl_chip_data$target )[1] / 3000 )
  this_width <- 15
  this_height <- 8

  pdf( file = ht_acc_plot, width = this_width, height = this_height )
  
  
  ht_main <- ComplexHeatmap::Heatmap( matrix = as.matrix( dna_acc_vs_tags[, keep_dna_bins ] ),
                                      cluster_columns = F, cluster_rows = F,
                                      show_row_names = F, col = col_fun, 
                                      name = "DNA acc. by bins",
                                      raster_device = "png",
                                      use_raster = T, raster_quality = 2 )
  
  ht_tags <- ComplexHeatmap::Heatmap( matrix = as.matrix( dna_acc_vs_tags[, tags_list ] ),
                                      cluster_columns = F, cluster_rows = F,
                                      show_row_names = F, col = col_fun_tags,
                                      name = "Tags ",
                                      raster_device = "png",
                                      use_raster = T, raster_quality = 2 )
  
  ht_tags2 <- ComplexHeatmap::Heatmap( matrix = as.matrix( dna_acc_vs_tags[, "tags_ht" ] ),
                                       cluster_columns = F, cluster_rows = F,
                                       show_row_names = F, col = col_fun_tags2,
                                       name = "log2( Pulldown tags) ",
                                       raster_device = "png",
                                       use_raster = T, raster_quality = 2 )
  
  ht_list <- ht_main + ht_tags + ht_tags2
  draw( ht_list, column_title = "DNA acc. around motif vs ChIP-seq tags" )
  dev.off()
}

plot_dna_acc_around_motif <- function( dna_acc_by_bins = dna_acc_by_bins, 
                                       methyl_chip_data = methyl_chip_data,
                                       acc_motif_plot_name = acc_motif_plot_name ){
  
  ### Plot DNA accessibility by position
  cat(paste0("Ploting: DNA acc. by position ...  \n"))
  means <- colMeans(methyl_chip_data$acc[,-c(1,2)])
  center <- floor( ncol(methyl_chip_data$acc[,-c(1,2)]) / 2 )
  even <- ncol(methyl_chip_data$acc[,-c(1,2)]) / 2
  
  if (even%%1==0) { pos <- -center:(center-1) } else{ pos <- -center:(center) }
  
  acc_data <- data.frame( pos = pos, means = means)
  
  plot_acc <- ggplot(data = acc_data, aes(x = pos, y = means) ) +
    geom_line(size = 0.5, colour = "black") +
    xlab("TF motif") + ylab("Normalized DNA accessibility") + 
    ggtitle("DNA accessibility around the TF motif") 
  
  # + theme(plot.title = element_text(hjust = 0.5))
  
  acc_mod <- reshape2::melt(data = dna_acc_by_bins[,-1], 
                            measure.vars = colnames(dna_acc_by_bins[,-1]) )
  
  plot_violin_acc <- ggplot(acc_mod, aes(x=variable, y=value)) + 
    geom_violin() + xlab("") + ylab("log(Norm. DNA acc.)")
  
  my_fig1 <- ggarrange(plot_acc, plot_violin_acc, nrow = 2, ncol = 1)
  
  ggsave( filename = acc_motif_plot_name, plot = my_fig1, height = 7, width = 10)
}

get_height <- function( df = df, this_width = this_width ){
  
  num_rows <- nrow( df )
  if( num_rows <= 1000 ){ this_height <- this_width*0.65 }
  if( num_rows > 1000 & num_rows <= 50000 ){ this_height <- this_width }
  if( num_rows > 50000 & num_rows <= 200000 ){ this_height <- this_width*1.5 }
  if( num_rows > 200000 ){ this_height <- this_width*3 }
  
  return( this_height ) }


create_hm_all_data <- function( methyl_chip_data = methyl_chip_data, 
                                meth_data = fit_CpG_only$x.Met,
                                acc_data = dna_acc_by_bins,
                                pulldown_motif = model.mCpG_only$motif_pulldown,
                                out_ht_path = out_ht_path,
                                experiment =  opt$experiment,
                                top_acc = top_acc ){

  # methyl_chip_data <- methyl_chip_data
  # meth_data <- fit_CpG_only$x.Met
  # acc_data <- dna_acc_by_bins
  # pulldown_motif <- model.mCpG_only$motif_pulldown
  # out_ht_path <- out_ht_path
  # experiment <- opt$experiment
  # top_acc <- 1
  
  # meth_data = fit_CpG_only$x.Met,
  ### Prepare all data for heatmap.
  colnames( meth_data ) <- paste0( "av_meth_", colnames( meth_data ) )
  meth_data$Name <- rownames( meth_data )

  pos_col <- c( "Name", paste0( "pos(", 0 : ( methyl_chip_data$pfm_length -1 ), ")" ) )
  
  seqs <- methyl_chip_data$seq[ , pos_col ]
  colnames( seqs ) <- c( "Name", paste0( "seq_", colnames( seqs[ , -1 ] ) ) )
  
  acc_col <- colnames( acc_data[, -1 ] )
  
  affimx_col <- colnames(  methyl_chip_data$affimx[, -1 ] )
  
  seq_col <- paste0( "seq_pos(",  0 : ( methyl_chip_data$pfm_length -1 ), ")" )
  
  meth_col <- paste0( "av_meth_pos(",  0 : ( methyl_chip_data$pfm_length -1 ), ")" )
  
  target_col <- c( "ctrl.tag.old", "pulldown.tag.old" )
  
  dfs <- list( methyl_chip_data$affimx, seqs, meth_data, acc_data, methyl_chip_data$target )
  
  
  all_dat <- Reduce( function( x, y ) merge( x, y, all = TRUE, by = "Name" ), dfs )
  rownames( all_dat ) <- all_dat$Name
  
  
  all_dat$log_pdwn <- log2( as.double( as.character( all_dat$pulldown.tag.old  ) ) )
  all_dat$ctrl.tag.old <- as.double( as.character( all_dat$ctrl.tag.old  ) )
  all_dat$pulldown.tag.old <- as.double( as.character( all_dat$pulldown.tag.old  ) )
  
  
  
  
  ############################################## #
  ### Mod matrix
  pulldown_motif <- model.mCpG_only$motif_pulldown
  pulldown_motif <- t( pulldown_motif[1:4,] )
  
  pulldown_motif[,1] <- as.double( as.character( pulldown_motif[,1] ) )
  pulldown_motif[,2] <- as.double( as.character( pulldown_motif[,2] ) )
  pulldown_motif[,3] <- as.double( as.character( pulldown_motif[,3] ) )
  pulldown_motif[,4] <- as.double( as.character( pulldown_motif[,4] ) )
  
  pulldown_motif <- as.data.frame( pulldown_motif )
  pulldown_motif[ pulldown_motif < 0 ] <- 0
  
  pulldown_motif$colMax <- apply( pulldown_motif, 1, function( x ) max( x ) )
  
  pulldown_motif[ pulldown_motif$colMax == pulldown_motif$A, "seq" ] <- "A"
  pulldown_motif[ pulldown_motif$colMax == pulldown_motif$G, "seq" ] <- "G"
  pulldown_motif[ pulldown_motif$colMax == pulldown_motif$C, "seq" ] <- "C"
  pulldown_motif[ pulldown_motif$colMax == pulldown_motif$T, "seq" ] <- "T"
  pulldown_motif[ pulldown_motif$colMax < ( max( pulldown_motif$colMax )*0.3 ) , "seq" ] <- "any"


  for( i in 1:length( pulldown_motif$seq ) ){
    if( pulldown_motif$seq[i] != "any" ){
      pos_name <- paste0( "seq_pos(", ( i - 1 ) , ")" )
      all_dat <- all_dat[ all_dat[ , pos_name ]  == pulldown_motif$seq[i] , ] } }
  
  
  if (  nrow( all_dat ) == 0 ) {
    
    cat( "No peak with all high coefficient bases\n" )
    
  } else{
    
    ############################################## #
    ### high accessibility peaks
    # bin_col <- c( "bin_up_2", "bin_up_1", "motif_bin", "bin_down_1", "bin_down_2" )
    bin_col <- c( "bin_up_1", "motif_bin", "bin_down_1" )
    all_dat <- all_dat[ order( rowMeans( all_dat[, bin_col] ), decreasing = TRUE ), ]
    
    ############################################## #
    ### Top % peaks (by DNA acc. )
    all_dat <- all_dat[ 1:ceiling( nrow( all_dat )*top_acc ), ]
    
    ############################################## #
    ### Sort by log( tags )
    all_dat <- all_dat[ order( all_dat$log_pdwn, decreasing = TRUE ), ]
    
    
    ############################################## #
    ### Colors 
    colors_affimx <- colorRamp2( c( min( all_dat$logAffinity ), 0 ), c( "blue", "white" ) )
    # colors_affimx <- colorRamp2( c( -20, 0 ), c( "blue", "white" ) )
    
    colors_bases <- structure( c("green", "orange", "blue", "red"), names = c( "A", "G", "C", "T" ) )
    colors_meth <- colorRamp2( c( 0, 0.5, 1 ), c( "white", "orange", "red" ) )
    
    colors_acc <- colorRamp2( c( min( all_dat[, acc_col ] ) - 0.001, max( all_dat[, acc_col ] ) + 0.001 ), 
                              c( "white", "red" ) )
    
    colors_count_tags <- colorRamp2( c( 0,  max( all_dat[, target_col] ) ), 
                                     c( "white", "red" ) )
    
    colors_log_tags <- colorRamp2( c( min( all_dat$log_pdwn ) - 0.001, max( all_dat$log_pdwn ) + 0.001 ),
                                   c( "white", "red" ) )
    
    # min( all_dat$log_pdwn )
    # max( all_dat$log_pdwn )
    
    ############################################## #
    ### Plot
    pdf( out_ht_path, width = 15, height =  get_height( df = all_dat, this_width = 15 ) )
    ## Affimx
    ht1 <- ComplexHeatmap::Heatmap( matrix = as.matrix( all_dat[ , "logAffinity" ] ) ,
                                    cluster_columns = F, cluster_rows = F,
                                    col = colors_affimx,
                                    column_title_rot = 90,
                                    column_title = "logAffimx",
                                    column_title_side = "bottom",
                                    show_row_names = F, show_column_names = T,
                                    name = "logAffimx", raster_device = "png",
                                    use_raster = T, raster_quality = 2 )
    
    ## sequence
    ht2 <- ComplexHeatmap::Heatmap( matrix = as.matrix( all_dat[ , seq_col ] ),
                                    column_title_rot = 0,
                                    column_title = "Motif sequence by position",
                                    column_title_side = "bottom",
                                    cluster_columns = F, cluster_rows = F,
                                    col = colors_bases,
                                    show_row_names = F, show_column_names = T,
                                    name = "Sequence", raster_device = "png",
                                    use_raster = T, raster_quality = 2 )
    ## DNA acc.
    ht3 <- ComplexHeatmap::Heatmap( matrix = as.matrix( all_dat[, acc_col] ),
                                    column_title_rot = 0,
                                    column_title = "Acc",
                                    column_title_side = "bottom",
                                    cluster_columns = F, cluster_rows = F,
                                    col = colors_acc,
                                    show_row_names = F, show_column_names = T,
                                    name = "DNA acc.", raster_device = "png",
                                    use_raster = T, raster_quality = 2 )
    ## Methylation
    ht4 <- ComplexHeatmap::Heatmap( matrix = as.matrix( all_dat[, meth_col] ),
                                    column_title_rot = 0,
                                    column_title = "Methylation",
                                    column_title_side = "bottom",
                                    cluster_columns = F, cluster_rows = F,
                                    col = colors_meth,
                                    show_row_names = F, show_column_names = T,
                                    name = "Methylation", raster_device = "png",
                                    use_raster = T, raster_quality = 2 )
    ## Count tags
    ht5 <- ComplexHeatmap::Heatmap( matrix = as.matrix( all_dat[, target_col] ),
                                    column_title_rot = 90,
                                    column_title = "tag",
                                    column_title_side = "bottom",
                                    cluster_columns = F, cluster_rows = F,
                                    col = colors_count_tags,
                                    show_row_names = F, show_column_names = T,
                                    name = "count tags", raster_device = "png",
                                    use_raster = T, raster_quality = 2 )
    ## log2( Tags )
    ht6 <- ComplexHeatmap::Heatmap( matrix = as.matrix( all_dat$log_pdwn ),
                                    column_title_rot = 90,
                                    column_title = "log2(tag)",
                                    column_title_side = "bottom",
                                    cluster_columns = F, cluster_rows = F,
                                    col = colors_log_tags,
                                    show_row_names = F, show_column_names = T,
                                    name = "log2(tag)", raster_device = "png",
                                    use_raster = T, raster_quality = 2 )
    
    ht_list <- ht1 + ht2 + ht3 + ht4 + ht5 + ht6
    
    draw( ht_list, column_title = paste0( experiment, "\nn = ", nrow(all_dat) ) )
    
    dev.off( )
    # acc_col  # affimx_col  # seq_col  # meth_col    # target_col
  
  }

}





heatmap_meth_affinity <- function( all_dat = all_dat,
                                   fit_object = fit_CpG_only,
                                   outdir = outdir, 
                                   experiment_name = opt$experiment ){
  # all_dat <- all_dat
  ## Read coefficients & correct for FDR
  fit_object <- fit_CpG_only
  coefs_summary <- as.data.frame( fit_object$coefs_summary )
  
  # ## Example comment in deployment
  # coefs_summary <- read.csv( file = in_path, sep = "\t" )
  # rownames( coefs_summary ) <- coefs_summary$X
  # coefs_summary <- coefs_summary[, 2:5]
  # coefs_summary <- as.data.frame( coefs_summary )
  # ## 
  
  coefs_summary <- cbind( coefs_summary, p.adjust( coefs_summary[,4], 
                                    method = "fdr" ) ) 
  colnames( coefs_summary )[5] <- "FDR"
  coefs_summary[ coefs_summary$FDR >= 0.00001, ]$FDR <- NA
  
  ## Decide methylation position to focus in.
  met_coefs <- coefs_summary[ grepl( pattern = "x\\.Met\\.pos\\..*\\.:t$", 
                                     x = rownames(coefs_summary)), ]

  met_coefs <- met_coefs[ ! is.na( met_coefs$FDR ), ]
  
  
  if ( nrow(met_coefs) >= 1 ){
    ## Sort by FDR
    met_coefs <- met_coefs[ order( met_coefs$FDR, decreasing = FALSE ), ]    
    
    target_positions <- gsub( "x\\.Met\\.pos\\.", "", rownames(met_coefs)[1] )
    target_positions <- as.integer( gsub( "\\.:t$", "", target_positions ) )
    
    ##########################################################  ##
    ## Delete rows with 1 in x.CG.pos.( position ).
    CG_pos <- paste0( "x.CG.pos.", target_positions, "." )

    CG_peaks <- all_dat[ all_dat[, CG_pos] == 1 , ]
    meth_coef <- paste0( "x.Met.pos.", target_positions, "." )
    
    
    ## Define data to use in glm
    t_indictator <- as.data.frame(c(rep(0, nrow( CG_peaks )), rep(1, nrow( CG_peaks ))))
    colnames( t_indictator ) <- "t"
    c_tags <- as.data.frame(c(CG_peaks$tags_ctrl, CG_peaks$tags))
    colnames( c_tags ) <- "c"
    X <- subset ( CG_peaks, select = -c( tags, tags_ctrl ) )
    XX <- rbind( X, X )
    
    glm.data <- cbind( XX, c_tags, t_indictator)
    rownames( glm.data ) <- c( paste0( "control.", rownames( X ) ),
                               paste0( "pulldown.", rownames( X ) ) ) 
    
    #### Define formula to use in glm
    ## Define control coeffs
    X.ctl <- colnames( X )
    ## Remove position of interest in X and define X:t (X pulldown)
    X.pdn <- setdiff( X.ctl, meth_coef )
    
    XX.c.formula <- paste0( " ", X.ctl, " +", collapse = "" ) 
    XX.p.formula <- paste0( " ", X.pdn, ":t +", collapse = "" ) 
    
    my.formula <- paste0( "c ~ ", XX.c.formula, XX.p.formula, " t", collapse = "" )
    # my.formula <- paste0( "c ~ acc.motif_bin:t"   , collapse = "" ) 
    
    ## Regression without pos in interest and only CpG sites.
    fit.CG <- glm.nb( my.formula, data = glm.data )
    CG_fit_path <- paste0( outdir, "/", experiment_name, "_CG_fit.rda" )
    save( fit.CG, file = CG_fit_path )
    
  }
  }






 

prepare_dat <- function( seq, CpGs, methyl, nonmethyl, dna_acc, 
                       target, pfm_length, flanking = 0, CpG_only = T) {


  ## Preparation
  # extract the sequence matrix and convert it to one-hot encoding
  x.A <- ( seq[,-1] == "A" ) * 1
  x.C <- ( seq[,-1] == "C" ) * 1
  x.G <- ( seq[,-1] == "G" ) * 1
  x.T <- ( seq[,-1] == "T" ) * 1
  rownames(CpGs) <- CpGs$Name
  x.CpG <- CpGs[,-1] * 1
  
  ncolx <- ncol(x.A)
  x.CG <- x.C[,-ncolx] * x.G[,-1]
  
  ## Calculate the fraction of methylated C's in the reverse strand (in front of every G)
  x.W <- x.G * ( methyl[,-1] + 1 ) / ( methyl[,-1] + nonmethyl[,-1] + 2 ) 
  ## Calculate the fraction of methylated C's in the forward strand
  x.M <- x.C * ( methyl[,-1] + 1 ) / ( methyl[,-1] + nonmethyl[,-1] + 2 ) 
  
  ## Methylation in reverse strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.W[ is.na(x.W) ] <- 0
  # only keep the W's that are within a CpG
  if( CpG_only ) { x.W <- x.W * x.CpG }
  ## Methylation in forward strand
  # replace NA values by zero (these correspond to A/T nucleotides)
  x.M[ is.na(x.M) ] <- 0
  # only keep the W's that are within a CpG
  if( CpG_only ) { x.M <- x.M * x.CpG }
  
  
  ## Methylation average M (forward) W (reverse).
  x.Met <- ( ( x.M[,-ncolx] + x.W[,-1] ) / 2 )
  
  
  motif <- paste0( "pos(", 0:(pfm_length - 1) , ")" )
  upsteam_flank <- paste0( "pos(", (-flanking):-1 , ")" )
  downsteam_flank <- paste0( "pos(", (pfm_length ):(pfm_length + flanking - 1), ")" )

  
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
  acc <- calculate_mean_by_bins(acc = dna_acc, pfm_length = pfm_length)

  ## Tags
  tags <- as.factor(target$pulldown.tag.old)
  tags_ctrl <- as.factor(target$ctrl.tag.old)
  
  dat <- data.frame( acc = acc[,-1], 
                     x.C = x.C, x.G = x.G, x.T = x.T,
                     x.CG = x.CG, x.Met = x.Met,
                     x.T_up = x.T_up, x.C_up = x.C_up, x.G_up = x.G_up,
                     x.M_up = x.M_up, x.W_up = x.W_up,
                     x.T_down = x.T_down, x.C_down = x.C_down,x.G_down = x.G_down,
                     x.M_down = x.M_down, x.W_down = x.W_down, 
                     tags = tags, tags_ctrl = tags_ctrl )
  
  return( dat )
  
}
  

