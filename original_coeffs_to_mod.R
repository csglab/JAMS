require(data.table)
require(optparse)
require(ggplot2)
require(patchwork)

options( error = traceback, nwarnings = 10000 )
source(  "/home/ahcorcha/repos/tools/JAMS/src/Methyl_ChIP_ftns.R" )
setwd( "/home/ahcorcha/repos/tools/JAMS" )

plot_dna_acc_coefficients <- function( coefs, plot_name ){
  
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
  y_min <- min( c( ( dna_coefs_ctrl$Estimate - dna_coefs_ctrl$Std..Error ), 
                   ( dna_coefs_pdwn$Estimate - dna_coefs_pdwn$Std..Error ) ) )
  
  y_max <- max( c( ( dna_coefs_ctrl$Estimate + dna_coefs_ctrl$Std..Error ), 
                   ( dna_coefs_pdwn$Estimate + dna_coefs_pdwn$Std..Error ) ) )
  
  mag <- abs( y_max - y_min )
  center <- y_min + ( mag / 2 )
  
  y_min <- center - ( mag / 2 )*1.2
  y_max <- center + ( mag / 2 )*1.2
  #################################################################### ## 
  
  dna_pdwn_plot <- ggplot( data = dna_coefs_pdwn, aes( x = Name, y = Estimate, color = FDR_color ) ) +
                           geom_point( size = 1 ) +
                           geom_errorbar( aes( x = Name, width=0.5,
                                               ymin = Estimate - Std..Error, 
                                               ymax = Estimate + Std..Error ) ) +
                           theme_minimal() + labs( x = "Pulldown", y = "Estimate" ) +
                           theme( axis.text.x = element_text( angle = 290, hjust = 0 ),
                                  legend.position = "none" ) + 
                           ylim( y_min, y_max ) + 
                           scale_color_manual( values = c( "grey", "black" ) )
  
  dna_ctrl_plot <- ggplot( data = dna_coefs_ctrl, aes( x = Name, y = Estimate, color = FDR_color ) ) +
                           geom_point( size = 1 ) +
                           geom_errorbar( aes( x = Name, width=0.5,
                                               ymin = Estimate - Std..Error, 
                                               ymax = Estimate + Std..Error ) ) +
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

option_list = list(
  make_option(c("-i", "--input"), type="character",
              default="./data/deprecated/test/CTCF_HEK293_GSM2026781_small_mCpG_only_model_coefficients_with_FDR.txt",
              help="original coefficients with FDR", metavar="character"),

  make_option(c("-o", "--output_dir"), type="character",
              default="./data/deprecated",
              help="", metavar="character") );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)

original_coeffs <- read.table( file = opt$input, header = TRUE, sep = "\t" )


original_coeffs$Coefficient_name <- rownames( original_coeffs )
original_coeffs$Predictor <- original_coeffs$Coefficient_name


original_coeffs[grepl(pattern="^acc.*$", x=original_coeffs$Coefficient_name), "Predictor" ] <- "DNA_accessibility"
original_coeffs[grepl(pattern="^x\\..*\\.pos\\..*$", x=original_coeffs$Coefficient_name), "Predictor" ] <- "Motif_Sequence"
original_coeffs[grepl(pattern="^x\\.CG\\.pos\\..*$", x=original_coeffs$Coefficient_name), "Predictor" ] <- "Motif_CpG_state"
original_coeffs[grepl(pattern="^x\\.Met\\.pos\\..*$", x=original_coeffs$Coefficient_name), "Predictor" ] <- "Motif_methylation"

original_coeffs[grepl(pattern="^x\\..*_up.*$", x=original_coeffs$Coefficient_name), "Predictor" ] <- "Upstream_base_composition"
original_coeffs[grepl(pattern="^x\\.M_up.*$", x=original_coeffs$Coefficient_name), "Predictor" ] <- "Upstream_methylation_composition"
original_coeffs[grepl(pattern="^x\\.W_up.*$", x=original_coeffs$Coefficient_name), "Predictor" ] <- "Upstream_methylation_composition"

original_coeffs[grepl(pattern="^x\\..*_down.*$", x=original_coeffs$Coefficient_name), "Predictor" ] <- "Downstream_base_composition"
original_coeffs[grepl(pattern="^x\\.M_down.*$", x=original_coeffs$Coefficient_name), "Predictor" ] <- "Downstream_methylation_composition"
original_coeffs[grepl(pattern="^x\\.W_down.*$", x=original_coeffs$Coefficient_name), "Predictor" ] <- "Downstream_methylation_composition"


original_coeffs <- original_coeffs[ !original_coeffs$Coefficient_name %in% c("(Intercept)", "t"), ]

original_coeffs <- original_coeffs[,c( "Predictor", "Coefficient_name", 
                                       "Estimate", "Std..Error", "z.value", 
                                       "Pr...z..", "FDR" )]

background_coeffs <- original_coeffs[ !grepl( pattern = "^.*:t$", x = original_coeffs$Coefficient_name ), ]

TF_binding_coeffs <- original_coeffs[ grepl( pattern = "^.*:t$", x = original_coeffs$Coefficient_name ), ]

TF_binding_coeffs$Coefficient_name <- gsub( pattern = ":t$", replacement = "", x = TF_binding_coeffs$Coefficient_name )


write.table( x = TF_binding_coeffs, 
           file = paste0( opt$output_dir, "/", gsub( ".txt", "", basename(opt$input) ), "_TF_binding.txt"), 
           sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE )

write.table( x = background_coeffs, 
           file = paste0( opt$output_dir, "/", gsub( ".txt", "", basename(opt$input) ), "_background.txt"), 
           sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE )


original_coeffs_ <- as.data.frame( read.table( file = opt$input, header = TRUE, sep = "\t" ) )

plot_dna_acc_coefficients( original_coeffs_, 
                           paste0( opt$output_dir, "/", gsub( ".txt", "", basename(opt$input) ), "_dna_coefficients.pdf" ) )


