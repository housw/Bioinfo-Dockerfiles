# Visualization of REVIGO csv files
# 
#call: R --slave -f REVIGO_plotter.R --args workdir=workdir REVIGO_BP=REVIGO_BP REVIGO_CC=REVIGO_CC REVIGO_MF=REVIGO_MF out_pdf=out_pdf

args <- commandArgs(trailingOnly = TRUE)
workdir <- "~/MEGA/GOstats/"
REVIGO_BP <- "HG001_00009__HG001_03232_CopraRNA_result_locustags_GO_enrichment_BP_revigo.csv"
REVIGO_CC <- "HG001_00009__HG001_03232_CopraRNA_result_locustags_GO_enrichment_CC_revigo.csv"
REVIGO_MF <- "HG001_00009__HG001_03232_CopraRNA_result_locustags_GO_enrichment_MF_revigo.csv"
out_pdf <- "HG001_00009__HG001_03232_CopraRNA_result_locustags_GO_enrichment.pdf"

for(i in 1:length(args)){
  temp<-strsplit(args[i],"=")
  temp<-temp[[1]]
  temp1<-temp[1]
  temp2<-temp[2]
  assign(as.character(temp1),temp2)
}

setwd(workdir)


# ----------------------------------------------------------------
# read in csv from REVIGO
# ----------------------------------------------------------------

# read in GO enrichment results
read_REVIGO_results <- function(REVIGO_BP, REVIGO_CC, REVIGO_MF){
  
  # case for "GO \t p-value" table as input for REVIGO, 11 output columns
  revigo.names_11 <- c("term_ID","description","frequency", 
                    "plot_X","plot_Y","plot_size",
                    "log10_p_value","uniqueness","dispensability", 
                    "representative", "eliminated");
  # case for "GO \t p-value \t count" table as input for REVIGO, 12 output columns
  revigo.names_12 <- c("term_ID","description","frequency", 
                       "plot_X","plot_Y","plot_size",
                       "log10_p_value","gene_counts", "uniqueness","dispensability", 
                       "representative", "eliminated");
  
  # read in REVIGO_BP
  df.REVIGO_BP <- read.csv(REVIGO_BP, header = T, stringsAsFactors = F)
  if (ncol(df.REVIGO_BP) == 11) {
    revigo.names = revigo.names_11
  } else {
    revigo.names = revigo.names_12
  }
  colnames(df.REVIGO_BP) <- revigo.names    
  cat_BP <- rep("BP", nrow(df.REVIGO_BP))
  df.REVIGO_BP <- cbind("Category"=cat_BP, df.REVIGO_BP)
  #head(df.REVIGO_BP)
  
  # read in REVIGO_CC
  df.REVIGO_CC <- read.csv(REVIGO_CC, header = T, stringsAsFactors = F)
  if (ncol(df.REVIGO_CC) == 11) {
    revigo.names = revigo.names_11
  } else {
    revigo.names = revigo.names_12
  }
  colnames(df.REVIGO_CC) <- revigo.names    
  cat_CC <- rep("CC", nrow(df.REVIGO_CC))
  df.REVIGO_CC <- cbind("Category"=cat_CC, df.REVIGO_CC)
  #head(df.REVIGO_CC)
  
  # read in REVIGO_MF
  df.REVIGO_MF <- read.csv(REVIGO_MF, header = T, stringsAsFactors = F)
  if (ncol(df.REVIGO_MF) == 11) {
    revigo.names = revigo.names_11
  } else {
    revigo.names = revigo.names_12
  }
  colnames(df.REVIGO_MF) <- revigo.names    
  cat_MF <- rep("MF", nrow(df.REVIGO_MF))
  df.REVIGO_MF <- cbind("Category"=cat_MF, df.REVIGO_MF)
  #head(df.REVIGO_MF)
  
  revigo.data <- rbind(df.REVIGO_BP, df.REVIGO_CC, df.REVIGO_MF)
  
  revigo.data <- revigo.data [(revigo.data$plot_X != "null" & revigo.data$plot_Y != "null"), ];
  revigo.data$plot_X <- as.numeric( as.character(revigo.data$plot_X) );
  revigo.data$plot_Y <- as.numeric( as.character(revigo.data$plot_Y) );
  revigo.data$plot_size <- as.numeric( as.character(revigo.data$plot_size) );
  revigo.data$log10_p_value <- as.numeric( as.character(revigo.data$log10_p_value) );
  if (ncol(revigo.data) == 12) {
    revigo.data$gene_counts <- as.numeric( as.character(revigo.data$gene_counts));
  }
  #revigo.data$frequency <- as.numeric( as.character(revigo.data$frequency) );
  revigo.data$uniqueness <- as.numeric( as.character(revigo.data$uniqueness) );
  revigo.data$dispensability <- as.numeric( as.character(revigo.data$dispensability) );
  #head(revigo.data);
  
  return(revigo.data)
}


ggRevigo <- function(revigo.data) {
  
  library( ggplot2 );
  library( scales );
  
  # --------------------------------------------------------------------------
  # Names of the axes, sizes of the numbers and letters, names of the columns,
  # etc. can be changed below
  
  p1 <- ggplot( data = revigo.data );
  if (ncol(revigo.data) == 13) {
    p1 <- p1 + geom_point( aes( plot_X, plot_Y, 
                                colour = log10_p_value, 
                                size = gene_counts), 
                           alpha = I(0.6) ) +
      scale_size_area();
  } else {
    p1 <- p1 + geom_point( aes( plot_X, plot_Y, 
                                colour = log10_p_value, 
                                size = plot_size), 
                           alpha = I(0.6) ) + 
      scale_size_area();
  }
  
  # scale color for log_p_value
  p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), 
                                     limits = c( min(revigo.data$log10_p_value), 0) );
  
  # add black circles to points
  if (ncol(revigo.data) == 13) {
    p1 <- p1 + geom_point( aes(plot_X, plot_Y, 
                               size = gene_counts), 
                           shape = 21, 
                           fill = "transparent", 
                           colour = I (alpha ("black", 0.6) )) + 
      scale_size_area();
  } else {
    p1 <- p1 + geom_point( aes(plot_X, plot_Y, 
                               size = plot_size), 
                           shape = 21, 
                           fill = "transparent", 
                           colour = I (alpha ("black", 0.6) )) + 
      scale_size_area();
  }
  
  p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  
  # add GO term descrition to plot, make dispensability less than 0.15 to black, others to gray
  black <- revigo.data [ revigo.data$dispensability < 0.15, ]; 
  gray <- revigo.data [ revigo.data$dispensability >= 0.15, ]; 
  p1 <- p1 + geom_text( data = black, 
                        aes(plot_X, plot_Y, label = description), 
                        colour = I(alpha("black", 0.85)), size = 3 );
  p1 <- p1 + geom_text( data = gray, 
                        aes(plot_X, plot_Y, label = description), 
                        colour = I(alpha("grey50", 0.85)), size = 3 );
  
  
  # add xlab and ylab
  p1 <- p1 + labs (y = "semantic space y", x = "semantic space x");
  
  # clean legend, remove key
  p1 <- p1 + 
    theme(#panel.border = element_blank(),
      panel.border = element_rect(fill=NA, colour = "black", size=1), 
      legend.key = element_blank(), 
      panel.grid.major = element_line(colour = "grey80", size=0.2), 
      panel.grid.minor = element_line(colour = "grey80", linetype = "dashed"), 
      axis.text=element_text(size=12, face="bold"),
      axis.title = element_text(size = rel(1), face="bold"),
      axis.title.y = element_text(angle = 90),
      axis.title.x = element_text(angle = 0),
      axis.ticks = element_line(size = 0.5), 
      axis.ticks.length=unit(0.2,'cm'), 
      plot.margin = unit(c(1,1,1,1), "cm")
    ) ;
  
  # set xlim and ylim
  revigo.x_range = max(revigo.data$plot_X) - min(revigo.data$plot_X);
  revigo.y_range = max(revigo.data$plot_Y) - min(revigo.data$plot_Y);
  p1 <- p1 + xlim(min(revigo.data$plot_X)-revigo.x_range/10,max(revigo.data$plot_X)+revigo.x_range/10);
  p1 <- p1 + ylim(min(revigo.data$plot_Y)-revigo.y_range/10,max(revigo.data$plot_Y)+revigo.y_range/10);
  
  # facet
  p1 <- p1 + facet_grid(~ Category);
  p1  
}


# read in BP, CC and MF data, return a revigo.data dataframe
revigo.data <- read_REVIGO_results(REVIGO_BP, REVIGO_CC, REVIGO_MF)
#revigo.data

# plot revigo.data
p <- ggRevigo(revigo.data)
#p

# save figure
# The file type depends on the extension (default=pdf).
ggsave(filename=out_pdf, plot=p, width=29.7, height=21.0, units="cm");

