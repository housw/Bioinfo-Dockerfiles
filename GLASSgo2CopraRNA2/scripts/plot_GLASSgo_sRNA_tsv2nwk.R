#!/usr/bin/env Rscript

# load packages
suppressPackageStartupMessages({
    require("argparse")
    require("ape")
    require("ggplot2")
    require("ggtree")
    require("grid")
})

# parse arguments
parser <- ArgumentParser(description='plot GLAssgo sRNA count, length and gc-content')
parser$add_argument('nwk_tree_file', type="character", help='input tree in newick format')
parser$add_argument('tsv_value_file', type="character", help='input sRNA associated values in tsv format')
parser$add_argument('-p', '--output_prefix', type="character", default='GLASSgo_tsv2nwk_plot', help='output prefix [default \"%(default)s\"]')
parser$add_argument('-o', '--output_dir', type="character", default="./", help='output directory [default \"%(default)s\"]')

# parse args
args <- parser$parse_args()
if (length(args) < 2) {
  stop("ERROR:not enough arguments are provided!", call.=FALSE)
  parser$print_help()
}
input_tree_file <- args$nwk_tree_file
input_value_file <- args$tsv_value_file
output_dir <- args$output_dir
output_file <- paste(args$output_prefix, ".pdf", sep="")
dir.create(output_dir, showWarnings = FALSE)

# read in tree
tree <- read.tree(input_tree_file)

# read in data
df <- read.table(input_value_file, header = TRUE, sep = '\t', as.is = TRUE, stringsAsFactors = TRUE)
colnames(df) <- c("id", "Count", "Identity", "Length", "GC_Content", "class", "order", "family", "genus", "species")
rank_cols <- c("class", "order", "family", "genus", "species") 
df[,rank_cols] <- lapply(df[,rank_cols], factor)


# visualize using ggtree
p <- ggtree(tree) + theme_tree2() #+ xlim(NA, 21) # this xlim changes all facet plots
p <- p %<+% df + geom_tiplab(size=1, offset=0.5, align=TRUE, linesize=.2, hjust=-0.1, vjust=-1) +
                 geom_tippoint(aes(size=Count, color=Identity), alpha=0.5)
data_df <- data.frame(id=tree$tip.label, data=df[df$id%in%tree$tip.label, c('Length', 'GC_Content', 'Identity')])
# https://guangchuangyu.github.io/2016/10/xlim_tree-set-x-axis-limits-for-only-tree-panel/
p2 <- facet_plot(p + xlim_tree(18), panel="GC Content", data=data_df, geom=geom_point, aes(x=GC_Content, color=Identity), alpha=.5)
p3 <- facet_plot(p2, panel="Length", data=data_df, geom=geom_segment, aes(x=0, xend=Length, y=y, yend=y, color=Identity), size=1, alpha=.5)
p4 <- p3 + theme(legend.position = c(0.1,0.7)) + scale_colour_gradient(low = "blue", high = "red") 

# grid draw
gt = ggplot_gtable(ggplot_build(p4))
gt$widths[5] = 4*gt$widths[5] # increase the width of tree panel
#grid.draw(gt)

# save plot
ggsave(filename=output_file, plot = gt, path=output_dir, 
       scale = 1, width = 297, height = 210, units = "mm",
       dpi = 300, limitsize = FALSE)

# remove Rplots.pdf opened by ggplot2
# see https://stackoverflow.com/questions/17348359/how-to-stop-r-from-creating-empty-rplots-pdf-file-when-using-ggsave-and-rscript
unlink("Rplots.pdf")
