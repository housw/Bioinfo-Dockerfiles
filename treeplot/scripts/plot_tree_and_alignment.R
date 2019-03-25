#!/usr/bin/env Rscript

# load packages
suppressPackageStartupMessages({
    require("argparse")
    require("ape")
    require("ggplot2")
    require("ggtree")
})

# parse arguments
parser <- ArgumentParser(description='plot sRNA tree and alignment')
parser$add_argument('nwk_tree', type="character", help='input tree in newick format')
parser$add_argument('fasta_alignment', type="character", help='sequence alignment in fasta format')
parser$add_argument('-p', '--output_prefix', type="character", default='treeplot', help='output prefix [default \"%(default)s\"]')
parser$add_argument('-o', '--output_dir', type="character", default="./", help='output directory [default \"%(default)s\"]')

# parse args
args <- parser$parse_args()
if (length(args) < 2) {
  stop("ERROR:not enough arguments are provided!", call.=FALSE)
  parser$print_help()
}
input_tree <- args$nwk_tree
input_aln <- args$fasta_alignment
output_dir <- args$output_dir
output_file <- paste(args$output_prefix, ".pdf", sep="")
dir.create(output_dir, showWarnings = FALSE)

# msaplot using ggtree
tree <- read.tree(input_tree)
#p = ggtree(tree, layout = 'circular')
p <- msaplot(p=ggtree(tree), 
             fasta=input_aln, 
             width=2, 
             offset = .1, 
             bg_line=TRUE, 
             height=0.5,
             color = c("grey90", "blue", "red", "yellow", "green"))
p + geom_tippoint(size=0.5) + geom_tiplab(size = 3, align = TRUE, linesize=.5) #+ ggtitle("test_sRNA1")

# save plot
ggsave(filename=output_file, plot = last_plot(), path=output_dir, 
       scale = 1, width = 50, height = 50, units = "cm",
       dpi = 300, limitsize = FALSE)
