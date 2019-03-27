#!/usr/bin/env Rscript

# load packages
suppressPackageStartupMessages({
    require("argparse")
    require("ape")
    require("ggplot2")
    require("ggtree")
    require("seqinr")
})

# parse arguments
parser <- ArgumentParser(description='plot sRNA tree and alignment')
parser$add_argument('nwk_tree', type="character", help='input tree in newick format')
parser$add_argument('fasta_alignment', type="character", help='sequence alignment in fasta format')
parser$add_argument('-b', '--bootstrap', type="double", default=0.1, help='show bootstrap values with this minimum threshold [default \"%(default)s\"]')
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
bt_threshold <- args$bootstrap
output_dir <- args$output_dir
output_file <- paste(args$output_prefix, ".pdf", sep="")
dir.create(output_dir, showWarnings = FALSE)

# get alignment length 
alignment <- seqinr::read.fasta(input_aln)
aln_length <- length(alignment[[1]])
nr_seqs <- length(alignment)

# msaplot using ggtree
tree <- read.tree(input_tree)
#p = ggtree(tree, layout = 'circular')
gt = ggtree(tree, ladderize = TRUE) + geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > bt_threshold), size=3)
# nucl color order (gap, 'A', 'C', 'G', 'T')
p <- msaplot(p=gt, 
             fasta=input_aln, 
             width=2, 
             offset = .0001 * aln_length, 
             bg_line=TRUE, 
             height=0.4,
             color = c("white", "royalblue2", "salmon", "goldenrod2", "green3"))
p + geom_tippoint(size=0.5) + geom_tiplab(size = 3, align = TRUE, linesize=.5, hjust=-0.1, vjust=-(1+1/nr_seqs)) #+ ggtitle("test_sRNA1")

# save plot
ggsave(filename=output_file, plot = last_plot(), path=output_dir, 
       scale = 1, width = 297, height = 210, units = "mm",
       dpi = 300, limitsize = FALSE)

# remove Rplots.pdf opened by ggplot2
# see https://stackoverflow.com/questions/17348359/how-to-stop-r-from-creating-empty-rplots-pdf-file-when-using-ggsave-and-rscript
unlink("Rplots.pdf")
