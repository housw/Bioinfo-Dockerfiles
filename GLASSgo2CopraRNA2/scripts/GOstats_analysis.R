# GO enrichment analysis on CopraRNA targets
# 

#call: R --slave -f GOstats_analysis.R --args workdir=workdir organism=organism target_gene_file=target_gene_file ipr2go_file=ipr2go_file 

args <- commandArgs(trailingOnly = TRUE)
organism_name <- "Staphylococcus_HG001"
workdir <- "~/MEGA/GOstats/"
target_gene_file <- "HG001_00009__HG001_03232_CopraRNA_result_locustags.tsv"
ipr2go_file <- "NZ_CP018205_ipr2go.tsv"

for(i in 1:length(args)){
  temp<-strsplit(args[i],"=")
  temp<-temp[[1]]
  temp1<-temp[1]
  temp2<-temp[2]
  assign(as.character(temp1),temp2)
}

setwd(workdir)

# load library
library(GOstats)
library(GSEABase)
#library(qvalue)
library(dplyr)
library(tools)

# install AnnotationForge
#source("https://bioconductor.org/biocLite.R")
#biocLite("AnnotationForge")
#biocLite("GOstats")

# check available organism in the AnnotationForge package
#library("AnnotationForge")
#available.dbschemas()


################
#              #
#  FUNCTIONS   #
#              #
################
# prepare GeneSetCollection from ipr2go file
# modified from https://github.com/davfre/GOstatsPlus/blob/master/R/b2g-to-gsc.R

get_ipr2go <- function(ipr2go_file){
    ipr2go <- read.table(file=ipr2go_file, sep="\t", fill=TRUE)
    colnames(ipr2go)<-c("GeneID", "GO")
    return(ipr2go)
}


ipr2go_to_gsc <- function(file, organism="organism"){

    # read in ipr2go_file    
    ipr2go <- get_ipr2go(file)
    
    # Create GO package
    GOdata <- data.frame(GO.id = ipr2go$GO,
                       evidence.code = rep("ISA", dim(ipr2go)[1]),
                       gene.id = as.character(ipr2go$GeneID))
    goFrame <- GOFrame(GOdata, organism = organism)
    goAllFrame <- GOAllFrame(goFrame)
    
    GO_gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
    #  GO_mappings<-getGOFrameData(goAllFrame)
    
    return(GO_gsc)
}


test_GO<-function(genes, ontology, gsc, direction = "over", universe, pval = 0.01){
  Test<-NA
  
  # default to test against all genes with any annotation
  if(missing(universe)){
    universe<-unique(unlist(geneIds(gsc)))
  }
  
  if(length(genes[genes %in% universe])>0){
    params <- GSEAGOHyperGParams(name = "GSEA based annotation Parameters",
                                 geneSetCollection = gsc,
                                 geneIds = genes,
                                 universeGeneIds = universe,
                                 ontology = ontology,
                                 pvalueCutoff = pval,
                                 conditional = FALSE,
                                 testDirection = direction)
    Test <- hyperGTest(params)
  }
  return(Test)
}


significant_terms<-function(GO_OBJ, cutoff){
  terms_pvals = pvalues(GO_OBJ)
  return(terms_pvals[terms_pvals<cutoff])
}


write2file <- function(data_frame, file){
  write.table(data_frame, file=file, sep="\t", row.names = F, quote = F)
}

# wrapper of ipr2go_to_gsc, to load gsc file from rda data or prepare it from ipr2go file
load_gsc <- function(ipr2go_file, organism_name){
  rda_file <- paste(organism_name, "_gsc.rda", sep="")
  if(file.exists(rda_file)){
    load(file=rda_file)
  }else{
    GO_gsc = ipr2go_to_gsc(file = ipr2go_file, organism = organism_name)
    save(GO_gsc, file = rda_file)
  }
  return(GO_gsc)
}


# do pvalue_adjustment and convert to dataframe
add_qvalue <- function(GO){
  Pvalue <- pvalues(GO)
  df <- summary(GO)
  Qvalue <- p.adjust(Pvalue, method = "BH")
  #Qvalue <- qvalue(Pvalue)
  #head(Qvalue$qvalues)
  #Qvalue$qvalues
  #df$qvalue <- Qvalue$qvalues[df[, 1]]
  df$qvalue <- Qvalue[df[, 1]]
  return(df)
}


# ----------------------------------------------
# 1) prepare gsc object from ipr2go_file
# ----------------------------------------------
# make GeneSetCollection object from ipr2go_file
GOstats_gsc <- load_gsc(ipr2go_file = ipr2go_file, organism_name = organism_name)

# 2) read in target genes 
target_geneIDs <- read.table(file=target_gene_file, header = F, stringsAsFactors = F) 
target_geneIDs <- unique(as.vector(target_geneIDs$V1))

# 3) population or background geneID
ipr2go <- get_ipr2go(ipr2go_file)
gene_universe <- unique(ipr2go$GeneID)
  

# 4) perform go enrichment test
# when using only expressed genes as universe, one can do: 
# > GO_BP <- test_GO(gene_IDs_of_interest, ontology="BP", gsc=Ana315_gsc, pval=0.05, universe=gene_universe)
target_gene_file_basename <- file_path_sans_ext(target_gene_file)
output_basename <- paste(target_gene_file_basename, "_GO_enrichment", sep="")
GO_BP <- test_GO(target_geneIDs, ontology="BP", gsc=GOstats_gsc, pval=0.05)
GO_MF <- test_GO(target_geneIDs, ontology="MF", gsc=GOstats_gsc, pval=0.05)
GO_CC <- test_GO(target_geneIDs, ontology="CC", gsc=GOstats_gsc, pval=0.05)

# 5) add qvalue and write to csv and html
df.GO_BP <- add_qvalue(GO_BP)
htmlReport(r=GO_BP, file = paste(output_basename, "_BP.html", sep=""), summary.args=list("htmlLinks"=TRUE))
write2file(df.GO_BP, paste(output_basename, "_BP.tsv", sep=""))
df.GO_MF <- add_qvalue(GO_MF)
write2file(df.GO_MF, paste(output_basename, "_MF.tsv", sep=""))
htmlReport(r=GO_MF, file = paste(output_basename, "_MF.html", sep=""), summary.args=list("htmlLinks"=TRUE))
df.GO_CC <- add_qvalue(GO_CC)
write2file(df.GO_CC, paste(output_basename, "_CC.tsv", sep=""))
htmlReport(r=GO_CC, file = paste(output_basename, "_CC.html", sep=""), summary.args=list("htmlLinks"=TRUE))
