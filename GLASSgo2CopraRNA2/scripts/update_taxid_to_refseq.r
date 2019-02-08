#!/usr/bin/env Rscript

# ----------------------------------------------
# Download and read in prokaryotic genome report
# ----------------------------------------------

cat("=> Fetching prokaryotic genome report from NCBI ...\n")
if (file.exists("prokaryotes.txt")){
  file.remove("prokaryotes.txt")
}
system("wget -q ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt")
cat("Done!\n")

cat("=> Reading prokaryotes.txt ...\n")
genome_report <- read.table(file="prokaryotes.txt", header=FALSE, fill = TRUE, stringsAsFactors = FALSE,
                            sep = "\t", row.names = NULL, quote = "")
col_names <- c("Organism_Name", "TaxID", "BioProject_Accession", "BioProject_ID", "Group",	
               "SubGroup", "Size_in_Mb", "GC_Content", "Replicons", "WGS",	"Scaffolds", 
               "Genes", "Proteins", "Release_Date", "Modify Date", "Status",	"Center", 
               "BioSample_Accession", "Assembly_Accession", "Reference", "FTP Path", 
               "Pubmed_ID", "Strain")
colnames(genome_report) <- col_names
cat("Done!\n")


# ----------------------------------------------
# Download and read in prokaryotic genome report
# ----------------------------------------------

cat("=> Preparing taxid_to_refseq lookup table ...\n")

# get refseq for each taxid, return a dataframe
get_row <- function(genome_report_row){
  taxid <- genome_report_row$TaxID
  replicon_str <- genome_report_row$Replicons
  full_genome_entry <- strsplit(replicon_str, ";")[[1]][1]
  if (grepl("chromosome", full_genome_entry)){
    full_accession <- strsplit(full_genome_entry, ":")[[1]][2]
    Chromosomes.RefSeq <- strsplit(full_accession, "\\/")[[1]][1]
  } else {
    Chromosomes.RefSeq <- "-"
  }  
  return(data.frame(t(c(taxid, full_genome_entry, Chromosomes.RefSeq))))
}

# collect all the taxid2refseq dataframes into a row_list
row_list <- as.list(seq_len(nrow(genome_report)))
for (i in 1:nrow(genome_report)) {
  taxid <- genome_report[i, "TaxID"]
  replicon_str <- genome_report[i, "Replicons"]
  full_genome_entry <- strsplit(replicon_str, ";")[[1]][1]
  if (grepl("chromosome", full_genome_entry)){
    full_accession <- strsplit(full_genome_entry, ":")[[1]][2]
    Chromosomes.RefSeq <- strsplit(full_accession, "\\/")[[1]][1]
  } else {
    Chromosomes.RefSeq <- "-"
  }
  df <- data.frame(t(c(taxid, full_genome_entry, Chromosomes.RefSeq)), stringsAsFactors = FALSE)
  row_list[[i]] <- df
}

# merge all the taxid2refseq dataframes
ref <- do.call(rbind, row_list)
colnames(ref) <- c("TaxID", "full_genome_entry", "Chromosomes.RefSeq")
#head(ref)

save(ref, file = "taxid_to_refseq")
cat("Done!\n")
