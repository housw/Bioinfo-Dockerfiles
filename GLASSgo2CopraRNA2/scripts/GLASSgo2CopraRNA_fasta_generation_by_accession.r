# author Jens Georg
# This script extracts all sequences which are available for CopraRNA and writes them in a fasta file
# Dependencies: "taxid_to_refseq" and "CopraRNA_available_organisms.txt"

#call: R --slave -f  GLASSgo_2_copraRNA_fasta_generation.r --args input_file=4083138.result refpath=taxid_to_refseq cop_path=copra_refseq_positivelist.txt output_file=coprarna_candidates.txt
args <- commandArgs(trailingOnly = TRUE)

input_file <- "GLASSgo_output.fa"
refpath <- "taxid_to_refseq"
cop_path <- "CopraRNA_available_organisms.txt"
output_file <- "coprarna_candidates.txt"

for(i in 1:length(args)){
  temp<-strsplit(args[i],"=")
  temp<-temp[[1]]
  argument<-temp[1]
  value<-temp[2]
  assign(as.character(argument), value)
}


to_refseq2<-function(result){
  load(refpath)
  result<-gsub("\\..*","",result)
  temp<-c()
  
  for(i in 1:length(result)){
    temp1<-grep(result[i],ref[,"full_genome_entry"])
    temp<-c(temp,temp1[1])
  }
  temp<-ref[temp,"Chromosomes.RefSeq"]
  refseq<-temp
  refseq
}

fill<-function(x){
  temp<-which(x[,ncol(x)]=="no_name")
  out<-x[,ncol(x)]
  count<-1
  while(length(temp)>0){
    out[temp]<-x[temp,(ncol(x)-count)]
    temp<-which(out=="no_name")
  }
  out
}

# rewrited to include accession number and coordinates
export_ncRNA_coordinates<-function(x){ # x = copy pasted text from GLASSgo fasta output
  #x <- fasta
  header_row <- grep(">", x)
  headers <- as.character(x[header_row])
  seqs<-as.character(x[header_row+1])
  
  # slice and clean strain_name from header
  strain_names <- headers
  strain_names <- gsub(".*[0-9]{1,}-[0-9]{1,} ","", strain_names)
  strain_names <- gsub(",.*","", strain_names)
  strain_names <- gsub("genome assembly","", strain_names)
  strain_names <- gsub("complete genome","", strain_names)
  strain_names <- gsub("\\sgenome-p.c.VAL.*","", strain_names)
  strain_names <- gsub("\\ssequence-p.c.VAL.*","", strain_names)
  strain_names <- gsub("\\s-p.c.VAL.*","", strain_names)
  strain_names <- gsub("\\sDNA.*","", strain_names)
  #strain_names
  
  # get locations contain start, end and strand infromation
  locations <- c()
  accessions <- c()

  for (i in 1:length(headers)){
    loc <- strsplit(headers, " ")[[i]][1] # strip genome information
    acc <- strsplit(loc, ":")[[1]][1] # get accession
    #acc <- strsplit(loc, "\\.")[[1]][1] # remove .1 suffix of accession
    loc <- strsplit(loc, ":")[[1]][2] # get location
    #cat(acc, loc, "\n")
    locations <- c(locations, loc)
    accessions <- c(accessions, gsub(">", "", acc))
  }
  #locations
  #accessions
  
  # summarize strand, start, end, name, header and sequence information
  out<-matrix(, length(locations), 7)
  colnames(out)<-c("Accesion_number", "Strand","start","end","name","Full_header","sequence")
  out[, 1] <- accessions
  out[, 2]<-"+"
  for(i in 1:length(locations)){
    temp<-grep("c",locations[i])
    if(length(temp)==1){
      out[i, 2]<-"-"
    }
    temp<-gsub("c","",locations[i])
    temp<-strsplit(temp,"-")
    out[i, 3] <- temp[[1]][1]
    out[i, 4] <- temp[[1]][2]
  }
  out[, 5] <- strain_names
  out[, 6] <- as.character(headers)
  out[, 7] <- as.character(seqs)
  out
}

fasta<-read.delim(input_file,header=F)[,1]

# for the query of GLASSgo, no taxID in the header, so skip the first two lines
direct<-grep("taxID:", fasta[1])
if(length(direct)==0){
  fasta<-fasta[-c(1,2)]
}

# remove sRNAs contain "N"s
nn<-grep(">",fasta)
nn1<-grep("n",fasta[nn+1], ignore.case=T)
if(length(nn1)>0){
  fasta<-fasta[-c(nn[nn1],nn[nn1]+1)]
}

# combine accession number and coordinates 
coor<-export_ncRNA_coordinates(fasta)
#head(coor)

# ----------------
# if more than 1 homolog is detected for one organism or the same Refseq ID, keep only the homolog with the highest identity to the input
# ----------------
# parse and add identity column
iden<-gsub(".*VAL:","",coor[,"Full_header"])
iden<-as.numeric(gsub("%.*","",iden))
coor<-cbind(coor,iden)
coor<-coor[order(iden, decreasing=T), ]

# find duplicates with the same accessions 
dup<-which(duplicated(coor[,1]))
if(length(dup)>0){
  coor<-coor[-dup,]
}
coor2<-coor
ac<-coor[,1]
coor[,1]<-paste(coor[,1],coor[,3], sep="_")

taxi<-gsub(".*taxID:","",coor[,"Full_header"])

# get refseq accessions
fin<-to_refseq2(ac)
coor<-cbind(coor,taxi,fin)

dup<-which(duplicated(fin))

if(length(dup)>0){
  coor<-coor[-dup,]
  coor2<-coor2[-dup,]
}

# remove rows with refseq == "-"
if(is.matrix(coor)==T){
  na<-which((coor[,"fin"])=="-")
  if(length(na)>0){
    coor<-coor[-na,]
  }
}

# remove rows with refseq == "NA"
if(is.matrix(coor)==T){
  na<-which(is.na(coor[,"fin"]))
  if(length(na)>0){
    coor<-coor[-na,]
  }
}

# keep only homologs represented in the coprarna reference file	
copref<-read.delim(cop_path, sep="\t", header=T,comment.char = "#")	
se<-function(x){
  out<-grep(x, copref[,1])[1]
  if(length(out)==0){
    out<-NA
  }
  out
}
notinlist<-which(is.na(unlist(lapply(gsub("\\..*","",coor[,"fin"]),se))))

if(length(notinlist)>0){
  coor<-coor[-notinlist,]	
}


fasta3<-c()
for(i in 1:nrow(coor)){
  fasta3<-c(fasta3, paste(">",coor[i,"fin"],"|",gsub(">","",coor[i,"Full_header"]),sep=""))
  fasta3<-c(fasta3, as.character(coor[i,"sequence"]))
}

write.table(fasta3, file=output_file, row.names=F, col.names=F, quote=F)

copref<-read.delim(cop_path, sep="\t", header=T,comment.char = "#")	
nam<-c()
for(i in 1:nrow(coor)){
  tnam<-grep(gsub("\\..*","",coor[i,"fin"]),copref[,1])
  nam<-c(nam,as.character(copref[tnam,2]))
}

nam2<-c()
for(i in 1:length(nam)){
	temp1<-substr(nam[i],1,3)
	temp2<-strsplit(nam[i],"_")[[1]]
	temp1<-paste(temp1,"_",temp2[2], sep="")
	if(length(temp2)>2){
		temp1<-paste(temp1, temp2[length(temp2)], sep="_")
	}
	nam2<-c(nam2,temp1)
}
nam2<-paste(nam2,coor[,"fin"], sep="_")

coor<-cbind(coor,nam2)
system("rm -f full_GLASSgo_table.Rdata")
save(coor, file="full_GLASSgo_table.Rdata")
