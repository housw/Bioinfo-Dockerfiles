#call: R --slave -f  synteny_pdf.r --args input_sRNA=GLASSgo_output.fa outfile_prefix=sRNA_synteny


rang<-function(s,e,strand,windo=3000, anno){
  temp<-which(as.numeric(anno[,2])>=max(0,(s-windo)) & as.numeric(anno[,3])<=(e+windo))
  
  temp2<-which(as.numeric(anno[,2])<=max(0,(s-windo)) & as.numeric(anno[,3])>=(s-windo))
  temp3<-which(as.numeric(anno[,2])<=max(0,(e+windo)) & as.numeric(anno[,3])>=(e+windo))
  temp<-c(temp, temp2, temp3)
  temp<-anno[temp,]
  temp
  
  
}


export_accesions<-function(x){ # x = copy pasted text from GLASSgo fasta output
  temp<-grep(">", x)
  x<-x[temp]
  x<-strsplit(x,"\\|")
  out<-c()
  for(i in 1:length(x)){
    if(length(x[[i]])>2){
      #out<-c(out,x[[i]][2]) # for(RNAlien)
	  out<-c(out,x[[i]][4]) # forGLASgo
    }
  }
  if(is.matrix(out)==F){
	out<-t(as.matrix(out))
  }
 # write.table(out, file="Accesions.txt", col.names = F, row.names = F, quote=F) # list of all Accesion numbers
 
  out<-matrix(out,length(out),1)
  colnames(out)<-"Accesion_number"
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



generate_feattable<-function(idvector, genomepath=FALSE, featpath=FALSE){
  #idvector <- idvect[i]
  #idvector
  #genomepath=FALSE
  wd<-getwd()
	if(featpath==FALSE){
		dir.create("feattables")
	}
	
	if(genomepath!=FALSE){
		setwd(genomepath)
	}
	
	
	#di1<-paste(wd,"/feattables", sep="") #path precomputed feattables
	#gbk<-#"path_genomefiles"
	#gbk1<-dir(gbk)
	#di1<-dir(di1)
	#di<-dir()
	#x<-as.data.frame(x)
	#ids<-unique(x[,1])
	#ex<-match(di1,ids)
	#ex<-na.omit(ex)
	#naex<-ids[ex]
	#if(length(ex)>0){
	#	ids<-ids[-ex]
	#}
	ids<-idvector
  #ids
	if(length(ids)>0){
		for(i in 1:length(ids)){
		  #i <- 1
			if(length(ids)>1){
				print(paste(i,"/",length(ids),sep=""))
			}
			#name1<-ids[i]
			
			#temp<-which(di==ids[i])
			#temp<-paste(gbk,"/",temp, sep="")
			temp_full<-read.delim(as.character(ids)[i])[,1]
			temp_full<-gsub("\\\n *","",temp_full)
			head(temp_full)
			anfang<-grep("FEATURES   ",temp_full)
			ende<-grep("ORIGIN   ",temp_full)
			name<-grep("VERSION  ",temp_full)
			for(jj in 1:length(anfang)){
			temp<-temp_full[anfang[jj]:ende[jj]]
			tax<-grep("taxon:", temp)
			tax<-gsub(".*taxon:","",temp[tax])
			name1<-temp_full[name[jj]]
			name1<-gsub("VERSION *","",name1)
			name1<-strsplit(name1," ")
			name1<-name1[[1]][1]
			#temp_gene1<-grep("gene +[0123456789]+..",temp)
			#temp_gene2<-grep("gene +complement\\([0123456789]+..",temp)
			#temp_gene<-c(temp_gene1, temp_gene2)
			temp_gene<-grep("gene   ",temp)
			temp_cds<-grep("CDS   ",temp)
			
			
			if(length(temp_gene)>0){
				#temp_cds<-grep("CDS ", temp)
				accession_info<-matrix(,length(temp_gene),7)
				colnames(accession_info)<-c("strand","start","end","gene_name","locus_tag","product","AA_sequence")
				accession_info[,1]<-"+"
				for(j in 1:length(temp_gene)){
					#print(paste(j,"/",length(temp_gene),sep=""))
					if(j<length(temp_gene)){
						temp_y<-temp[temp_gene[j]:(temp_gene[j+1]-1)]
					}
					if(j==length(temp_gene)){
						temp_y<-temp[temp_gene[j]:(length(temp)-1)]
					}
					comp<-grep("complement",temp_y[1])
					if(length(comp)==0){
					
					coor<-gsub(" ","",temp_y[1])
					coor<-gsub("gene","",coor)
					coor1<-coor
					coor<-strsplit(coor,"\\.\\.")
					s<-as.numeric(coor[[1]][1])
					e<-as.numeric(coor[[1]][2])
					if(is.na(s)==FALSE & is.na(e)==FALSE){
						accession_info[j,2]<-coor[[1]][1]
						accession_info[j,3]<-coor[[1]][2]
					}
					}
					if(length(comp)==1){
					coor1<-gsub("gene","",temp_y[1])
					coor1<-gsub(" ","",coor1)
					coor<-gsub("gene.*complement\\(","",temp_y[1])
					coor<-gsub("\\)","",coor)
					coor<-gsub(" ","",coor)
					
					coor<-strsplit(coor,"\\.\\.")
					s<-as.numeric(coor[[1]][1])
					e<-as.numeric(coor[[1]][2])
					if(is.na(s)==FALSE & is.na(e)==FALSE){
						accession_info[j,2]<-coor[[1]][1]
						accession_info[j,3]<-coor[[1]][2]
						accession_info[j,1]<-"-"
					}
					}
					tempg<-grep("/gene=",temp_y)
					if(length(tempg)>0){
						tempg<-gsub("/gene=","",temp_y[tempg[1]])
						tempg<-gsub(" ","",tempg)
						accession_info[j,4]<-tempg
					}
					templ<-grep("/locus_tag=",temp_y)
					if(length(temp)>0){
						templ<-gsub("/locus_tag=","",temp_y[templ[1]])
						templ<-gsub(" ","",templ)
						accession_info[j,5]<-templ
					}
					cds<-grep("CDS", temp_y)
					if(length(cds)==1){
						coor2<-gsub("CDS","",temp_y[cds])
						coor2<-gsub(" ","",coor2)
						if(coor2==coor1){
							tempp<-grep("/product=", temp_y)
							if(length(tempp)==1){
								tempp<-gsub(" */product=","",temp_y[tempp])
								#tempp<-gsub("\\\n *","",tempp)
								accession_info[j,6]<-tempp
							}
							tempp<-grep("/translation=", temp_y)
							if(length(tempp)==1){
								tempp<-gsub(" */translation=","",temp_y[tempp])
								#tempp<-gsub("\\\n *","",tempp)
								accession_info[j,7]<-tempp
							}
						}
					}
					
					
				}
			
		feat<-accession_info
		if(length(tax)==0){
				tax<-"no_tax_info"
		}
		feat<-list(feat,tax)
        #names(feat)<-name1
        if(featpath!=FALSE){
			#setwd(featpath)
			#dir.create("feattables")
			di<-paste(featpath,"/",name1, sep="")
			print(di)
		}
        if(featpath==FALSE){
			di<-paste(wd,"/feattables/",name1, sep="")
			print(di)
		}
		
		
		nas1<-which(is.na(feat[[1]][,5]))
		feat[[1]][nas1,5]<-paste("orf_",seq(1,length(nas1)),sep="")
        save(feat, file=di)
		if(genomepath!=FALSE){
			#setwd(genomepath)
		}
		}
		
		
		
		if(length(temp_gene)==0){
				temp_gene<-temp_cds
				if(length(temp_gene)>0){
				#temp_cds<-grep("CDS ", temp)
				accession_info<-matrix(,length(temp_gene),7)
				colnames(accession_info)<-c("strand","start","end","gene_name","locus_tag","product","AA_sequence")
				accession_info[,1]<-"+"
				for(j in 1:length(temp_gene)){
					#print(paste(j,"/",length(temp_gene),sep=""))
					if(j<length(temp_gene)){
						temp_y<-temp[temp_gene[j]:(temp_gene[j+1]-1)]
					}
					if(j==length(temp_gene)){
						temp_y<-temp[temp_gene[j]:(length(temp)-1)]
					}
					comp<-grep("complement",temp_y[1])
					if(length(comp)==0){
					
					coor<-gsub(" ","",temp_y[1])
					coor<-gsub("CDS","",coor)
					coor1<-coor
					coor<-strsplit(coor,"\\.\\.")
					s<-as.numeric(coor[[1]][1])
					e<-as.numeric(coor[[1]][2])
					if(is.na(s)==FALSE & is.na(e)==FALSE){
						accession_info[j,2]<-coor[[1]][1]
						accession_info[j,3]<-coor[[1]][2]
					}
					}
					if(length(comp)==1){
					coor1<-gsub("CDS","",temp_y[1])
					coor1<-gsub(" ","",coor1)
					coor<-gsub("CDS.*complement\\(","",temp_y[1])
					coor<-gsub("\\)","",coor)
					coor<-gsub(" ","",coor)
					
					coor<-strsplit(coor,"\\.\\.")
					s<-as.numeric(coor[[1]][1])
					e<-as.numeric(coor[[1]][2])
					if(is.na(s)==FALSE & is.na(e)==FALSE){
						accession_info[j,2]<-coor[[1]][1]
						accession_info[j,3]<-coor[[1]][2]
						accession_info[j,1]<-"-"
					}
					}
					tempg<-grep("/gene=",temp_y)
					if(length(tempg)>0){
						tempg<-gsub("/gene=","",temp_y[tempg[1]])
						tempg<-gsub(" ","",tempg)
						accession_info[j,4]<-tempg
					}
					templ<-grep("/locus_tag=",temp_y)
					if(length(temp)>0){
						templ<-gsub("/locus_tag=","",temp_y[templ[1]])
						templ<-gsub(" ","",templ)
						accession_info[j,5]<-templ
					}
					cds<-grep("CDS", temp_y)
					
							tempp<-grep("/product=", temp_y)
							if(length(tempp)==1){
								tempp<-gsub(" */product=","",temp_y[tempp])
								#tempp<-gsub("\\\n *","",tempp)
								accession_info[j,6]<-tempp
							}
							tempp<-grep("/translation=", temp_y)
							if(length(tempp)==1){
								tempp<-gsub(" */translation=","",temp_y[tempp])
								#tempp<-gsub("\\\n *","",tempp)
								accession_info[j,7]<-tempp
							}
						
					
					
					
				}
			
		feat<-accession_info
		
		
		
		
		
		if(length(tax)==0){
				tax<-"no_tax_info"
		}
		feat<-list(feat,tax)
        #names(feat)<-name1
        if(featpath!=FALSE){
			#setwd(featpath)
			#dir.create("feattables")
			di<-paste(featpath,"/",name1, sep="")
			print(di)
		}
        if(featpath==FALSE){
			di<-paste(wd,"/feattables/",name1, sep="")
			print(di)
		}
		
		
		
		nas1<-which(is.na(feat[[1]][,5]))
		feat[[1]][nas1,5]<-paste("orf_",seq(1,length(nas1)),sep="")
        save(feat, file=di)
		if(genomepath!=FALSE){
			#setwd(genomepath)
		}
		}
				
				
			}
		
		
		
		if(length(temp_gene)==0){
			if(length(tax)==0){
				tax<-"no_tax_info"
			}
			feat<-"no_info"
			feat<-list(feat,tax)
		
		if(featpath!=FALSE){
			#setwd(featpath)
			#dir.create("feattables")
			di<-paste(featpath,"/",name1, sep="")
			print(di)
		}
        if(featpath==FALSE){
			di<-paste(wd,"/feattables/",name1, sep="")
			print(di)
		}
		
		
		
		
        save(feat, file=di)
		if(genomepath!=FALSE){
			#setwd(genomepath)
		}
		}
		#setwd(wd)
		}
		}
		
	}
	setwd(wd)
}

generate_feattable_fetch<-function(idvect,featpath=FALSE){
	#wd<-getwd()
  #idvect <- ids
	require(seqinr)
	for(i in 1:length(idvect)){
		print(paste(i,"/",length(idvect)))
		gen <- entrez_fetch(db="nucleotide", id=idvect[i],rettype="gb",retmode="text")
		write.table(gen, file=as.character(idvect[i]), sep="\t", quote=FALSE, row.names=F, col.names=F)
		generate_feattable(idvect[i], featpath=featpath)
		unlink(idvect[i])
	}
	#setwd(wd)
}

readFeatures2<-function(coor, featpath=FALSE){ #x = coor object
  #coor <- coor2
  #featpath="/mnt/data/db/GLASSgo/test_synteny"
	require(seqinr)
	wd<-getwd()
	 di1<-paste(wd,"/feattables", sep="") #path precomputed feattables
	 if(featpath!=FALSE){
		di1<-featpath
	 }
	di2<-dir(di1)
	x<-as.data.frame(coor)
	head(x)
	ids<-unique(x[,1])
	ex<-match(di2,ids)
	ex<-na.omit(ex)
	naex<-ids[ex]
	if(length(ex)>0){
		ids<-ids[-ex]
	}
	na<-naex
	taxids<-rep(NA, nrow(coor))
	if(length(ids)>0){
		print(ids)
		generate_feattable_fetch(ids, featpath=featpath)
		#di1<-paste(wd,"/feattables", sep="") #path precomputed feattables
		if(featpath!=FALSE){
		d1<-featpath
		}
		di2<-dir(di1)
		ids<-unique(x[,1])
		ex<-match(di2,ids)
		ex<-na.omit(ex)
		naex<-ids[ex]
		na<-naex
	}
	na<-as.character(na)
	cou<-0
	genomes<-list()
	for(i in 1:length(na)){
		test<-which(di2==na[i]) 
		if(length(test)>0){
			test<-test[1]
			di<-paste(wd,"/feattables/",na[i],sep="")
			if(featpath!=FALSE){
				di<-paste(featpath,"/",na[i],sep="")
			}
			load(di)
			tax<-feat[[2]]
			pos<-which(coor[,1]==na[i])
			taxids[pos]<-tax[1]
			if(feat[1]!="no_info"){
				genomes[[length(genomes)+1]]<-feat[[1]]
				names(genomes)[length(genomes)]<-as.character(na[i])
				print( names(genomes)[length(genomes)])
			}
		}
	}
	genomes<-list(genomes,taxids)
}





data_preparation<-function(x, coor){ # x = output of readFeatures, coor = output of wrapper 
  require(rentrez)
  out<-list()
  exist<-matrix(,nrow(coor),2)
  colnames(exist)<-c("GLASSgo","genomes")
  for(i in 1:nrow(coor)){
    temp<-grep(coor[i,1],names(x))
    #print(paste(i, temp, sep="_"))
    if(length(temp)>0){
      exist[i,1]<-i
      exist[i,2]<-temp[1]
      
    }
  }
  if(length(which(is.na(exist[,1])))>0){	
    print("d")
    exist<-exist[-which(is.na(exist[,1])),]	
  }		
  print(paste(length(names(x)),nrow(coor),nrow(exist), sep="/"))
  man<-c()
  mas<-c()		
  
  for(i in 1:nrow(exist)){
    #for(i in 1:45){  
   # print(i)
    temp<-coor[exist[i,1],]
    stra<--1
    if(temp[2]=="+"){
      stra<-1
    }
    na<-as.character(temp[5])
    nam<-paste(temp[1],temp[3], sep="_")
    s<-as.numeric(temp[3])
    e<-as.numeric(temp[4])
    anno<-x[[exist[i,2]]]
    temp_out<-rang(s,e, anno=anno)
    if(length(dim(temp_out))<2){
      temp_out<-t(as.matrix(temp_out))
      
    }
    if(nrow(temp_out)==0){
      g<-c(stra,s,e,"NA","NA","NA","NA")
      temp_out<-rbind(temp_out,g)
      
    }
    if(nrow(temp_out)>0){
      
      
      ma<-max(as.numeric(temp_out[,3]))
      mi<-min(as.numeric(temp_out[,2]))
      aa<-as.numeric(temp_out[,2])-mi
      #print(aa)
      bb<-as.numeric(temp_out[,3])-mi
      s_srna<-min(as.numeric(temp[3:4]))-mi
      e_srna<-max(as.numeric(temp[3:4]))-mi
      temp_out<-cbind(temp_out,aa,bb,rep(s_srna,nrow(temp_out)),rep(e_srna,nrow(temp_out)),rep(stra,nrow(temp_out)),rep(na,nrow(temp_out)))
	  
	  orf<-grep("orf_", temp_out[,"locus_tag"])
	  if(length(orf)>0){
			temp_out[orf,"locus_tag"]<-paste(temp[1],temp_out[orf,"locus_tag"], sep="_")
	  }
	  
	  
      d<-ma-mi
      man<-c(man,d)
      mas<-c(mas,s_srna)
      out[[length(out)+1]]<-temp_out
      names(out)[length(out)]<-nam
      #print(c(i, length(out)))
    }
    
  }		
  
  d<-max(d)
  mas<-max(mas)	
  
  
  s_vect<-c()
  #for(i in 1:length(out)){
  #print(i)
  #	temp<-out[[i]]
  #	s<-mas-as.numeric(temp[1,9])
  #	s_vect<-c(s_vect,s)
  #	temp[,7]<-as.numeric(temp[,7])+s
  #	temp[,8]<-as.numeric(temp[,8])+s
  #	temp[,9]<-as.numeric(temp[,9])+s
  #	temp[,10]<-as.numeric(temp[,10])+s
  #	out[[i]]<-temp
  #}
  
  out<-list(out, exist,  d, coor[exist[,1],], s_vect)	
  out	
}


get_prot_fasta<-function(out){# outpit from data_preparation function first list entry
  require(rentrez)
  fasta<-c()
  for(i in 1:length(out)){
    for(j in 1:nrow(out[[i]])){
      #print(j)
      if(is.na(out[[i]][j,5])==F){
        temp<-entrez_search(db="protein",term=out[[i]][j,5])
        #print(temp)
        if(temp[[1]][1]!="NULL"){
          
          temp<-entrez_fetch(db="protein", id=temp[[1]][1], rettype="fasta",retmode="text")
          temp<-unlist(strsplit(temp,"\n"))
          temp[1]<-paste(">",out[[i]][j,5], sep="" )
          #print(out[[i]][j,5])
          #print(temp[1])
          fasta<-c(fasta,temp)
        }
      }
    }
  }
  write.table(fasta, file="protein_fasta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  
}

locus_tag2org<-function(out2){
	tag<-c()
	org<-c()
	for(i in 1:length(out2[[1]])){
		temp_tag<-out2[[1]][[i]][,5]
		temp_org<-rep(names(out2[[1]])[i], length(temp_tag))
		tag<-c(tag,temp_tag)
		org<-c(org,temp_org)
	}
	out<-cbind(tag,org)
	out
}


get_prot_fasta3<-function(out){# outpit from data_preparation function first list entry
  
  fasta<-c()
  for(i in 1:length(out)){
    for(j in 1:nrow(out[[i]])){
      #print(j)
      if(is.na(out[[i]][j,7])==F){
        #temp<-entrez_search(db="protein",term=out[[i]][j,5])
        #print(temp)
        temp<-out[[i]][j,7]
         #print(out[[i]][j,5])
         #print(temp[1])
		 na<-paste(">",out[[i]][j,5],sep="")
		 temp<-c(na,temp)
         fasta<-c(fasta,temp)
        }
      }
    }
  
  write.table(fasta, file="protein_fasta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  
}



cdhit_run<-function(fasta="protein_fasta.txt", outname="psi", thres=0.3, psi=T){
  wd<-getwd()
  di<-paste(wd,"/", "psi_out",sep="")
  dir.create(di)
  if(psi==T){
    inp<-paste("./psi-cd-hit.pl -i ", fasta,  " -o ", di,"/",outname, " -c ", thres, sep="")
  }
  if(psi==F){
    inp<-paste("cd-hit -i ", fasta,  " -o ", di,"/",outname, " -c ",  thres ," -n 2", " -aL 0.6", sep="")
  }
  print(inp)
  system(inp)
  cd<-paste(di, "/", outname, ".clstr", sep="")
  cd<-read.delim(cd, header=F, sep="?")
  cd<-as.character(cd[,1])
  cd<-gsub("\t"," ", cd)
  cd
}

proc_cdhit<-function(x){ #x= pasted sorted cdhit cluster file
  clustlist<-list()
  numb<-grep(">Cluster", x)
  
  for(i in 1:length(numb)){
    if(i<length(numb)){
      end<-numb[i+1]-1
    }
    if(i==length(numb)){
      end<-length(x)
    }
    
    temp<-x[(numb[i]+1):end]
    temp<-gsub(".*aa, >","",temp)
    temp<-gsub("\\.\\.\\..*","",temp)
    clustlist[[i]]<-temp
    
  }
  clustlist	
}



plot_function4<-function(out, cdhit_result, fasta="sRNA"){ #objects from data_preparation function , cdhit_result output from proc_cdhit function
  
  #library(randomcoloR)
  
  #collist<-c("dodgerblue2","darkolivegreen2","darkorange1","darkorchid1","tan3","goldenrod2","lightcyan1","antiquewhite1","mediumpurple1","seagreen2")
  cdl<-unlist(lapply(cdhit_result,length))
  one<-which(cdl==1)
  more<-which(cdl>1)
  nl<-length(more)
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  while(length(more)>length(color)){
	color<-c(color,color)
  }
  if(length(more)>0){
    #collist<- distinctColorPalette(nl)
	collist<- sample(color, nl)
	
    for(i in 1:length(more)){
      names(cdhit_result)[more[i]]<-collist[i]
      
    }
  }
  names(cdhit_result)[one]<-"grey"
  #count<-1
  
  
  numb<-ceiling(length(out)/5)
  nam<-paste(fasta,"_synteny.pdf",sep="")
  pdf(nam,paper = "a4r", width = 0, height = 0)
  
  d_vect<-c()
  s_vect<-c()
  for(i in 1:numb){
    
    
    #for(i in 1:nrow(exist)){
    
    #print(i)
    #	temp<-out[[i]]
    #	s<-mas-as.numeric(temp[1,9])
    #	s_vect<-c(s_vect,s)
    #	temp[,7]<-as.numeric(temp[,7])+s
    #	temp[,8]<-as.numeric(temp[,8])+s
    #	temp[,9]<-as.numeric(temp[,9])+s
    #	temp[,10]<-as.numeric(temp[,10])+s
    #	out[[i]]<-temp
    #}
    count2<-i*5
    count2<-count2-5
    
    mas<-c()
    mi<-c()
    for(ii in 1:5){
      #print(c(i,ii))
      tryCatch({
        te<-(out[[count2+ii]])
        mas<-c(mas,max(as.numeric(te[,9])))
        mi<-c(mi,min(as.numeric(te[,8])))
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
    }
    mas<-max(mas)
    mi<-min(mi)
    s<-c()
    se<-c()
    for(ii in 1:5){
      #print(c(i,ii))
      tryCatch({
        te<-(out[[count2+ii]])
        s2<-mas-as.numeric(te[1,10])[1]
		print(as.numeric(te[1,10])[1])
        s<-c(s,s2+mas)
        se<-c(se,s2+mi)
        te[,8]<-as.numeric(te[,8])+s2
        te[,9]<-as.numeric(te[,9])+s2
        te[,10]<-as.numeric(te[,10])+s2
        te[,11]<-as.numeric(te[,11])+s2
        out[[count2+ii]]<-te
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
    }
    
    s<-max(s)
	null<-which(se==0)
	#if(length(null)>0){
	#	se<-se[-null]
	#}
	
    se<-min(se)
    d_vect<-c(d_vect,s)	
    s_vect<-c(s_vect,se)
	
  }
  
  

  
  for(jj in 1:numb){
    
    d<-d_vect[jj]
    s<-s_vect[jj]
    count<-0
    plot(1,1, type="n",xlim=c(s,d*1),ylim=c(0,5*5),xaxt="n",yaxt="n",xlab="",ylab="")	
    
    for(i in 1:5){
      if((i+(jj-1)*5)<=length(out)){
        temp<-out[[i+(jj-1)*5]]
        
        mi<-min(as.numeric(temp[,8]))
        ma<-max(as.numeric(temp[,9]))
        
        lines(c(mi,ma),c(i+1+count,i+1+count))
        lines(c(mi,ma),c(i+2+count,i+2+count))
        for(j in 1:nrow(temp)){
          
          m<-(as.numeric(temp[j,9])-as.numeric(temp[j,8]))/2+as.numeric(temp[j,8])
          mm<-as.numeric(temp[j,9])-as.numeric(temp[j,8])
          n<--1
          if(temp[j,1]=="+"){
            n<-1
          }
          #if(as.numeric(temp[1,11]==-1)){
          #n<-n*-1
          #}
          if(n==-1){
            n<-0
          }
          #symbols(m,n+count, rectangles=matrix(c(mm,10),1,2), add=TRUE)
          color<-"white"
          
          tcolor<-na.omit(grep(temp[j,5],cdhit_result))
          tcolor2<-na.omit(grep(temp[j,4],cdhit_result))
          tcolor<-c(tcolor,tcolor2)
          if(length(tcolor)>0){
            #if(is.na(temp[j,5])==F){
            color<-names(cdhit_result)[tcolor]
            #}
          }
          rect(as.numeric(temp[j,8]),i+0.5+n+count,as.numeric(temp[j,9]),i+1.5+n+count, col=color)
          text(m,i+n+count+1+0.25,gsub("NA ","",paste(temp[j,4],temp[j,5],sep=" ")),cex=0.4)
		  
          text(m,i+n+count+1-0.25,temp[j,6],cex=0.4)
          
        }
        n<-as.numeric(temp[1,12])
        if(n==-1){
          n<-0
        }
        rect(as.numeric(temp[j,10]),i+0.5+n+count,as.numeric(temp[j,11]),i+1.5+n+count, col=2)
        text(d*1.5/2,i+count,temp[j,13], cex=0.6)
        
        #abline(h=3+count,col=2)
        count<-count+4
      }		
    }	
  }	
  dev.off()
}


read_glassgo<-function(fasta){
  #fasta <- "GLASSgo_output_HG001_02861.fa"
  require(seqinr)
  fast<-read.fasta(fasta)
  x<-read.delim(fasta, header=F, sep="\t")
  nam_pos <- grep(">",x[,1])
  nam <- x[nam_pos,1]
  write.fasta(fast,file.out=fasta, names=names(fast),nbchar=100000 )
  x<-read.delim(fasta, header=F, sep="\t")
  x<-as.character(x[,1])
  nam2_pos <-grep(">",x)
  x[nam2_pos]<-as.character(nam)
  fast<-read.fasta(fasta)
  write.fasta(fast,file.out=paste(fasta, "_old", sep=""), names=names(fast) )
  write.fasta(fast,file.out=fasta, names=seq(1,length(fast)))
  x
}


###############
#    main      #
###############


args <- commandArgs(trailingOnly = TRUE)

input_sRNA <- "GLASSgo_output_HG001_02347__HG001_03272.fa"
output_prefix <- "sRNA_synteny"

for(i in 1:length(args)){
  temp<-strsplit(args[i],"=")
  temp<-temp[[1]]
  temp1<-temp[1]
  temp2<-temp[2]
  assign(as.character(temp1),temp2)
}

require(rentrez)
x<-read_glassgo(input_sRNA)
coor<-export_ncRNA_coordinates(x)
coor <- coor[-is.na(coor['start']), ]
#head(coor)
coor2<-coor
#View(coor2)
coor[,1]<-paste(coor[,1],coor[,3], sep="_")
out<-readFeatures2(coor2, featpath="/mnt/data/db/GLASSgo/data")
coor<-cbind(coor,out[[2]])
colnames(coor)[ncol(coor)]<-"taxid"
out2<-data_preparation(out[[1]],coor2)
tagtable<-locus_tag2org(out2)
get_prot_fasta3(out2[[1]])
cd<-cdhit_run(psi=F,thres=0.4)
cd<-proc_cdhit(cd)
plot_function4(out2[[1]], cd, fasta=output_prefix)

