

#call: R --slave -f  promoter_sequence_fetcing_from_coor.r --args datapath=refined_GLASSgo_table.Rdata output_prefix="sRNA_promoter" 


args <- commandArgs(trailingOnly = TRUE)

datapath <- "refined_GLASSgo_table.Rdata"
output_prefix <- "sRNA_promoter" 

for(i in 1:length(args)){
  temp<-strsplit(args[i],"=")
  temp<-temp[[1]]
  temp1<-temp[1]
  temp2<-temp[2]
  assign(as.character(temp1),temp2)
}

load(datapath)

require(rentrez)


get_sequence<-function(coor, up=0, down=0, name="sequences", from_start=FALSE){
  out<-c()
  for(i in 1:nrow(coor)){
    if(from_start==TRUE){
      s<-as.numeric(coor[i,3])-up
      e<-as.numeric(coor[i,3])+down
    }
    if(from_start==FALSE){
      s<-as.numeric(coor[i,3])-up
      e<-as.numeric(coor[i,4])+down
    }
    st<-1
    if(coor[i,2]=="-"){
      st<-2
      if(from_start==TRUE){
        s<-as.numeric(coor[i,4])-up
        e<-as.numeric(coor[i,4])+down
      }
     if(from_start==FALSE){
        s<-as.numeric(coor[i,4])-down
        e<-as.numeric(coor[i,3])+up
      }
    }
    prom<-entrez_fetch(db="nucleotide", id=coor[i,"fin"], rettype="fasta", retmode="text", seq_start=s, seq_stop=e, strand=st)
    prom<-unlist(strsplit(prom, "\n"))
    print(paste(i,"/",nrow(coor), "_", coor[i,"fin"], sep=""))
    name1<-paste(">",coor[i,"fin"], "_",coor[i,6], sep="")
    name1<-gsub(" ","_",name1)
    prom[1]<-name1
    out<-c(out,prom)
  }
  write.table(out, file=paste(name,".fasta",sep=""), sep="\t", row.names=F,col.names=F,quote=F)
  out
}

#coor2<-coor[na.omit(inp2),]
#coor2[,1]<-gsub("_.*","",coor2[,1])
prom<-get_sequence(coor2,up=100,down=50,name=output_prefix)

