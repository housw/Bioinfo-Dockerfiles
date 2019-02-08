# author Jens Georg
# removes 100% identical sequences and keeps pre-selected organsisms based on their Refseq ID
# 

#call: R --slave -f  GLASSgo_2_copraRNA_exclusion_script.r --args datapath=full_GLASSgo_table.Rdata 


args <- commandArgs(trailingOnly = TRUE)

datapath<-"full_GLASSgo_table.Rdata"
# Staphylococcus aureus, epidermidis, lugdunensis,  warneri,      capitis,    pseudintermedius, saprophyticus
wildcard<-c("NC_007795","NC_004461", "NC_013893", "NC_020164", "NZ_CP007601", "NC_014925", "NC_007350")
ooi<-"NZ_CP018205"

for(i in 1:length(args)){
  temp<-strsplit(args[i],"=")
  temp<-temp[[1]]
  temp1<-temp[1]
  temp2<-temp[2]
  assign(as.character(temp1),temp2)
}

load(datapath)

wildcard<-unique(c(ooi,wildcard))
pre<-c()
for(i in 1:length(wildcard)){
  pre_temp<-grep(wildcard[i],coor[,"fin"] )
  if(length(pre_temp)>0){
    pre<-c(pre,pre_temp)
  }
}

cop_pre<-c()
prestring<-c()
if(length(pre)>0){
  cop_pre<-coor[pre,"fin"]
  prestring<-coor[pre,]
  #coor<-coor[-pre,]
}

if(length(prestring)>0){	
  if(is.matrix(prestring)==F){
    prestring<-t(as.matrix(prestring))
  }
}

fasta3<-c()
for(i in 1:nrow(coor)){
  fasta3<-c(fasta3, paste(">",coor[i,"fin"],sep=""))
  fasta3<-c(fasta3, as.character(coor[i,"sequence"]))
}


write.table(fasta3, file="temp_fasta", row.names=F, col.names=F, quote=F)

# calculate percent identity matrix

clustalo2<-function(coor){
  wd<-getwd()
  
  command<-paste("clustalo -i ", "temp_fasta", " --distmat-out=distmatout.txt --full --percent-id --output-order=input-order --force --max-hmm-iterations=-1", sep="")
  system(command)
  na<-c(coor[,1])
  
  #na<-c("query",coor[,1])
  temp<-read.delim("distmatout.txt",sep="",header=F, , skip=1)
  unlink("distmatout.txt")
  temp<-temp[,2:ncol(temp)]
  colnames(temp)<-na
  rownames(temp)<-na
  temp
}
dis<-clustalo2(coor)
unlink("temp_fasta")
temp<-na.omit(match(colnames(dis)[2:ncol(dis)],coor[,1]))


### exclude 100% identicals to wildcard organisms
if(length(prestring)>0){
  for(i in 1:length(pre)){
    g<-match(prestring[,1],colnames(dis))
    temp<-prestring[i,1]
    temp<-grep(temp, colnames(dis))
    temp1<-which(dis[temp,]==100)
    ex<-na.omit(match(g,temp1))
    if(length(ex)>0){
      temp1<-temp1[-ex]
    }
    
    
    if(length(temp1)>0){
      dis<-dis[-temp1,-temp1]
      
    }
  }
}

### exclude remaining 100% identicals 
i<-1
while(nrow(dis)>i){
  
  g<-c()
  if(length(prestring)>0){
    g<-match(prestring[,1],colnames(dis))
  }
  temp<-which(dis[i,]==100)
  ex<-na.omit(match(c(i,g), temp))
  if(length(ex)>0){
    temp<-temp[-ex]
  }
  
  if(length(temp)>0){
    if(i==1){
      
      if(colnames(dis)[1]=="query"){
        temp<-sort(temp)
        if(length(temp)>1){
          temp<-temp[2:length(temp)]
        }
        if(length(temp)==1){
          temp<-c()
        }
      }
    }
    
  }
  if(length(temp)>0){
    dis<-dis[-temp,-temp]
    
  }
  i<-i+1
}


temp<-na.omit(match(colnames(dis),coor[,1]))
temp_cop_ref<-na.omit(match(cop_pre, coor[temp,"fin"]))
cop_pre_w_ooi<-c()
if(length(cop_pre)>1){
  cop_pre_w_ooi<-cop_pre[2:length(cop_pre)]
}
temp_cop_ref_w_ooi<-na.omit(match(cop_pre_w_ooi, coor[temp,"fin"]))




temp2<-temp
if(length(temp_cop_ref)>0){
  temp2<-temp[-temp_cop_ref]
}
temp_5<-temp
if(length(temp_cop_ref_w_ooi)>0){
  temp_5<-temp[-temp_cop_ref_w_ooi]
}

fasta_wo_wild_w_ooi<-c()
for(i in 1:length(temp_5)){
  fasta_wo_wild_w_ooi<-c(fasta_wo_wild_w_ooi, paste(">",coor[temp_5[i],"fin"],sep=""))
  fasta_wo_wild_w_ooi<-c(fasta_wo_wild_w_ooi, as.character(coor[temp_5[i],"sequence"]))
  
}



fasta3<-c()
for(i in 1:length(temp)){
  fasta3<-c(fasta3, paste(">",coor[temp[i],"fin"],sep=""))
  fasta3<-c(fasta3, as.character(coor[temp[i],"sequence"]))
  
}

fasta_wo_wild<-c()
for(i in 1:length(temp2)){
  fasta_wo_wild<-c(fasta_wo_wild, paste(">",coor[temp2[i],"fin"],sep=""))
  fasta_wo_wild<-c(fasta_wo_wild, as.character(coor[temp2[i],"sequence"]))
  
}
write.table(fasta_wo_wild_w_ooi, file="temp_fasta_wo_wild_w_ooi.fasta", row.names=F, col.names=F, quote=F)
write.table(fasta_wo_wild, file="temp_fasta_wo_wild.fasta", row.names=F, col.names=F, quote=F)
write.table(fasta3, file="temp_fasta", row.names=F, col.names=F, quote=F)

fasta_unique<-fasta3 # unique CopraRNA "positive" sequences 
write.table(fasta_unique, file="un_fasta", row.names=F, col.names=F, quote=F)


coor2<-coor[temp,]
system("rm -f refined_GLASSgo_table.Rdata")
save(coor2, file="refined_GLASSgo_table.Rdata")

