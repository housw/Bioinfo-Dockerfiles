# author Jens Georg
# selects candidates for copraRNA based on a phylogenetic tree of the input sRNAs
# selects a pre-defined set of organisms close to the organism of interest

#call: R --slave -f  GLASSgo2CopraRNA_balanced_ooi_selection.r --args wildcard=NC_000913,NC_003197 exclude=NC_020260 max_number=15 outfile_prefix=sRNA ooi=NC_000913 sim=3

args <- commandArgs(trailingOnly = TRUE)

ooi<-"NZ_CP018205"
# Staphylococcus aureus, epidermidis, lugdunensis,  warneri,      capitis,    pseudintermedius, saprophyticus
wildcard<-c("NC_007795", "NC_004461", "NC_013893", "NC_020164", "NZ_CP007601", "NC_014925", "NC_007350")
max_number<-20
outfile_prefix<-"sRNA"
exclude<-c()
sim<-3

for(i in 1:length(args)){
	temp<-strsplit(args[i],"=")
	temp<-temp[[1]]
	temp1<-temp[1]
	temp2<-temp[2]
	assign(as.character(temp1),temp2)
 }

wil<-grep(",",wildcard)
if(length(wil)>0){
	wildcard<-strsplit(wildcard,",")[[1]]
} 
 
max_number<-as.numeric(max_number)
sim<-as.numeric(sim)
require(ape)
load("refined_GLASSgo_table.Rdata")


temp<-coor2
if(length(exclude)>0){
	temp_ex<-c()
	for(i in 1:length(exclude)){
		temp_ex1<-grep(exclude[i], coor2[,"fin"])
		if(length(temp_ex1)>0){
			temp_ex<-c(temp_ex,temp_ex1)
		}
	}
	if(length(temp_ex)>0){
	
	temp<-temp[-temp_ex,]
	}
}
coor2<-temp



# clustal omega call with kimura distance matrix for tree generation
clustalo3<-function(coor, positions){
	fasta<-c()
	for(i in 1:length(positions)){
		fasta<-c(fasta, paste(">",coor[positions[i],"fin"],sep=""))
		fasta<-c(fasta, as.character(coor[positions[i],"sequence"]))
	}
	write.table(fasta, file="temp_fasta", row.names=F, col.names=F, quote=F)
	wd<-getwd()
	command<-paste("clustalo -i ", "temp_fasta", " --distmat-out=distmatout.txt --full --output-order=input-order --use-kimura --force --max-hmm-iterations=-1", sep="")
	system(command)
	na<-grep(">", fasta)
	na<-gsub(">","",fasta[na])
	temp<-read.delim("distmatout.txt",sep="",header=F, skip=1)
	unlink("distmatout.txt")
	unlink("temp_fasta")
	temp<-temp[,2:ncol(temp)]
	colnames(temp)<-na
	rownames(temp)<-na
	temp
}

# clustal omega call with percent identity matrix 
clustalo4<-function(coor, positions){
	wd<-getwd()
	fasta<-c()
	for(i in 1:length(positions)){
		fasta<-c(fasta, paste(">",coor[positions[i],"fin"],sep=""))
		fasta<-c(fasta, as.character(coor[positions[i],"sequence"]))
	}
	write.table(fasta, file="temp_fasta", row.names=F, col.names=F, quote=F)
	command<-paste("clustalo -i ", "temp_fasta", " --distmat-out=distmatout.txt --full --percent-id --output-order=input-order --force --max-hmm-iterations=-1", sep="")
	system(command)
	na<-grep(">", fasta)
	na<-gsub(">","",fasta[na])
	temp<-read.delim("distmatout.txt",sep="",header=F, skip=1)
	unlink("distmatout.txt")
	temp<-temp[,2:ncol(temp)]
	colnames(temp)<-na
	rownames(temp)<-na
	temp
}

if(nrow(coor2)<max_number){
	ooi_pos<-grep(ooi, coor2[,"fin"])
	fasta<-c()
	if(length(ooi_pos)>0){
                acc <- strsplit(coor2[ooi_pos,"fin"], "\\.")[[1]][1]
		fasta<-c(paste(">",acc,sep=""), as.character(coor2[ooi_pos,"sequence"]))
		coor2<-coor2[-ooi_pos,]
	}
	for(i in 1:nrow(coor2)){
                acc <- strsplit(coor2[i,"fin"], "\\.")[[1]][1]
                fasta<-c(fasta, paste(">",acc,sep=""))
                fasta<-c(fasta, as.character(coor2[i,"sequence"]))
	}
	nam<-paste(outfile_prefix,"CopraRNA_input_balanced.fasta", sep="_" )
	write.table(fasta, file=nam, row.names=F, col.names=F, quote=F)
}

if(nrow(coor2)>max_number){
	
	
	wildcard<-c(ooi,wildcard)
	pos_wild<-c()
	pos_ooi<-grep(ooi,coor2[,"fin"])[1]
	for(i in 1:length(wildcard)){
		pos_wild<-c(pos_wild,grep(wildcard[i],coor2[,"fin"])[1])
	}
	pos_wild<-unique(na.omit(pos_wild))
	
	
	
	max_number2<-max_number-sim-length(pos_wild)+1
	pos<-seq(1,nrow(coor2))
	if(length(pos_wild)>0){
		pos<-pos[-pos_wild]
	}
	pos<-unique(c(pos_ooi,pos))
	dis<-clustalo3(coor2, pos)
	dis2<-clustalo4(coor2, pos)
	dis<-as.dist(dis)
	clus<-(hclust(dis,method="average"))
	plot(clus)
	knum<-min(max_number2,length(clus$labels)-1)
	if(knum<2){
		knum<-2
	}
	clus2<-rect.hclust(clus,k=knum)
	
	out<-c()
	for(i in 1:length(clus2)){
		temp<-clus2[[i]]
		temp_ooi<-grep(ooi,names(temp))
		if(length(temp_ooi)==0){
			temp2<-sample(length(temp),1)
			out<-c(out, names(temp)[temp2])
		}
	}
	
	
	
	ooil<-sort(dis2[grep(ooi, colnames(dis2)),], decreasing=T)
	
	sel_ooi<-na.omit(match(out,names(ooil)))
	if(length(sel_ooi)>0){
		ooil<-ooil[-sel_ooi]
	}
	ident<-which(ooil==100)
	if(length(ident)>0){
		ooil<-ooil[-ident]
	}
	sel<-c()

	close_orgs<-c()
	iii<-0
	while(length(close_orgs)<sim){
	knum2<-min(max_number2,length(clus$labels)-1)-iii
	if(knum2<2){
	break
	}
	clus3<-rect.hclust(clus,k=knum2)
	ooi1<-grep(ooi, names(unlist(clus3)))
	len<-as.numeric(summary(clus3)[,1])
	su<-0
	ii<-1
	while(su<ooi1){
		su<-su+len[ii]
		ii<-ii+1
	}
	ii<-ii-1
	close_orgs<-sort(ooil[intersect(names(clus3[[ii]]),names(ooil))],decreasing=T)

		
	iii<-iii+1
	}


	if(length(close_orgs)>sim){
		n<-length(close_orgs)%/%sim
		n<-seq(1,n*sim,by=n)
		n<-names(close_orgs)[n]
		sel<-unique(c(sel,n))
	}


	if(length(close_orgs)<=sim){
	sel<-unique(c(sel,names(close_orgs)))

	}
	out_old<-out
	out<-c(coor2[pos_wild,"fin"],sel,out)
	out<-match(out,coor2[,"fin"])
	fasta<-c()
	for(i in 1:length(out)){
            acc <- strsplit(coor2[out[i],"fin"], "\\.")[[1]][1]
            fasta<-c(fasta, paste(">",acc,sep=""))
            fasta<-c(fasta, as.character(coor2[out[i],"sequence"]))
	}
	nam<-paste(outfile_prefix,"CopraRNA_input_balanced_ooi_neighbourhood.fasta", sep="_" )
	fasta<-gsub("\\..*","",fasta)
	write.table(fasta, file=nam, row.names=F, col.names=F, quote=F)
	dis<-clustalo3(coor2, seq(1,nrow(coor2)))
	dis<-as.dist(dis)
	clus<-(hclust(dis,method="average"))
	clus<-as.phylo(clus)
	lab<-clus$tip.label

	nam_selected<-match(out_old,lab)
	nam_wildcard<-match(coor2[pos_wild,"fin"],lab)
	nam_neighbourhood<-match(sel,lab)
	nam_ooi<-grep(ooi,lab)

	lab<-match(lab,coor2[,"fin"])
	lab<-coor2[lab,"nam2"]
	clus$tip.label<-lab
	nam<-paste(outfile_prefix,"tree_coprarna_candidates_balanced_ooi_neighbourhood.pdf", sep="_" )
	pdf(nam)
	colo<-rep("1",length(lab))
	colo[nam_selected]<-"dodgerblue1"
	colo[nam_wildcard]<-"olivedrab2"
	colo[nam_neighbourhood]<-"orangered"
	colo[nam_ooi]<-"purple1"
	par(mar=c(3, 1, 1, 1), xpd=TRUE)
	plot(clus,tip.color=colo, cex=0.5 )

	legend("bottom",  inset=c(-0.05),bty="n", legend=c("organism of interst (ooi)","pre-selected organisms","selected organisms","close to ooi"), text.col=c("purple1","olivedrab2","dodgerblue1","orangered"),cex=0.6)
	par(xpd=FALSE)
	dev.off()
	
}

unlink("Rplots.pdf")

