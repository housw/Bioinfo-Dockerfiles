data_a <- read.table(commandArgs()[5], header=TRUE)
data_r <- read.table(commandArgs()[6], header=TRUE)
data_all<-na.exclude(data_a)
if (nrow(data_all) < 1) stop("No control sequence has a valid CUB index: sequences may be too short")
data_rp<-na.exclude(data_r)
if (nrow(data_rp) < 1) stop("No highly expressed sequence has a valid CUB index: sequences may be too short")
geterrmessage()
OGT<-as.numeric(commandArgs()[7])
all_encs <- data_all$ENCp
all_sis <- data_all$Si
rp_encs <- data_rp$ENCp
rp_sis <- data_rp$Si
deltaE_lst <- ((all_encs - rp_encs)/all_encs)
deltaS_lst <- log((rp_sis*(1-all_sis)/all_sis)/(1-rp_sis)) 
F_sara_lst <- (6.747*deltaE_lst + 1.184*deltaS_lst - 1.438)
dbox_lst <- (1.743-0.0226*OGT-0.7372*F_sara_lst)
d_lst <- (-0.1664*dbox_lst+1)^(-1/0.1664)
res<-paste("Predicted minimum generation time: ",round(mean(d_lst, na.rm=TRUE),2),"hours +/-",round(sd(d_lst, na.rm=TRUE),2), seq=" ")
write(res,file=commandArgs()[8],append=FALSE)



