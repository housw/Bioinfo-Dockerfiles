data_a <- read.table(commandArgs()[5], header=TRUE)
data_r <- read.table(commandArgs()[6], header=TRUE)
data_all<-na.exclude(data_a)
if (nrow(data_all) < 1) stop("No control sequence has a valid CUB index: sequences may be too short")
data_rp<-na.exclude(data_r)
if (nrow(data_rp) < 1) stop("No highly expressed sequence has a valid CUB index: sequences may be too short")
geterrmessage()
OGT<-as.numeric(commandArgs()[7])
res<-1
ref<-1
for(i in 1:100){
sample_all <- sample(data_all$name,nrow(data_all),replace = TRUE)
sample_rp <- sample(data_rp$name,nrow(data_rp),replace = TRUE)
m_all <- data_all$name %in% sample_all
m_rp <- data_rp$name %in% sample_rp
pos_all <- which(m_all=="TRUE")
pos_rp <- which(m_rp=="TRUE")
all_encs_TOT <- data_all$ENCp[pos_all]
all_sis_TOT <- data_all$Si[pos_all]
rp_encs_TOT <- data_rp$ENCp[pos_rp]
rp_sis_TOT <- data_rp$Si[pos_rp]
deltaE_TOT <- ((mean(all_encs_TOT, na.rm=TRUE)-mean(rp_encs_TOT, na.rm=TRUE))/mean(all_encs_TOT, na.rm=TRUE))
deltaS_TOT <- log((mean(rp_sis_TOT, na.rm=TRUE)*(1-mean(all_sis_TOT, na.rm=TRUE))/mean(all_sis_TOT, na.rm=TRUE))/(1-mean(rp_sis_TOT, na.rm=TRUE))) 
F_sara_TOT <- (6.747*deltaE_TOT + 1.184*deltaS_TOT - 1.438)
dbox_TOT <- (1.743-0.0226*OGT-0.7372*F_sara_TOT)
d_TOT <- (-0.1664*dbox_TOT+1)^(-1/0.1664)
res <- append(res,d_TOT)
ref <- append(ref,F_sara_TOT)
}
res2<-paste("Predicted minimum generation time: ",round(mean(res[2:101], na.rm=TRUE),2),"hours +/-",round(sd(res[2:101], na.rm=TRUE),2), seq=" ")
write(res2,file=commandArgs()[8],append=FALSE)



