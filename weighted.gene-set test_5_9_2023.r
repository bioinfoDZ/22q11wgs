
if(!require("car")){
  install.packages("car")
}
library(car)

if(!require("qqman")){
  install.packages("qqman")
}
library(qqman)

if(!require("CompQuadForm")){
  install.packages("CompQuadForm")
}
library(CompQuadForm)


adapt.null.sim<-function(m,num.rep=50000){ 
  nulldistr<-numeric() 
  nullW<-rchisq(num.rep,df=m-1)
  nullC<-rchisq(num.rep,df=1)
  min.pweighted <- numeric()
  for (r in 1:num.rep){
    null.weighted <- numeric()
    for (gamma in seq(0,1,0.05)){
      null.weighted <- c(null.weighted, davies(nullC[r]+(1-gamma)*nullW[r],lambda=c(1,1-gamma), h=c(1,m-1))$Qq)
    }
    min.pweighted <-c(min.pweighted,min(null.weighted))
  }
  return(min.pweighted)
}

adapt.test <- function(z,weight,m){
#define T0
  T0<-(weight%*%z)^2/(t(weight)%*%weight)  
  pt0<-pchisq(T0,df=1,lower.tail=FALSE)
  
  ## T1
  Q=z%*%z 
  pt1<-pchisq(Q,df=length(z),lower.tail=FALSE)
  
  
  padpted <- numeric()
  for (gamma in seq(0,1,0.05)){
    padpted <- c(padpted,davies((1-gamma)*T0+gamma*Q,lambda=c(1,gamma), h=c(1,length(z)-1))$Qq)
  }
  
  
  min.p <-min(padpted)
  loc <- seq(0,1,0.05)[which.min(padpted)]
  
  null = adapt.null.sim(m=m, num.rep=50000)
  padp.sim<-sum(null<min.p)/length(null) 
  
  return(list(format(padp.sim,digits=10),
              sum(null<min.p),
              length(null) ,
              format(pt0,digits = 6), 
              format(pt1,digits = 6), 
              loc))
}



setwd("/Users/yingjiezhao/Desktop/Rdata")
#readin data,CTD_all.6_cadd10_loP1_ccrs80.fet.no.filter.subset
headers = read.table("CTD_all.6_cadd10_phast.5_ccrs80.fet.no.filter.txt", header = F, nrows = 1, as.is = T)
gene_out = read.table("CTD_all.6_cadd10_phast.5_ccrs80.fet.no.filter.txt",skip = 1, header = F,na.strings=c("NA",".", "None"))
colnames(gene_out)= headers

for ( i in c(51:51)) { ##39:130 total gene sets##
  gene_set_name<-headers[1,i]
  ### extract one gene-set
  gene_out_set = gene_out[gene_out[,i]==1 ,]
  
  ### get z-score  
  sign = ((gene_out_set$odds_ratio1>=1)*1-0.5)*2
  sign1=((gene_out_set$odds_ratio1>=1)-0.5)*2
  z = sign*qnorm(gene_out_set$p_value1/2,lower.tail=FALSE)
  zlen<-length(z)
  
  ### weight by gene expression 
  weight =  (gene_out_set$mouse_mean_WT_OFT_TPM_2+gene_out_set$mouse_mean_WT_LV_TPM_3+gene_out_set$mouse_mean_WT_PA_TPM_3)/3
  weight[is.na(weight)] = 1
  weight = log(weight)
  
  
  ### weight by gene connectivity
  #weight =  gene_out_set$Interactions_BioGRID
  #weight =  gene_out_set$Interactions_ConsensusPathDB
  #weight =  gene_out_set$Interactions_IntAct
  #weight =  (gene_out_set$Interactions_IntAct+gene_out_set$Interactions_ConsensusPathDB+gene_out_set$Interactions_BioGRID)/3
  #weight[is.na(weight)] = 1
  #weight =  log(weight)
  
  ### weight by gene constraint scores
  #weight =  gene_out_set$virlof_percentile
  #weight[is.na(weight)] = 100
  #weight=100-weight
  
   tryCatch({
    ret1<-adapt.test (z,
                      weight,
                      dim(gene_out_set)[1])
    writeLines(paste(c(gene_set_name, zlen, paste(ret1,collapse = "\t")), collapse = "\t"))
  }, error=function(e){
    writeLines(paste(c(gene_set_name, zlen, "NA\tNA\tNA\tNA\tNA\tNA"), collapse = "\t"))
  })
  
} 








