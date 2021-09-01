################# cancer plot suvr curve IGG
rm(list = ls())
#library(OIsurv)
library(survival)
library(survminer)
source("/Users/sunhy/Documents/PA/Fenton_basic_functions.R")
genelist=as.vector(read.table("/Users/sunhy/Desktop/hy_data/gene_symbol_list_exceptSLC35E2")$V1)

data_path=c("/Users/sunhy/Documents/all_project/gbm_survival/fenton and survival/TCGA_RData/")
cancer_type=c('BRCA','COAD','LUAD','BLCA','LIHC','LUSC','PRAD','ESCA','STAD','HNSC','KICH','THCA','KIRC','KIRP')

#fe=c("TFR2","TFRC","HMOX1","SLC25A37","STEAP3","HAMP","FTH1","FTL","SLC40A1")
#fe=c("TFR2","TFRC","SOD1","SOD2","NOX3","SLC40A1")
load("/Users/sunhy/Documents/PA/genesets.RData")

#candi=as.vector(genesets[[which(names(genesets)=="ER_STRESS")]]$V1); weight=rep(1,12)
#candi=as.vector(genesets[[which(names(genesets)=="DNA_integrity_maintain")]]$V1); 
#candi=c("TFR2","TFRC","STEAP3","SLC40A1"); 
candi=c("NOX3","TFRC","TXNRD1","SOD2")

mid1=match(intersect(candi,genelist),genelist);          weight=c(rep(1,(length(mid1)-1)),-1) ###################
pdf(file="/Users/sunhy/Desktop/FR_surv_iGG.pdf")
survdiff.list=list()

for (ci in 1:length(cancer_type))
{
  cancer=cancer_type[ci]
  load(paste(data_path,cancer,"_tumor_RNASeqV2.genes.normalized.RData",sep=""))
  tdata[which(is.na(tdata))]=0  
  surv_data <- read.delim(paste("/Users/sunhy/Documents/PA/survival_analysis/TCGA_file/",cancer,"_follow_up_status.txt",sep=""))
  
  load(paste(data_path,cancer,"_tumor_RNASeqV2.genes.normalized.RData",sep=""))
  t_name=colnames(tdata)
  tdata[which(is.na(tdata))]=0  
  newname=rewriteTCGAname2(t_name,st="\\-")
  C_mid1=match(intersect(surv_data$TCGA_barcode,newname),newname)
  C_mid2=match(intersect(surv_data$TCGA_barcode,newname),surv_data$TCGA_barcode)
  tsubdata1=log(tdata[mid1,C_mid1]+1)
  tu=apply(tsubdata1,1,mean)
  tv=apply(tsubdata1,1,var)
  
  subdata1=(tsubdata1-tu)*weight/tv
  
  protea_infor=colSums(subdata1)
  
  ex_level=rep("high-expression",length(protea_infor))
  ex_level[which(protea_infor<0)]=c("low-expression")
  t_surv_infor=cbind(surv_data$follow_up[C_mid2], surv_data$vital_status[C_mid2],ex_level)
  
  tsurv<- Surv(surv_data$follow_up[C_mid2], surv_data$vital_status[C_mid2])  
  
  surv.by.aml.rx = survfit(tsurv ~ t_surv_infor[,3])
  surv.diff = survdiff(tsurv ~ t_surv_infor[,3])
  
  survdiff.list[[ci]]=surv.diff
  names(survdiff.list)[ci]=cancer
  
  LALA=summary(surv.by.aml.rx)
  plot(surv.by.aml.rx, xlab = "Time", ylab="Survival",col=c("black", "red"), lty = 1:2, main=paste("Kaplan-Meier Survival of ER ",cancer,sep=""))
  legend(100, .6, c("high-expression", "low-expression"), lty = 1:2, col=c("black", "red"))  
}
dev.off()

################################
############# cancer surv difference 
################################
rm(list = ls())
library(OIsurv)
source("/Users/sunhy/Documents/PA/Fenton_basic_functions.R")
genelist=as.vector(read.table("/Users/sunhy/Documents/essential code/gene_symbol_list_exceptSLC35E2")$V1)
load("/Users/sunhy/Desktop/fenton and survival/FR_list.RData")
#Proteasome=FR_list[[which(names(FR_list)=="Proteasome")]]
#Proteasome=c("KIF4B","KIF4A","PGK1","AURKA","KPNA2","BIRC5","PLK1")
#Proteasome=c("PGK1","LDHA")
#Proteasome=c("PGK1","AURKA","KPNA2","BIRC5","PLK1")
#Proteasome=c("NOX3","TFRC","TXNRD1","SOD3")
Proteasome=c("COL1A2","TRAM2")

data_path=c("/Users/sunhy/Desktop/fenton and survival/TCGA_RData/")
cancer_type=c('BRCA','COAD','LUAD','BLCA','LIHC','LUSC','PRAD','ESCA','STAD','HNSC','KICH','THCA','KIRC','KIRP')
load("/Users/sunhy/Documents/PA/genesets.RData")
load("/Users/sunhy/Desktop/fenton and survival/TCGA_up_down_gene_list.RData")

#pdf(file="/Users/sunhy/Desktop/DNA_integrity_surv_iGG.pdf")
survdiff.list=list()
pvalues=c()
for (ci in 1:length(cancer_type))
{
  cancer=cancer_type[ci]
  load(paste(data_path,cancer,"_tumor_RNASeqV2.genes.normalized.RData",sep=""))
  tdata[which(is.na(tdata))]=0  
  cid=which(names(up_down_gene_list)==paste(cancer,"_upgene",sep=""))
  upgene=up_down_gene_list[[cid]]#####修改index 
  
  mid1=match(intersect(Proteasome,genelist),genelist)  ##########whether choosing up regulated genes
  
  surv_data <- read.delim(paste("/Users/sunhy/Desktop/fenton and survival/survival_analysis/TCGA_file/",
                                cancer,"_follow_up_status.txt",sep=""))
  load(paste(data_path,cancer,"_tumor_RNASeqV2.genes.normalized.RData",sep=""))
  t_name=colnames(tdata)
  tdata[which(is.na(tdata))]=0  
  newname=rewriteTCGAname2(t_name,st="\\-")
  C_mid1=match(intersect(surv_data$TCGA_barcode,newname),newname)
  C_mid2=match(intersect(surv_data$TCGA_barcode,newname),surv_data$TCGA_barcode)
  tsubdata1=log(tdata[mid1,C_mid1]+1)
  tu=apply(tsubdata1,1,mean)
  tv=apply(tsubdata1,1,var)
  
  #weight=c(rep(1,(length(mid1)-1)),-1)##########
  weight=rep(1,length(mid1))
  
  subdata1=(tsubdata1-tu)*weight/tv
  
  protea_infor=colSums(subdata1)
  
  ex_level=rep("high-expression",length(protea_infor))
  ex_level[which(protea_infor<0)]=c("low-expression")
  t_surv_infor=cbind(surv_data$follow_up[C_mid2], surv_data$vital_status[C_mid2],ex_level)
  
  tsurv<- Surv(surv_data$follow_up[C_mid2], surv_data$vital_status[C_mid2])  
  
  surv.diff = survdiff(tsurv ~ t_surv_infor[,3])
  pvalue=1-pchisq(surv.diff$chisq,length(surv.diff$n)-1)
  pvalues=c(pvalues,pvalue)
}