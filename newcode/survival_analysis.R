##**画出各个肿瘤的生存曲线
##**TCGA match samples and gene related Harzard_coef
##**summary gene_Harzard_coef
##**check gene_Harzard_coef in different cancers



#########画出各个肿瘤的生存曲线
library(survival)
cancer_type=sort(c('BRCA','COAD','LUAD','BLCA','LIHC','LUSC','PRAD','ESCA','STAD','HNSC','KICH','THCA','KIRC','KIRP','PAAD'))
#cancer_type=c('BRCA','COAD','LUAD','BLCA','LIHC','LUSC',
pdf("/Users/sunhy/Desktop/survival curve.pdf")
par(mfrow=c(2,3))
for (ci in 1:length(cancer_type))
{
  cancer=cancer_type[ci]
  data <- read.delim(paste("/Users/sunhy/Documents/PA/survival_analysis/TCGA_file/",cancer,"_follow_up_status.txt",sep=""))
  tsurv<- Surv(data$follow_up, data$vital_status)
  plot(survfit(tsurv ~ 1),xlab="time-day",
        ylab="Estimated Survival rate",
        main=paste(cancer,"_survival curve"))
}
dev.off()

#####
##### TCGA match samples and gene related Harzard_coef
library(OIsurv)
source("/Users/sunhy/Documents/deskfile2/Fenton_basic_functions.R")
cancer_type=c('BRCA','COAD','LUAD','BLCA','LIHC','LUSC','PRAD','ESCA','STAD','HNSC','KICH','THCA','KIRC','KIRP')
genelist=read.table("/Users/sunhy/Documents/deskfile/sep_RNASEQ/gene_symbol_list_exceptSLC35E2")
g1=c("SLC2A1","SLC2A4","HK1","HK2","GCK","G6PD","PRPS1","PRPS2","GPI","PFKM","PFKP","PFKL","ALDOA","ALDOB","ALDOC","TKT","TPI1",
            "GPD1","GPD2","GPD1L","GAPDH","PGK1","PGK2","PHGDH","SPTLC1","SPTLC2","SPTLC3","PGAM1","PGAM2","PGAM4","PGAM5","ENO1","ENO2","ENO3","PKM2","PKM","PKLR",
            "PCK1","PCK2","FBP1","FBP2","G6PC",
            "ME1","ME2","ME3","MDH1","MDH2","BRP44L","BRP44","ACLY","PDHA1","PDHA2","PDHB","DLAT","CS","ACO1","ACO2","SLC25A1","IDH1","IDH2","IDH3A","IDH3B","IDH3C","OGDH","DLST","DLD",
            "SUCLG1","SUCLG2","SDHA","SDHB","SDHC","SDHD","FH","SLC25A11","SLC1A5","GLS","GLS2","GLUD1","GLUD2","GOT1","GOT2","GPT","GPT2","PC","PSAT1","PDXP","PSPH",
            "ACACA","ACACB","FASN","HMGCR","ACAT1","ACAT2","ACAA1","ACAA2",
            "LDHA","LDHB","SLC1A1","SLC1A2","SLC1A3","SLC1A6","SLC1A7","SLC17A6","SLC17A7","SLC17A8")
mid1=match(intersect(g1,genelist$V1),genelist$V1)

coef_list=list()
for (ci in 1:length(cancer_type))
{
  cancer=cancer_type[ci]
  surv_data <- read.delim(paste("/Users/sunhy/Desktop/survival_analysis/TCGA_file/",cancer,"_follow_up_status.txt",sep=""))
  load(paste("/Users/sunhy/Documents/deskfile/sep_RNASEQ/TCGA_RData/",cancer,"_tumor_RNASeqV2.genes.normalized.RData",sep=""))
  t_name=colnames(tdata)
  tdata[which(is.na(tdata))]=0  
  log_t=log(tdata+1)
  
  newname=rewriteTCGAname2(t_name,st="\\-")
  C_mid1=match(intersect(surv_data$TCGA_barcode,newname),newname)
  C_mid2=match(intersect(surv_data$TCGA_barcode,newname),surv_data$TCGA_barcode)
  
  LLL=cbind(as.vector(surv_data$TCGA_barcode[C_mid2]),as.vector(newname[C_mid1]))
  
  tsurv<- Surv(surv_data$follow_up[C_mid2], surv_data$vital_status[C_mid2])
  #coxph.fit <- coxph(tsurv ~ t(log_t[mid1,C_mid1]), method="breslow")
  #coxph.fit <- coxph(tsurv ~ t(log_t[mid1,C_mid1]))
  coef=c()
  for (mi in 1:length(mid1))
  {
    coxph.fit <- coxph(tsurv ~ log_t[mid1[mi],C_mid1])
    lll=summary(coxph.fit)
    t_coef=c(lll$coefficients[1],lll$coefficients[5])
    coef=rbind(coef,as.vector(t_coef))
  } 
  row.names(coef)=genelist$V1[mid1]
  coef_list[[ci]]=coef
  names(coef_list)[ci]=cancer
}
save(coef_list,file="/Users/sunhy/Desktop/survival_analysis/result/Harzard_coef_list.RData")

########
########summary gene_Harzard_coef
load("/Users/sunhy/Desktop/survival_analysis/result/Harzard_coef_list.RData")
uplist=c()
downlist=c()
for (ci in 1:length(coef_list))
{
  tcoef=coef_list[[ci]]
  tid1=intersect(which(tcoef[,2]<0.05),which(tcoef[,1]>0))
  tg_up=row.names(tcoef)[tid1]
  tid2=intersect(which(tcoef[,2]<0.05),which(tcoef[,1]<0))
  tg_down=row.names(tcoef)[tid2]
  uplist=c(uplist,as.vector(tg_up))
  downlist=c(downlist,as.vector(tg_down))
}
suplist=sort(table(uplist),decreasing=T)
sdownlist=sort(table(downlist),decreasing=T)

#########check gene_Harzard_coef in different cancers
PGK1_cancers=c()
for (coi in 1:length(coef_list))
{
  tco=coef_list[[coi]][22,]
  PGK1_cancers=rbind(PGK1_cancers,as.vector(tco))
}
row.names(PGK1_cancers)=names(coef_list)



###########
###########
#install.packages("survival")
library(survival)
aml=aml[order(aml$time),]
with(aml, plot(time, type="h"))
aml.survfit = survfit(Surv(time, status == 1) ~ 1, data=aml)

tsurv=Surv(aml$time, aml$status)
plot(survfit(tsurv ~ 1),xlab="time",
     ylab="Estimated Survival Function",
     main=paste("survival curve"))  ###conf.int=FALSE 没有虚线 confidence bounds.
plot(aml.survfit, xlab = "Time (weeks)", ylab="Proportion surviving", conf.int=FALSE, main="Survival in AML")

summary(aml.survfit)


surv.by.aml.rx = survfit(Surv(time, status == 1) ~ x, data = aml)
summary(surv.by.aml.rx)
plot(surv.by.aml.rx, xlab = "Time", ylab="Survival",col=c("black", "red"), lty = 1:2, main="Kaplan-Meier Survival vs. Maintenance in AML")
legend(100, .6, c("Maintained", "Not maintained"), lty = 1:2, col=c("black", "red"))


surv.diff.aml= survdiff(Surv(time, status == 1) ~ x, data=aml)
surv.diff.aml

plot(survfit(Surv(time, status == 1) ~ x, data=aml))

#install.packages("ISwR")
library(ISwR)
help(melanom) # description of the melanoma data

surv.diff.sex = survdiff(Surv(days, status == 1) ~ sex, data = melanom)
surv.diff.sex

coxph.sex = coxph(Surv(days, status == 1) ~ sex, data = melanom)
coxph.sex.thick = coxph(Surv(days, status == 1) ~ sex + log(thick), data = melanom)

cox.zph(coxph.sex)



#install.packages("rpart")
library(rpart)
head(stagec)
fit <- rpart(Surv(pgtime, pgstat) ~ age + eet + g2 + grade + gleason + ploidy, data=stagec)
plot(fit, uniform=T, branch=.4, compress=T)
text(fit, use.n=T)
print(fit)


fit <- rpart(Surv(survtime, status == 1) ~ lll[,1] + lll[,2] + lll[,3])

