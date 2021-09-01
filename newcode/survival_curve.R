library(survival)
library(survminer)
allset_path=paste0("D:/WORK/样本多分类/BRCAcompare/BRCAall.csv")
allset=read.csv(allset_path,row.names = 1)
sample_lable=allset$BRCA_lable
t_lable_index=which(sample_lable!="1_ndata")
t_Data=allset[t_lable_index,]
t_lable=sample_lable[t_lable_index]
sample_name=rownames(t_Data) #长名字

sur_file_path=paste0("D:/WORK/样本多分类/survival_ana/BRCA_follow_up_status.txt")
sur=read.table(sur_file_path,header = T)
samples=sur$TCGA_barcode #短名字

match_index=c()
remain_index=c()
for(i in 1:length(samples)){
  t_index=grep(samples[i],sample_name)
  print(paste("i:",i,"t_index:",t_index))
  if(length(t_index)==1){
    remain_index=c(remain_index,i)
    match_index=c(match_index,t_index)
  }
}

remain_sur=sur[remain_index,]
subtype_t=t_lable[match_index]
sur_information=cbind(remain_sur,subtype_t)
SurvObject<-with(sur_information,Surv(follow_up,vital_status))
fit=survfit(Surv(follow_up,vital_status)~subtype_t,data=sur_information)
surv_summary(fit)
png(file = "D:/WORK/样本多分类/survival_ana/BRCA_surcurve.png",
    width=700,height=700)
ggsurvplot(fit,pval=T,conf.int=T,risk.table = T,
           risk.table.col="strata",linetype = "strata",
           surv.median.line = "hv",ggtheme = theme_bw(),
           palette = c("red","blue","green"))
dev.off()
survdiff(Surv(follow_up,vital_status)~subtype_t,data=sur_information)
