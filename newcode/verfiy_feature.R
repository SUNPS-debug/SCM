library(Rcpp)
library(RSNNS)
pick_key_gene<-function(dif_gene_set,spe_edges_loc){
   features_gene=list()
   for(i in 1:ntype){
     unique_index=unique(as.numeric(spe_edges_loc[[i]])) 
     features_gene[[i]]=dif_gene_set[[i]][unique_index]
   }
   return(features_gene)
}
ntype=4
dif_gene_set=all.result$dif$dif_gene_set
spe_edges_loc=all.result$spe_edges
features_gene_list=pick_key_gene(dif_gene_set,spe_edges_loc)
features_genes_fre=table(unlist(features_gene_list))
com_genes=names(features_genes_fre[which(features_genes_fre==2)])
#全部的feature_gene
features_gene=Reduce(union,features_gene_list)
#去掉重叠后的feature_gene
features_gene2=setdiff(features_gene,com_genes)
feature_select=read.csv("D:/WORK/���������/feature_list.csv",header=0)
inter=intersect(feature_select,features_gene)
generate_NN_data<-function(new_cleandata,features_gene){
   class_data=list()
   ran_data=list()
   gene_num=dim(new_cleandata[[1]])[2]
   ran_index=sample(gene_num,size=length(features_gene))
   for(i in 1:ntype){
      sam_num=dim(new_cleandata[[i]])[1]
      class_data[[i]]=new_cleandata[[i]][,features_gene]
      ran_data[[i]]=new_cleandata[[i]][,ran_index]
      target=rep(cancertype[i],sam_num)
      class_data[[i]]=cbind(class_data[[i]],target)
      ran_data[[i]]=cbind(ran_data[[i]],target)
   }
   cla=do.call(rbind,class_data)
   ran=do.call(rbind,ran_data)
   df=list(cla=cla,ran=ran)
   return(df)
}

df=generate_NN_data(new_cleandata,features_gene2)

MPL_pre<-function(data){
   suf_df <- data[sample(nrow(data)),]
   dfValues <- suf_df[,-ncol(suf_df)]
   dfTargets <- decodeClassLabels(suf_df[,ncol(suf_df)])
   df_trainandtest <- splitForTrainingAndTest(dfValues, dfTargets, ratio=0.1)
   norm_data<- normTrainingAndTestSet(df_trainandtest)
   print("training")
   model <- mlp(norm_data$inputsTrain, norm_data$targetsTrain, size=60, learnFunc="BackpropBatch", learnFuncParams=c(10, 0.1), maxit=2000, inputsTest=norm_data$inputsTest, targetsTest=norm_data$targetsTest)
   print("testing")
   predictions = predict(model,norm_data$inputsTest)
   #生成混淆矩阵，观察预测精�?
   print("predict result")
   confusionMatrix(norm_data$targetsTest,predictions)
}
cla=df$cla
ran=df$ran
write.csv(cla,file=paste0(rootpath,"/Expdata.csv"))
write.csv(ran,file=paste0(rootpath,"/ranExpdata.csv"))