library(mvoutlier)
load("D:/WORK/样本多分类/test/trainset.RData")
load("D:/WORK/样本多分类/test/testset.RData")
load("D:/WORK/样本多分类/test/dif_gene_set.RData")
load("D:/WORK/样本多分类/test/remain_edge_set.RData")


check_by_ma<-function(presample,num_list){
  ma_times=array()
  for(j in 1:ntype){
    data=rbind(train_dataset[[j]],presample)
    pre_index=num_list[j]+1
    pre_index=as.character(pre_index)
    genelist=colnames(data)
    outliers_score=list()
    for(i in 1:100){
      gene1=dif_gene_set[[j]][remain_edge_set[[j]][i,1]]
      gene2=dif_gene_set[[j]][remain_edge_set[[j]][i,2]]
      index1=which(genelist==gene1)
      index2=which(genelist==gene2)
      X1=data[,index1]
      X2=data[,index2]
      z=cbind(X1,X2)
      res1 <- uni.plot(z,quan=3/4,alpha=0.025)
      outliers_score[[i]]=which(res1$outliers==T)
    }
    outliers=unlist(outliers_score)
    aa=names(table(outliers))
    pre_loc=which(aa==pre_index)
    if(length(pre_loc)==0){ma_times[j]=0}
    else{ma_times[j]=table(outliers)[pre_loc]}
  }
  return(ma_times)
}
check_by_ma_all<-function(test_dataset){
  test_num=array()
  pre_samples=data.frame()
  for(j in 1:ntype){
    test_num[j]=nrow(test_dataset[[j]])
    pre_samples=rbind.data.frame(pre_samples,test_dataset[[j]])
  }
  #pre_samples=do.call(rbind,test_dataset)
  true_type=rep(c(1:ntype),test_num)
  ma_outliers=matrix(nrow =length(true_type),ncol=ntype)
  for(i in 1:length(true_type)){
    print(i)
    pre_sample=pre_samples[i,]
    ma_outliers[i,]=check_by_ma(pre_sample,train_num_list)
  }
  return(ma_outliers)
}
ma_outliers=check_by_ma_all(test_dataset)
minout_times=apply(ma_outliers,1,min) 
minout_index=list()
for(i in 1:length(minout_times)){
  minout_index[[i]]=which(ma_outliers[i,]==minout_times[i])
}
txtname=names(dif_gene_set)
for(i in 1:5){
  outpath=paste0("D:/WORK/样本多分类/test/",txtname,".txt")
  kk=dif_gene_set[[i]]
  write.table(kk,outpath[i],row.names =F,quote=F ,col.names = F)
}