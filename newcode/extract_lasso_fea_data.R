setwd("D:/WORK/样本多分类/对比/kmeans/lasso/example/STAD_latter")
datapath='STAD_all.csv'
data=read.csv(datapath)
select_data_lasso<-function(fea_num,data){
  feapath=paste0(fea_num,'features.csv')
  features=read.csv(feapath,header = 0)
  features_gene=features[,1]
  lable=data$all_label
  select=data[,features_gene]
  Xselect=cbind(select,lable)
  print(dim(select))
  write.csv(Xselect,file=paste0(fea_num,"_Xselect.csv"),row.names = F)
}
fea_num_list=c(20,50,100,200)
for(i in fea_num_list){
  print(i)
  select_data_lasso(i,data)
}

