file_path="D:/WORK/样本多分类/BRCAsubtype/trainset/ORtrainandtest/train_data1.csv"
fea_select<-function(file_path){
  train_data=read.csv(file_path,check.names = F,row.names = NULL)
  shape=dim(train_data)
  data=train_data[,-1]
  data.shape=dim(data)
  label=data[,data.shape[2]]
}


fea_select_list=c("lasso","mutual_info_classif","f_classif","chi2","co-ex")
check_commengene<-function(fea_select_list,fea_num){
  fea_select_num=length(fea_select_list)
  common_gene=matrix(0,nrow = fea_select_num,ncol = fea_select_num,
                     dimnames=list(fea_select_list,fea_select_list))
  feature_list=list()
  for(i in 1:fea_select_num){
    print(fea_select_list[i])
    feature_path=paste0("C:/Users/孙佩硕的magic/Desktop/STAD_latter/",fea_select_list[i],
                    "/feature",fea_num,".csv")
    feature_file=read.csv(feature_path,header = 0)
    feature_list[[i]]=feature_file[,1]
  }
  for(a in 1:fea_select_num){
    for(b in 1:fea_select_num){
      com_ge=intersect(feature_list[[a]],feature_list[[b]])
      common_gene[a,b]=length(com_ge)
    }
  }
  out_path=paste0("C:/Users/孙佩硕的magic/Desktop/BRCA/common_gene_ana/",fea_num,"features.csv")
  write.csv(common_gene,out_path)
  return(common_gene)
}
check_commengene(fea_select_list,200)

