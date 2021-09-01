cvlist=list()
for(cv in 1:ntype){
  datasize=nrow(new_cleandata[[cv]])
  cvlist[[cv]]=CVgroup(10,datasize,cv)
}
sim_train<-function(dif_gene_set,traindata){
  sim_mat=list()
  for(i in 1:ntype){
    data=traindata[[i]][,dif_gene_set[[i]]]
    sim_mat[[i]]=cor(data,method = "spearman")
    rownames(sim_mat[[i]])=dif_gene_set[[i]]
    colnames(sim_mat[[i]])=dif_gene_set[[i]]
  }
  return(sim_mat)
} 


#输入一个样本，输出在每类特征边中的扰动值，暂定为一个向量
feautre_edge_value<-function(pre_sample,dif_gene_set,
                             ref_data,ref_net,spe_edges_loc){
  new_edge=list()
  delta=list()
  edge_delta=array()
  ref_edge=list()
  for(i in 1:ntype){
    new_mat=rbind(ref_data[[i]],pre_sample)
    new_mat=new_mat[,dif_gene_set[[i]]]
    new_net=cor(new_mat,method = "spearman")
    new_edge[[i]]=new_net[spe_edges_loc[[i]]]
    ref_edge[[i]]=ref_net[[i]][spe_edges_loc[[i]]]
    delta[[i]]=abs(abs(ref_edge[[i]])-abs(new_edge[[i]]))
  }
  all_feature=unlist(delta)
  return(all_feature)
}
generate_ref_set<-function(traindata,train_num,ref_size){
  ref_data=list()
  ref_index=list()
  for(re in 1:ntype){
    ref_index[[re]]=sample(train_num[re],size=ref_size)
    ref_data[[re]]=traindata[[re]][ref_index[[re]],]
  }
  return(ref_data)
}
# generate_train_set<-function(traindata,ref_size,dif_gene_set,spe_edges_loc){
#   train_ref=generate_ref_set(traindata,ref_size)
#   train_num=train_ref$all_num
#   ref_data=train_ref$ref_data
#   ref_net=sim_train(dif_gene_set,ref_data)
#   train_num_sum=sum(train_num)
#   train_target=rep(cancertype,train_num)
#   train_all=do.call(rbind,traindata)
#   train_sams_name=rownames(train_all)
#   fea_num=nrow(spe_edges_loc[[1]])
#   edge_feature_train=matrix(0,ncol=ntype*fea_num,nrow=train_num_sum)
#   rownames(edge_feature_train)=train_sams_name
#   for(sa in 1:train_num_sum){
#     print(paste0("training_sample:",sa))
#     pre_sam=train_all[sa,]
#     dis_value=feautre_edge_value(pre_sam,dif_gene_set,ref_data,ref_net,spe_edges_loc)
#     edge_feature_train[sa,]=dis_value
#   }
#   train_csv=cbind(edge_feature_train,train_target)
#   return(train_csv)
# }
#一组ref生成traindata,返回matrix形式，
generate_train_part<-function(train_ref_part,train_data_part,fea_num,train_num,ty){
  ref_net=sim_train(dif_gene_set,train_ref_part)
  train_name=rownames(train_data_part)
  edge_feature_train=matrix(0,ncol=ntype*fea_num,nrow=train_num[ty])
  rownames(edge_feature_train)=train_name
  for(sa in 1:train_num[ty]){
    print(paste0("type:",ty,"#training_sample:",sa))
    pre_sam=train_data_part[sa,]
    dis_value=feautre_edge_value(pre_sam,dif_gene_set,train_ref_part,ref_net,spe_edges_loc)
    edge_feature_train[sa,]=dis_value
  }
  return(edge_feature_train)
}

generate_train_set<-function(traindata,ref_size,base_num,dif_gene_set,spe_edges_loc){
  train_num=sapply(traindata,nrow)
  ref_num=ceiling(base_num/min(train_num))  
  train_ref=list()
  fea_num=nrow(spe_edges_loc[[1]])
  for(re_t in 1:ref_num){
    train_ref[[re_t]]=generate_ref_set(traindata,train_num,ref_size)
  }
  train_target=list()
  train_data=list()
  for(ty in 1:ntype){
    epoch=ceiling(base_num/train_num[ty])
    train_target[[ty]]=rep(cancertype[ty],train_num[ty]*epoch)
    train_data_part=list()
    for(ep in 1:epoch){
      train_data_part[[ep]]=generate_train_part(train_ref[[ep]],traindata[[ty]],fea_num,train_num,ty)
    }
    train_data_value=do.call(rbind,train_data_part)
    train_data[[ty]]=cbind(train_data_value,train_target[[ty]])
  }
  train_set=do.call(rbind,train_data)
  return(train_set)
}
 
 
generate_test_set<-function(testdata,traindata,ref_size,te_epoch,dif_gene_set,spe_edges_loc){
  test_num=sapply(testdata,nrow)
  test_num_sum=sum(test_num)
  test_target=rep(cancertype,test_num)
  train_num=sapply(traindata,nrow)
  test_all=do.call(rbind,testdata)
  test_sams_name=rownames(test_all)
  fea_num=nrow(spe_edges_loc[[1]])
  test_ref_data=list()
  test_ref_sim=list()
  test_csv=list()
  for(te_ep in 1:te_epoch){
    test_ref_data[[te_ep]]=generate_ref_set(traindata,train_num,ref_size)
    test_ref_sim[[te_ep]]=sim_train(dif_gene_set,test_ref_data[[te_ep]])
    edge_feature_test=matrix(0,ncol=ntype*fea_num,nrow=test_num_sum)
    rownames(edge_feature_test)=test_sams_name
    for(te in 1:test_num_sum){
      print(paste0("test_sample:",te))
      te_sam=test_all[te,]
      test_dis=feautre_edge_value(te_sam,dif_gene_set,test_ref_data[[te_ep]],test_ref_sim[[te_ep]],spe_edges_loc)
      edge_feature_test[te,]=test_dis
    }
    test_csv[[te_ep]]=cbind(edge_feature_test,test_target)
  }
  return(test_csv)  
}
  

fea_dir=paste0(rootpath,"/50feature")
  if(dir.exists(fea_dir)==F){
    dir.create(fea_dir)
  }

juduge <- function(x){
  if(x%%2 ==0){epoch=x+1}
  else{epoch=x+2}
  return(epoch)
}

for(pretime in 1:10){
  print(paste0("交叉验证第",pretime,"轮"))
  train_test=trainandtest(new_cleandata,cvlist,pretime)
  traindata=train_test$traindata
  testdata=train_test$testdata
  train_dir=paste0(fea_dir,"/trainset")
  if(dir.exists(train_dir)==F){
    dir.create(train_dir)
  }
  train_dataset=generate_train_set(traindata,20,200,dif_gene_set,spe_edges_loc)
  #train_loc=paste0(rootpath,"/trainset/train_data",pretime,".csv")
  train_loc=paste0(train_dir,"/train_data",pretime,".csv")
  write.csv(train_dataset,train_loc)

  ##生成test数据
  

  test_epoch=juduge(ntype)
  test_samples=generate_test_set(testdata,traindata,20,test_epoch,
                                 dif_gene_set,spe_edges_loc)
  for(te_dir in 1:test_epoch){
    test_dir=paste0(fea_dir,"/testset",te_dir)
    if(dir.exists(test_dir)==F){
      dir.create(test_dir)
  }
  test_loc=paste0(test_dir,"/test_data",pretime,".csv")
  write.csv(test_samples[[te_dir]],test_loc)
  }
  print(paste0("第",pretime,"轮数据生成完毕"))
}
  
  


