predict_sampling<-function(pre_sample,traindata,dif_gene_set,spe_edges_loc,train_num){
  new_edge=list()
  train_edge=list()
  delta=array()
  sam_train=list()
  sam_len=min(train_num)-1
  for(j in 1:ntype){
    data=traindata[[j]]
    index <- sample(1:train_num[j], size=sam_len)
    sam_train[[j]]=data[index,]
  }
  sam_train_cor=sim_train(dif_gene_set,sam_train)
  print("sam_train_cor finished")
  for(k in 1:ntype){
    print(paste0("testtype:",k))
    new_mat=rbind(sam_train[[k]],pre_sample)
    new_mat=new_mat[,dif_gene_set[[k]]]
    new_net=cor(new_mat,method = "spearman")
    new_edge[[k]]=new_net[spe_edges_loc[[k]]]
    train_edge[[k]]=sam_train_cor[[k]][spe_edges_loc[[k]]]
    delta[k]=mean(abs(train_edge[[k]]-new_edge[[k]]))
    #delta[k]=mean(abs(abs(train_edge[[k]])-abs(new_edge[[k]])))
  }
  return(delta)
}

predict_all2<-function(sam_times){
  train_num=array()
  for(i in 1:ntype){
  train_num[i]=dim(traindata[[i]])[1]
  print(paste0("train_num",train_num[i]))
  }
  pre_samples=do.call(rbind,testdata)
  sam_num=dim(pre_samples)[1]
  pre_type=array()
  delta_mat_list=list()
  for(time in 1:sam_times){
    print(paste0("time:",time))
    delta_mat=matrix(0,nrow=sam_num,ncol=ntype)
    for(sam in 1:sam_num){
      print(paste0("sample:",sam))
      presample=pre_samples[sam,]
      delta_mat[sam,]=predict_sampling(presample,traindata,
                                       dif_gene_set,
                                       spe_edges_loc,train_num)
    }
    delta_mat_list[[time]]=delta_mat
  }
  sum_delta=Reduce("+",delta_mat_list)
  pre_type=apply(sum_delta,1,which.min)
  pre_result=list(sum_delta=sum_delta,pre_type=pre_type)
  return(pre_result)
}
pre_result=predict_all2(12)
pre_type=pre_result$pre_type
for(i in 1:ntype){
  test_sam=testdata[[i]]
  test_num[i]=nrow(test_sam)
}
true_type=rep(1:ntype,test_num)

true_type_loc=paste0(rootpath,"/truetype.csv")
pre_type_loc=paste0(rootpath,"/predicttype.csv")
write.csv(true_type,true_type_loc)
write.csv(pre_type,pre_type_loc)