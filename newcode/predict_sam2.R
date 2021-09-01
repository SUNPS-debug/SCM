sampling<-function(sub_cleandata,num_sam,fold,typei){
  re_index=list()
  nsample=nrow(sub_cleandata)
  need_g=ceiling((fold*num_sam)/nsample)
  for(i in 1:need_g){
    set.seed(i)
    re_index[[i]]=sample(1:nsample,nsample)
  }
  print(re_index)
  new_order=unlist(re_index)
  sampling_index=list()
  samplingData=list()
  sampling_point=seq(1,length(new_order),num_sam)
  #sampling_data_dir=paste0(rootpath,"/samplingdata/",cancertype[typei])
  for(j in 1:fold){
    sampling_index[[j]]=new_order[sampling_point[j]:(sampling_point[j+1]-1)]
    samplingData[[j]]=sub_cleandata[sampling_index[[j]],]
    #   sampling_data_loc=paste0(sampling_data_dir,"/fold",j,".csv")
    #   write.csv(samplingData[[j]],sampling_data_loc)
  }
  return(samplingData)
}
sampling_all<-function(new_cleandata,num_list){
  min_sample=min(num_list)
  max_sample=max(num_list)
  if(max_sample/min_sample<=10){fold=10}
  if(max_sample/min_sample>10) {fold=ceiling(max_sample/min_sample)}
  num_sam=min_sample-1
  sampling_data_all=list()
  for(i in 1:ntype){
    sampling_data_all[[i]]=sampling(new_cleandata[[i]],num_sam,fold,i)
  }
  return(sampling_data_all) 
}

predict_sampling<-function(pre_sample,dif_gene_set,
                           sampling_train,sim_sam_train_all,
                           spe_edges_loc){
  new_edge_all=list()
  train_edge_all=list()
  delta_all=list()
  for(ty in 1:ntype){
    new_edge=list()
    train_edge=list()
    delta=array()
    fold=length(sampling_train[[ty]])
    for(f in 1:fold){
      new_mat=rbind(sampling_train[[ty]][[f]],pre_sample)
      new_mat=new_mat[,dif_gene_set[[ty]]]
      new_net=cor(new_mat,method = "spearman")
      new_edge[[f]]=new_net[spe_edges_loc[[ty]]]
      train_edge[[f]]=sim_sam_train_all[[ty]][[f]][spe_edges_loc[[ty]]]
      #delta[f]=mean(abs(abs(train_edge[[f]])-abs(new_edge[[f]])))
      delta[f]=mean(abs(train_edge[[f]]-new_edge[[f]]))
    }
    new_edge_all[[ty]]=new_edge
    train_edge_all[[ty]]=train_edge
    delta_all[[ty]]=delta
  }
   return(delta_all)
}
  


predict_all2<-function(traindata,testdata,spe_edges_loc){
  train_num=array()
  sim_sam_train_all=list()
  for(i in 1:ntype){
    train_num[i]=dim(traindata[[i]])[1]
    print(paste0("train_num",train_num[i]))
  }
  sampling_train=sampling_all(traindata,train_num)
  for(k in 1:ntype){
    sim_sam_train=list()
    fold=length(sampling_train[[k]])
    for(j in 1:fold){
      sub_sam_train=sampling_train[[k]][[j]][,dif_gene_set[[k]]]
      sim_sam_train[[j]]=cor(sub_sam_train,method = "spearman")
    }
    sim_sam_train_all[[k]]=sim_sam_train
  }
  pre_samples=do.call(rbind,testdata)
  sam_num=dim(pre_samples)[1]
  pre_type=array()
    for(sam in 1:sam_num){
      print(paste0("sample:",sam))
      delta_sum=array()
      presample=pre_samples[sam,]
      print(rownames(pre_samples)[sam])
      delta_mat=predict_sampling(presample,dif_gene_set,
                                       sampling_train,sim_sam_train_all,
                                       spe_edges_loc)
      
      for(t in 1:ntype){
        delta_sum[t]=sum(delta_mat[[t]])
      }
      print(delta_sum)
      pre_type[sam]=which.min(delta_sum)
    }
  return(pre_type)
  }
  
  
         

pre_result=predict_all2(traindata,testdata,spe_edges_loc)
pre_type=pre_result
test_num=array()
for(i in 1:ntype){
  test_sam=testdata[[i]]  
  test_num[i]=nrow(test_sam)
}
true_type=rep(1:ntype,test_num)

true_type_loc=paste0(rootpath,"/truetype.csv")
pre_type_loc=paste0(rootpath,"/predicttype.csv")
write.csv(true_type,true_type_loc)
write.csv(pre_type,pre_type_loc)