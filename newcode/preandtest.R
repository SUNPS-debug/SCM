#输入一个矩阵
pick_edge<-function(P_mat,dif_gene,cut){
  sub_p=P_mat[dif_gene,dif_gene]
  numrow=length(dif_gene)
  sub_p[upper.tri(sub_p)]=0
  diag(sub_p)=rep(0,numrow)
  sub_p=abs(sub_p)
  #sub_p[sub_p<=cut]<-0
  sub_p[sub_p>cut]<-0
  node1_index=list()
  node2_index=list()
  for(i in 1:numrow){
    node2_index[[i]]=which(sub_p[i,]>0)
    if(length(node2_index[[i]])!=0){
      node1_index[[i]]=rep(i,length(node2_index[[i]]))
    }
  }
  node1=unlist(node1_index)
  node2=unlist(node2_index)
  edge_loc=cbind(node1,node2)
  return(edge_loc)
}
#根据dif_gene取子mat，返回值为5*5的list
caculate_dif_sim<-function(mat,dif_gene_set,cutoff){
  mat_cut_all=list()
  for(i in 1:ntype){
    mat_cut=list()
    for(j in 1:ntype){
      submat=mat[[i]][dif_gene_set[[j]],dif_gene_set[[j]]]
      #submat[submat>cut]<-0
      #若mat为cor
      submat[abs(submat)<cutoff]<-0 
      submat[upper.tri(submat)]=0
      diag(submat)=rep(0,nrow(submat))
      mat_cut[[j]]=submat
    }
    mat_cut_all[[i]]=mat_cut
  }
  return(mat_cut_all)
}

max_edge<-function(mat_cut_all,cut_index){
    remain_edges=list()
    submat=list()
    for(i in 1:ntype){
       mod_size=length(dif_gene_set[[i]])
      submat[[i]]=mat_cut_all[[i]][[i]]
      edges_sort=sort(abs(submat[[i]]),decreasing = TRUE)
      cut=edges_sort[cut_index]
      submat[[i]][abs(submat[[i]])<cut]<-0
    }
    remain_edges_loc=pick_edge_loc(submat)
    remain_edges=list(mat=submat,edges_loc=remain_edges_loc)
    return(remain_edges)
}

check<-function(cor_cut_all,remain_edge_loc){
  subcor=list()
  edges_re_all=list()
  classn=rep(1:ntype)
  for(i in 1:ntype){
    edges_re=list()
    subcor[[i]]=cor_cut_all[[i]][[i]]
    edges1=subcor[[i]][remain_edge_loc[[i]]]
    for(j in classn[-i]){
      othercor=cor_cut_all[[j]][[i]]
      edges2=othercor[remain_edge_loc[[i]]]
      edges_re[j]=t.test(subcor[[i]],othercor,alternative ="greater")$p.value
    }
    edges_re_all[[i]]=edges_re
  }
  return(edges_re_all)
}

#输入：p_cut_all，输出：每个模块中保留特异边的矩阵
pick_edge_all<-function(mat_cut_all){
  cut_mat=list()
  class_n=rep(1:ntype)
  for(i in 1:ntype){
    sub_mat=list()
    for(j in 1:ntype){
      sub_mat_p=mat_cut_all[[j]][[i]]
      sub_mat_p[sub_mat_p!=0]<-1
      sub_mat[[j]]=sub_mat_p
    }
    mid_mat=sub_mat[[i]]
    for(k in class_n[-i]){
      mid_mat=mid_mat-sub_mat[[k]]
    }
    mid_mat[mid_mat<=0]<-0
    cut_mat[[i]]=mid_mat
  }
  return(cut_mat)
}
#记录特异性边所在位置，即边的值不为0的边，输出为边的坐标
pick_edge_loc<-function(mat_set){
  edge_loc_set=list()
  for(i in 1:ntype){
    mod_size=dim(mat_set[[i]])[1]
    node1_index=list()
    node2_index=list()
    for(j in 1:mod_size){
      node2_index[[j]]=which(mat_set[[i]][j,]!=0)
      if(length(node2_index[[j]])!=0){
        node1_index[[j]]=rep(j,length(node2_index[[j]]))
      }
    }
    node1=unlist(node1_index)
    node2=unlist(node2_index)
    edge_loc=cbind(node1,node2)
    edge_loc_set[[i]]=edge_loc
  }
  return(edge_loc_set)
}
dif_cor<-function(meancor,dif_gene_set){
  subcor=list()
  for (i in 1:ntype){
    subcor[[i]]= meancor[[i]][dif_gene_set[[i]],dif_gene_set[[i]]]
    subcor[[i]][upper.tri(subcor[[i]])]=0
    diag(subcor[[i]])=rep(0,nrow(subcor[[i]]))
  }
 return(subcor)
}

#找值最大的一系列特异性边
filter_edge<-function(subcor,edge_loc_set,cut_index){
  remain_edges=list()
   for(i in 1:ntype){
    mod_size=length(dif_gene_set[[i]])
    edge_num=dim(edge_loc_set[[i]])[1]
    edges=subcor[[i]][edge_loc_set[[i]]]
    edges_sort=sort(abs(edges),decreasing = TRUE)
    if(edge_num>=cut_index){cut=edges_sort[cut_index]}
    if(edge_num<cut_index) {cut=edges_sort[edge_num]}
    subcor[[i]][abs(subcor[[i]])<cut]<-0
    remain_edges[[i]]=subcor[[i]]
   }
  return(remain_edges)
}

select_edges<-function(remain_edges,p_cut_mat){
  max_edge=list()
  for(i in 1:ntype){
    max_edge[[i]]=remain_edges[[i]]*p_cut_mat[[i]]
  }
  return(max_edge)
}
#将原始数据分为测试集和训练集
trainandtest<-function(cleandata){
  set.seed(5)
  traindata=list()
  testdata=list()
  for(i in 1:ntype){
    data=cleandata[[i]]
    num_sample=nrow(data)
    index <- sample(1:num_sample, size=0.8*num_sample)
    traindata[[i]]<-data[index,]
    testdata[[i]]<-data[-index,]
  }
  train_test=list(traindata=traindata,testdata=testdata)
  return(train_test)
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
#输入一个样本判断其类型
predict<-function(pre_sample,traindata,train_sim_mat,max_edge_loc){
  new_edge=list()
  delta=array()
  edge_delta=array()
  train_edge=list()
  for(i in 1:ntype){
    nk=dim(traindata[[i]])[1]
    new_mat=rbind(traindata[[i]],pre_sample)
    new_mat=new_mat[,dif_gene_set[[i]]]
    new_net=cor(new_mat,method = "spearman")
    new_edge[[i]]=new_net[max_edge_loc[[i]]]
    train_edge[[i]]=train_sim_mat[[i]][max_edge_loc[[i]]]
    delta[i]=sum(abs(abs(train_edge[[i]])-abs(new_edge[[i]])))
    delta[i]=delta[i]*(nk^(2/3))
    edge_num=nrow(max_edge_loc[[i]])
    edge_delta[[i]]=delta[i]/edge_num
  }
  min_index=which.min(edge_delta)
  pre_result=list(edge_delta=edge_delta,pre=min_index)
  return(pre_result)
}
predict_all<-function(){
  pre_samples=do.call(rbind,testdata)
  sam_num=dim(pre_samples)[1]
  pre_type=list()
  for(j in 1:sam_num){
    print(j)
    pre_sample=pre_samples[j,]
    pre_result=predict(pre_sample,traindata,train_sim_mat,max_edges_loc)
    pre_type[[j]]=pre_result
    print(paste0("result:",pre_type[[j]]$pre))
  }
  return(pre_type)
}

max_edges_loc=list()
for(i in 1:ntype){
  max_edges_loc[[i]]=max_edges_list[[i]][[dif_mod[[i]]]]
}

############################################################
cor_cut_dif=caculate_dif_sim(meancor,dif_gene_set,0.6)

max_edges=max_edge(cor_cut_dif,100)
max_edges_loc=max_edges$edges_loc



cor_cut_mat=pick_edge_all(cor_cut_dif)
edge_loc_set=pick_edge_loc(cor_cut_mat)

dif_cor_set=dif_cor(meancor,dif_gene_set)  
#选取边的值最大的100条边
filter_edge_set=filter_edge(dif_cor_set,edge_loc_set,100) 
#选取特异边中最大的100条
max_mat=select_edges(filter_edge_set,cor_cut_mat)
#记录特异边中最大的100条的位置
max_edges_loc=pick_edge_loc(max_mat)
train_test=trainandtest(new_cleandata)
traindata=train_test$traindata
testdata=train_test$testdata
train_sim_mat=sim_train(dif_gene_set,traindata)
max_edges_loc=spe_edges_loc

pre_type=predict_all()
pre_ty=array()
for(i in 1:length(pre_type)){
  pre_ty[i]=pre_type[[i]]$pre
}
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