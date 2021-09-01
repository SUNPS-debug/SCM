library(WGCNA)
library(reshape2)
library(stringr)
library(plyr)
library(caret)
#在rootpath下新建文件夹
create_dir<-function(rootpath){
  dirname=c("samplingdata","cleandata","trainset","testset","data",cancertype)
  create_dir=paste0(rootpath,"/",dirname)
  n=length(dirname)
  for(i in 1:length(dirname)){
    if(dir.exists(create_dir[i])==F){
      dir.create(create_dir[i]) 
    }
  }
}
#将txt数据转换为csv 
txt_to_csv<-function(txtpath){
  genename=genename[,1]
  OR_dir_txt=paste0(txtpath,"/",txtname)
  OR_dir_csv=paste0(ORpath,"/",cancertype,".csv")
  for(i in 1:length(txtname)){
    txtdata=read.table(OR_dir_txt[i],check.names = F)
    rownames(txtdata)=genename
    write.csv(txtdata,OR_dir_csv[i])
  }
}
#读入数据并记录每类样本数
data_num <- function(ORpath) {                
  numlist = rep(0, ntype)
  datapath = paste0(ORpath, "/",cancertype,".csv")
  data = list()
  for (i in 1:ntype) {
    data[[i]]= read.csv(datapath[i], row.names = 1, check.names = F)
    numlist[i] <- ncol(data[[i]])
    print(numlist[i])
  }
  data_num=list(data=data,numlist=numlist)
  return(data_num)
}
#清理数据
datacleaning <- function(data,ORdata,cut) {
  cleandata=list()
  mean_gene=apply(data,1,mean)
  max_gene=which(mean_gene>=cut)
  dataExpr=data[max_gene,]
  for(or in 1:ntype){
    ORdata[[or]]=ORdata[[or]][max_gene,]
  }
  m.mad <- apply(dataExpr, 1, mad)
  remain_index=which(m.mad >max(quantile(m.mad, probs = seq(0, 1, 0.1))[2], 0))
  for(cl in 1:ntype){
    tmp=ORdata[[cl]][remain_index,]
    cleandata[[cl]]=log(tmp+1)
    cleandata[[cl]]=as.data.frame(t(cleandata[[cl]]))
    print(dim(cleandata[[cl]]))
  }
  return(cleandata)
   #dataExprVar <- dataExpr[remain_index, ]
  #dataExpr <- as.data.frame(t(dataExprVar))
  # gsg = goodSamplesGenes(dataExpr, verbose = 3)
  # if (!gsg$allOK) {
  #   # Optionally, print the gene and sample names that were removed:
  #   if (sum(!gsg$goodGenes) > 0)
  #     printFlush(paste("Removing genes:",
  #                      paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")))
  #   
  #   if (sum(!gsg$goodSamples) > 0)
  #     printFlush(paste("Removing samples:",
  #                      paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")))
  #   
  #   # Remove the offending genes and samples from the data:
  #   dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
  #   print(dim(dataExpr))
  }


#求五个数据集里共同的基因，数据行为样本列名为基因
common_gene<-function(datalist){
  genelist=list()
  for(i in 1:ntype){
    genelist[[i]]=colnames(datalist[[i]])
  }
  common_genelist=Reduce(intersect,genelist)
  print(length(common_genelist))
  return(common_genelist)
}
filter_cleandata<-function(cleandata,common_genelist){
  new_cleandata=list()
  for (i in 1:ntype) {
    #common_gene_index=which(colnames(cleandata[[i]])==common_genelist)
    new_cleandata[[i]]=cleandata[[i]][,common_genelist]
  }
  return(new_cleandata)
}
#对一类样本进行采样
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
#对所有类型的样本进行采样
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
#计算平均的cor和p矩阵
mean_corandp<-function(sub_sampling_data,genenum){
  new_cor=matrix(0,nrow = genenum,ncol = genenum )
  new_p=matrix(1,nrow = genenum,ncol = genenum )
  fold=length(sub_sampling_data)
  for(i in 1:fold){
    data=sub_sampling_data[[i]]
    nsamples=nrow(data)
    mat_cor=cor(data,method = "spearman")
    mat_p=corPvalueStudent(mat_cor,nsamples)
    new_cor=new_cor+mat_cor
    new_p=new_p*mat_p
    print(i)
  }
  mean_cor=new_cor/fold
  mean_p=new_p^(1/fold)
  mean_corp=list(mean_cor=mean_cor,mean_p=mean_p)
  return(mean_corp)
}
caculate_power<-function(Expdata,nSamples,type){
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(Expdata, RsquaredCut = 0.8, 
                          powerVector=powers, 
                          networkType="unsigned", verbose=5)
  power = sft$powerEstimate
  if (is.na(power)){
    power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                   ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                 ifelse(type == "unsigned", 6, 12))       
                   )
    )
  }
  return(power)
}

#根据meancor生成模块
generate_module2<-function(sub_meancor,power,datExpr,classn){
  adj=(abs(sub_meancor))^power
  TOM=TOMsimilarity(adj,TOMType = "unsigned")
  dissTOM=1-TOM
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = 10
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = FALSE, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  names(dynamicMods)=rownames(adj)
  #dynamicColors = labels2colors(dynamicMods)
  ### 合并表达图谱相似的模块
  #计算模块的主成分
  MEList = moduleEigengenes(datExpr, colors = dynamicMods,
                            nPC = 2, softPower = 6,)
  MEs = MEList$eigengenes
  varExp=MEList$varExplained
  npcpath=paste0(rootpath,"/",cancertype[classn],"/PCA.csv")
  write.csv(varExp,npcpath)
  
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  MEDissThres =0
  
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr,dynamicMods, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors
  # Eigengenes of the new merged
  newMEs=merge$newMEs
  modules_loc=paste0(rootpath,"/",cancertype[classn],"/modules.csv")
  write.csv(mergedColors,modules_loc)
  Mod_info=list(newMEs=newMEs,Mod=mergedColors)
  return(Mod_info)
}
#将模块存储成list形式
modtolist<-function(modules){
  MODlist=list()
  for(i in 1:ntype){
    mod_list=list()
    modulename=modules[[i]]$Mod
    rename_mod=paste0("Mod",modulename)
    gene=names(modules[[i]]$Mod)
    modname=names(table(rename_mod))
    modnum=length(modname)
    for(j in 1:modnum){
      mod_index=which(rename_mod==modname[[j]])
      mod_list[[j]] =gene[mod_index]
    }
    names(mod_list)=modname
    MODlist[[i]]=mod_list
  }
  names(MODlist)=cancertype
  return(MODlist)
}
mod_overlap<-function(mod1,mod2){
  clusters1=names(mod1)
  clusters2=names(mod2)
  clus1_num=length(clusters1)
  clus2_num=length(clusters2)
  overlap=matrix(0,nrow = clus1_num,ncol = clus2_num)
  over_rate=matrix(0,nrow = clus1_num,ncol = clus2_num)
  score=matrix(0,nrow = clus1_num,ncol=1)
  over_num=matrix(0,nrow = clus1_num,ncol=1)
  for(i in 1:clus1_num){
    for(j in 1:clus2_num){
      overgenes=intersect(mod1[[i]],mod2[[j]])
      overlap[i,j]=length(overgenes)
      over_rate[i,j]=overlap[i,j]/length(mod1[[i]])
    }
    one_index=which(overlap[i,-1]==1)
    over_loc=which(overlap[i,-1]>0)
    score1=length(one_index)
    over_num[i,1]=length(over_loc)
    if(overlap[i,1]!=0){score[i,1]=score1+overlap[i,1]}
    else score[i,1]=score1
  }
  colnames(overlap)=clusters2
  colnames(over_rate)=clusters2
  #rownames(overlap)=clusters1
  over=list(overlap=overlap,over_rate=over_rate,score=score,over_num=over_num)
  return(over)
}
#计算所有期的overlap，共五个大矩阵,作用即合并
mod_overlap_all<-function(modlist){
  type_index=rep(1:ntype)
  overlap_list=list()
  overrate_list=list()
  score_list=list()
  over_num_list=list()
  overlap_loc=paste0(rootpath,"/",cancertype,"/overlap.csv")
  overrate_loc=paste0(rootpath,"/",cancertype,"/overrate.csv")
  for(i in 1:ntype){
    overlap=list()
    overrate=list()
    scores=list()
    over_num=list()
    for (j in type_index[-i]) {
      overmat=mod_overlap(modlist[[i]],modlist[[j]])
      overlap[[j]]=overmat$overlap
      overrate[[j]]=overmat$over_rate
      scores[[j]]=overmat$score
      over_num[[j]]=overmat$over_num
    }
    overlap_list[[i]]=do.call(cbind,overlap)
    overrate_list[[i]]=do.call(cbind,overrate)
    score_list[[i]]=do.call(cbind,scores)
    over_num_list[[i]]=do.call(cbind,over_num)
    rownames(overlap_list[[i]])=names(modlist[[i]])
    rownames(overrate_list[[i]])=names(modlist[[i]])
    rownames(score_list[[i]])=names(modlist[[i]])
    colnames(score_list[[i]])=type_index[-i]
    rownames(over_num_list[[i]])=names(modlist[[i]])
    colnames(over_num_list[[i]])=type_index[-i]
    write.csv(overlap_list[[i]],overlap_loc[i])
    write.csv(overrate_list[[i]],overrate_loc[i])
  }
  names(overlap_list)=cancertype
  names(overrate_list)=cancertype
  names(score_list)=cancertype
  names(over_num_list)=cancertype
  over_list=list(overlap=overlap_list,
                 overrate=overrate_list,
                 score=score_list,
                 over_num=over_num_list)
  return(over_list)
}


select_every_mod<-function(meancor,modlist){
  mat_cut_all=list()
  for(i in 1:ntype){
    modnum=length(modlist[[i]])
    part_mat=list()
    part_bymod_all=list()
    for(j in 2:modnum){
      print(paste0("当前第",i,"期，第",j-1,"个模块"))
      part_bymod=list()
      geneset=modlist[[i]][[j]]
      for(k in 1:ntype){
        submat=meancor[[k]][geneset,geneset]
        submat[upper.tri(submat)]=0
        diag(submat)=rep(0,nrow(submat))
        part_bymod[[k]]=submat
      }
      part_bymod_all[[j]]= part_bymod
    }
    names(part_bymod_all)=names(modlist[[i]])
    mat_cut_all[[i]]=part_bymod_all
  }
  return(mat_cut_all)
}
#记录矩阵中非0值的位置
edge_loc<-function(Amat){
  mod_size=dim(Amat)[1]
  node1_index=list()
  node2_index=list()
  for(j in 1:mod_size){
    node2_index[[j]]=which(Amat[j,]!=0)
    if(length(node2_index[[j]])!=0){
      node1_index[[j]]=rep(j,length(node2_index[[j]]))
    }
  }
  node1=unlist(node1_index)
  node2=unlist(node2_index)
  edges_loc=cbind(node1,node2)
  return(edges_loc)
}
#取最大的前多少条边
mod_max_edges_index<-function(mat_cut_all,modlist,cut_index){
  max_edges_all=list()
  for(i in 1:ntype){
    mod_num=length(modlist[[i]])
    max_edges=list()
    for(j in 2:mod_num){
      ORcor=mat_cut_all[[i]][[j]][[i]]
      sort_edge=sort(abs(ORcor),decreasing = T)
      cut=sort_edge[cut_index]
      ORcor[abs(ORcor)<cut]<-0
      max_edges[[j]]=edge_loc(ORcor)
    }
    names(max_edges)=names(modlist[[i]])
    max_edges_all[[i]]=max_edges
  }
  return(max_edges_all)
}
#取值大于某个数值的边
mod_max_edges_value<-function(mat_cut_all,modlist,cut){
  max_edges_all=list()
  for(i in 1:ntype){
    mod_num=length(modlist[[i]])
    max_edges=list()
    for(j in 2:mod_num){
      ORcor=mat_cut_all[[i]][[j]][[i]]
      ORcor[abs(ORcor)<cut]<-0
      max_edges[[j]]=edge_loc(ORcor)
    }
    names(max_edges)=names(modlist[[i]])
    max_edges_all[[i]]=max_edges
  }
  return(max_edges_all)
}
#用mean衡量模块
eva_mod_by_mean<-function(mat_cut_all,modlist,max_edges_all){
  mean_all=list()
  for(i in 1:ntype){
    mod_num=length(modlist[[i]])
    mat_mean_mod=matrix(0,nrow=mod_num,ncol = ntype)
    for(j in 2:mod_num){
      mat_mean=array()
      edges_loc=max_edges_all[[i]][[j]]
      edges_num=dim(edges_loc)[1]
      weight_mean=array()
      for(k in 1:ntype){
        mod_mats=mat_cut_all[[i]][[j]][[k]]
        edge_value=mod_mats[edges_loc]
        #std=sd(abs(edge_value))
        mat_mean[k]=mean(abs(edge_value))
        #weight_mean[k]=mat_mean/std
      }
      mat_mean_mod[j,]=mat_mean
      #mat_mean_mod[j,]=weight_mean
    }
    rownames(mat_mean_mod)=names(modlist[[i]])
    mean_all[[i]]=mat_mean_mod
  }
  return(mean_all)
}
sort_modby_mean<-function(mean_all,modlist){
  sort_mean_mod=list()
  classn=rep(1:ntype)
  for(i in 1:ntype){
    modnum=dim(mean_all[[i]])[1]
    mod_names=names(modlist[[i]])
    maxmean=array()
    deltamean=array()
    for(k in 1:modnum){
      maxmean[k]=max(mean_all[[i]][k,-i])
      deltamean[k]=mean_all[[i]][k,i]-maxmean[k]
    }
    order_index=order(-deltamean)
    sort_mean_mod[[i]]=mod_names[order_index]
  }
  return(sort_mean_mod)
}
#三个指标评估模块
mod_evaluate<-function(over_list,modlist,sort_mean_mod){
  score=over_list$score
  overrate=over_list$overrate
  dif_mod=list()
  dif_gene_set=list()
  for(i in 1:ntype){
    mod_names=rownames(score[[i]])
    num =nrow(score[[i]]) 
    evaluate1=array()
    evaluate2=array()
    for(j in 2:num){
      evaluate1[j-1]=min(score[[i]][j,])  
      evaluate2[j-1]=max(overrate[[i]][j,])
    }
    names(evaluate1)=mod_names[-1]
    names(evaluate2)=mod_names[-1]
    sort1= order(evaluate1,decreasing = TRUE)
    sort2= order(evaluate2,decreasing = FALSE)
    modorder1=names(evaluate1)[sort1]
    modorder2=names(evaluate2)[sort2]
    modorder3=sort_mean_mod[[i]]
    weight=array()
    for(k in 1:(num-1)){
      weight2=which(modorder1==modorder2[k])  
      weight3=which(modorder3==modorder2[k]) 
      weight[k]=k+0*weight2+weight3
    }
    names(weight)=modorder2
    minmod_index=which.min(weight)
    dif_mod_name=modorder2[minmod_index]
    print(dif_mod_name)
    dif_mod[[i]]=dif_mod_name
    dif_gene_set[[i]]=modlist[[i]][[dif_mod_name]]
  }
  dif=list(dif_mod=dif_mod,dif_gene_set=dif_gene_set,weight=weight)
  return(dif)
}
pick_dif<-function(sort_mean_mod){
  dif_mod=list()
  dif_gene_set=list()
  for(m in 1:ntype){
    sort_mod=sort_mean_mod[[m]]
    mod_num=length(sort_mod)
    find=FALSE
    mo=1
    while(find==FALSE){
      if(length(modlist[[m]][[sort_mod[mo]]])<=100){
         mo=mo+1
         find==FALSE
      }
      else{
        find=TRUE
        dif_mod[[m]]=sort_mod[[mo]]
        dif_gene_set[[m]]=modlist[[m]][[sort_mod[[mo]]]]
        }
    }
    print(dif_mod[[m]])
  }
  dif=list(dif_mod=dif_mod,dif_gene_set=dif_gene_set)
  return(dif)
}

save_dif_mod<-function(mat_cut_all,dif_mod,dif_gene_set){
  mod_edges=list()
  for(i in 1:ntype){
    mod_edges[[i]]=mat_cut_all[[i]][[dif_mod[[i]]]]
    for(j in 1:ntype){
      mod_edges[[i]][[j]]=mod_edges[[i]][[j]][dif_gene_set[[i]],dif_gene_set[[i]]]
    }
  }
  return(mod_edges)
}
calculate_delta<-function(mod_edges){
  classn=rep(1:ntype)
  delta_mat_all=list()
  for(i in 1:ntype){
    OR_edges=mod_edges[[i]][[i]]
    row_num=dim(OR_edges)[1]
    delta_mat=list()
    for(j in classn[-i]){
      delta_mat[[j]]=abs(OR_edges)-abs(mod_edges[[i]][[j]])
    }
    delta_mat_all[[i]]=delta_mat
  }
  return(delta_mat_all)
}
pick_min_delta<-function(delta_mat_all,dif_gene_set){
  min_delta_list=list()
  classn=rep(1:ntype)
  for(i in 1:ntype){
    print(paste0("type:",i))
    mod_size=length(dif_gene_set[[i]])
    min_delta=matrix(1,nrow=mod_size,ncol=mod_size)
    for(j in 2:mod_size){
      print(paste0("row:",j))
      for(k in 1:(j-1)){
        print(paste0("col:",j))
        for(l in classn[-i]){
          print(paste0("type2:",i))
          if(delta_mat_all[[i]][[l]][j,k]<min_delta[j,k])
            min_delta[j,k]=delta_mat_all[[i]][[l]][j,k]
        }
      }
    }
    min_delta[min_delta==1]<-0
    min_delta_list[[i]]=min_delta
  }
  return(min_delta_list)
}

pick_spe_edges<-function(min_delta_list,cut_index){
  spe_edges_loc=list()
  for(i in 1:ntype){
    delta_mat=min_delta_list[[i]]
    sort_spe_edge=sort(delta_mat,decreasing = T)
    cut=sort_spe_edge[cut_index]
    delta_mat[delta_mat<cut]<-0 
    spe_edges_loc[[i]]=edge_loc(delta_mat)
  }
  return(spe_edges_loc)
}
#交叉验证分组
CVgroup <- function(k,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]    #将数据分成K份，并生成的完成数据集n
  temp <- sample(n,datasize)   #把n打乱
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}
#对应于pretime，每次对一个类进行划分
trainandtest<-function(cleandata,cvlist,pretime){
  traindata=list()
  testdata=list()
  for(i in 1:ntype){
    data=cleandata[[i]]
    test_index=cvlist[[i]][[pretime]]
    testdata[[i]]<-data[test_index,]
    traindata[[i]]<-data[-test_index,]
  }
  train_test=list(traindata=traindata,testdata=testdata)
  return(train_test)
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
      delta[f]=mean(abs(abs(train_edge[[f]])-abs(new_edge[[f]])))
      #delta[f]=mean(abs(train_edge[[f]]-new_edge[[f]]))
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


main_classification<-function(rootpath){
  # 打开多线程
  options(stringsAsFactors = FALSE)
  enableWGCNAThreads()
  #如果原始数据是txt格式则需要转换成csv
  #txtpath=paste0(rootpath,"/BRCA_Subtypes")
  #设置数据存储目录
  ORpath=paste0(rootpath,"/data")
  txtname<<-list.files(ORpath)
  #ntype代表类别数
  ntype<<-length(txtname)
  #cancertype代表类别的名字
  tmp<-strsplit(txtname,split=".",fixed=TRUE)
  cancertype<<-unlist(lapply(tmp,head,1)) 
  #基因名需要替换
  genename_loc<<-paste0(rootpath,"/gene_symbol_list_exceptSLC35E2")
  genename<<-read.table(genename_loc)
  #新建所需文件夹
  create_dir(rootpath)
  #txt转换成csv
  #txt_to_csv(txtpath)
  data_mes=data_num(ORpath)
  #读入数据并记录每类的数目
  ORdata=data_mes$data
  num_list=data_mes$numlist
  for(i in 1:ntype){
    ORdata[[i]]=ORdata[[i]][-(1:29),]
  }
  #清理每一类数据
  all_dataset=do.call(cbind,ORdata)
  new_cleandata=datacleaning(all_dataset,ORdata,cut=0)
  #new_cleandata=sapply(ORdata,t)
  #cleandata=list()
  #取每个类剩余基因的交集
  # common_genelist=common_gene(cleandata)
  # new_cleandata=filter_cleandata(cleandata,common_genelist)
  genenum=ncol(new_cleandata[[1]])
  sampling_data_all=sampling_all(new_cleandata,num_list)
  sampling_loc=paste0(rootpath,"/samplingdata/sampling_data_all.RData")
  save(sampling_data_all,file=sampling_loc)
  meancor=list()
  meanp=list()
  modules=list()
  for(j in 1:ntype){
    print(paste0("type:",j))
    corandp=mean_corandp(sampling_data_all[[j]],genenum)
    meancor[[j]]=corandp$mean_cor
    meanp[[j]]=corandp$mean_p
    power=caculate_power(new_cleandata[[j]],num_list[j],type="unsigned")
    print(power)
    modules[[j]]=generate_module2(meancor[[j]],power,new_cleandata[[j]],j)
  }
  meancor_loc=paste0(rootpath,"/meancor.RData")
  meanp_loc=paste0(rootpath,"/meanp.RData")
  save(meancor,file=meancor_loc)
  save(meanp,file=meanp_loc)
  #将模块存储在list里
  modlist<<-modtolist(modules)
  overlap_list=mod_overlap_all(modlist)
  overlap_list_loc=paste0(rootpath,"/overlap_list.RData")
  save(overlap_list,file=overlap_list_loc)
  #将每个模块的cormat单独取出来
  cor_cut_all=select_every_mod(meancor,modlist)
  #记录每个模块中大于0.4边的位置
  max_edges_list=mod_max_edges(cor_cut_all,modlist,0.4)
  #计算每个模块中值较大的边在本类和其他类的mean,估算边的特异性程度
  meanlist=eva_mod_by_mean(cor_cut_all,modlist,max_edges_list)
  #根据mean 将mod排序
  sort_mod_mean=sort_modby_mean(meanlist,modlist)
  ###参照overlap和mean给mod排出优先级，选出最特异的模块
  dif=mod_evaluate(overlap_list,modlist,sort_mod_mean)
  #dif=pick_dif(sort_mod_mean)
  dif_mod=dif$dif_mod
  dif_gene_set=dif$dif_gene_set
  #将每一类对应的特异性模块提取
  mod_edges=save_dif_mod(cor_cut_all,dif_mod,dif_gene_set)
  #寻找每个类特异性的边
  delta_mat=calculate_delta(mod_edges)
  min_delta_list=pick_min_delta(delta_mat,dif_gene_set)
  spe_edges_loc=pick_spe_edges(min_delta_list,50)
  #分出train和test
  #进行预测过程
  cvlist=list()
  for(cv in 1:ntype){
    datasize=nrow(new_cleandata[[cv]])
    cvlist[[cv]]=CVgroup(10,datasize,cv)
  }
  pre_type=list()
  true_type=list()
  
  for(pretime in 1:10){
    print(paste0("交叉验证第",pretime,"轮"))
    train_test=trainandtest(new_cleandata,cvlist,pretime)
    traindata=train_test$traindata
    testdata=train_test$testdata
    pre_result=predict_all2(traindata,testdata,spe_edges_loc)
    pre_type[[pretime]]=pre_result
    #pre_type_unlist=pre_result
    test_num=sapply(testdata,nrow)
    true_type[[pretime]]=rep(1:ntype,test_num)
    #true_type_unlist=rep(1:ntype,test_num)
  }
  pre_type_unlist=unlist(pre_type)
  true_type_unlist=unlist(true_type)
  pre_true=cbind(pre_type=pre_type_unlist,true_type=true_type_unlist)
  pre_true_loc=paste0(rootpath,"/pre_true.csv")
  write.csv(pre_true,pre_true_loc)
  con_mat=table(true_type_unlist,pre_type_unlist)
  
  result=list(true=true_type_unlist,
              predict=pre_type_unlist,
              confusionMatrix=con_mat
              )
  return(result)
}

result_analysis<-function(result){
  reference=result$true
  con.mat=result$confusionMatrix
  pre_eva=caret::confusionMatrix(con.mat)
  ana_table=pre_eva$byClass
  class_num=nrow(ana_table)
  eva_stan=colnames(ana_table)
  F1_index=which(eva_stan=="F1")
  precision_index=which(eva_stan=="Precision")
  F1_score=ana_table[,F1_index]
  precision=ana_table[,precision_index]
  each_type=table(reference)
  weight=F1_score*each_type
  weight_f1=sum(weight)/length(reference)
  Macro_precision=sum(precision)/class_num
  analysis=list(result.eva=pre_eva,weight_f1=weight_f1,Macro_precision=Macro_precision)
  return(analysis)
  ana_loc=paste0(rootpath,"/analysis.RData")
  save(analysis,file = ana_loc)
}
#主函数
rootpath="/home/rstudio/kitematic/STAD_latter"
result=main_classification(rootpath)
result_ana=result_analysis(result)
all.result=list(modlist=modlist,dif=dif,min_delta=min_delta_list,spe_edges=spe_edges_loc)
result_loc=paste0(fea_dir,"/result200_dif1.Data")
save(all.result,file=result_loc)
load("/home/rstudio/kitematic/BRCA_Subtypes/result.RData")