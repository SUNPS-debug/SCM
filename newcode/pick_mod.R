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
#输入为两个类的模块信息，包括模块名及所对应的基因
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

mod_evaluate<-function(over_list,modlist){
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
    weight=array()
    for(k in 1:(num-1)){
    weight2=which(modorder1==modorder2[k])  
    weight[k]=k+weight2
    }
    names(weight)=modorder2
    minmod_index=which.min(weight)
    dif_mod_name=modorder2[minmod_index]
    print(dif_mod_name)
    dif_mod[[i]]=dif_mod_name
    dif_gene_set[[i]]=modlist[[i]][[dif_mod_name]]
  }
  dif=list(dif_mod=dif_mod,dif_gene_set=dif_gene_set)
  return(dif)
}


##根据score评估mod,返回值是dif_mod和dif_gene_set
mod_eva<-function(score,modlist){
  dif_mod=list()
  dif_gene_set=list()
  for(i in 1:ntype){
    num =nrow(score[[i]]) 
    every_eva=array()
    for(j in 2:num){
      min_index=which.min(score[[i]][j,])  
      every_eva[j]=score[[i]][j,min_index]
      }
    maxmod_index=which.max(every_eva)
    print(maxmod_index)
    dif_mod[[i]]=names(modlist[[i]])[[maxmod_index]]
    dif_gene_set[[i]]=modlist[[i]][[maxmod_index]]
  }
  dif=list(dif_mod=dif_mod,dif_gene_set=dif_gene_set)
  return(dif)
}

del_overlap_gene<-function(dif_gene_set){
  classn=rep(1:ntype)
  overgene=list()
  for(i in 1:ntype){
    overgene_list=list()
    for(j in classn[-i]){
      overgene_list[[j]]=intersect(dif_gene_set[[i]],dif_gene_set[[j]])
    }
    overgene[[i]]=Reduce(union,overgene_list)
    overgene[[i]]=unlist(overgene[[i]])
    dif_gene_set[[i]]=setdiff(dif_gene_set[[i]],overgene[[i]])
  }
  return(dif_gene_set)
}
pick_test_dif<-function(modlist,random){
  dif_mod=list()
  dif_gene_set=list()
  for(i in 1:ntype){
    dif_mod[[i]]=names(modlist[[i]])[[random[i]]]
    dif_gene_set[[i]]=modlist[[i]][[random[i]]]
  }
  dif=list(dif_mod=dif_mod,dif_gene_set=dif_gene_set)
  return(dif)
}
random=c(2,2,2,2,2)
dif=pick_test_dif(modlist,random)

modlist=modtolist(modules)
overlap_list=mod_overlap_all(modlist)
overlap_list_loc=paste0(rootpath,"/overlap_list.RData")
save(overlap_list,file=overlap_list_loc)
dif=mod_evaluate(overlap_list,modlist)
#dif=mod_eva(scores,modlist)
dif_mod=dif$dif_mod
dif_gene_set=dif$dif_gene_set
dif_gene_set=del_overlap_gene(dif_gene_set)
