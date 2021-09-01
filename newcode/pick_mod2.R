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
cut_byvalue<-function(matcor,geneset,cut_value){
  remain_edge_loc=list()
  mod_size=length(geneset)
  subcor=matcor[geneset,geneset]
  subcor[abs(subcor)<cut]<-0
  remain_edge_loc=edge_loc(subcor)
  return(remain_edge_loc)
}

#返回每个模块的最大边的位置
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
#查看每个模块中比较大的边在不同期中的mean
evaluate_mod<-function(mat_cut_all,modlist,max_edges_all){
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
pick_dif_mod<-function(mean_all,modlist){
  dif_mod=list()
  dif_gene_set=list()
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
      dif_index=which.max(deltamean)
      dif_mod[[i]]=mod_names[dif_index]
      dif_gene_set[[i]]=modlist[[i]][[dif_index]]
  }
  dif=list(dif_mod=dif_mod,dif_gene_set=dif_gene_set)
  return(dif)
}
cor_cut_all=select_every_mod(meancor,modlist)
max_edges_list=mod_max_edges(cor_cut_all,modlist,0.4)

all_edge_loc=mod_max_edges_value(cor_cut_all,modlist,0)
#tlist=evaluate_mod(cor_cut_all,modlist,max_edges_list)
meanlist=evaluate_mod(cor_cut_all,modlist,max_edges_list)
dif2=pick_dif_mod(meanlist,modlist)
dif_mod=dif2$dif_mod 
dif_gene_set=dif2$dif_gene_set
