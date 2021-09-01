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

mod_edges=save_dif_mod(cor_cut_all,dif_mod,dif_gene_set)
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
delta_mat=calculate_delta(mod_edges)
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

min_delta_list=pick_min_delta(delta_mat,dif_gene_set)
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
spe_edges_loc=pick_spe_edges(min_delta_list,100)