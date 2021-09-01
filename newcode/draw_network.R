library(igraph)
change_index<-function(spe_edges,dif_gene_set){
  gene_num=length(dif_gene_set)
  gene_index=rep(1:gene_num)
  node1=spe_edges[,1]
  node2=spe_edges[,2]
  for(i in 1:gene_num){
    node1[node1==i]<-dif_gene_set[i]
    node2[node2==i]<-dif_gene_set[i]
  }
  edges_gene=cbind(node1,node2)
  return(edges_gene)
}

draw_network<-function(subtype,fea_num,index,col){
  #inputpath=paste0("D:/WORK/样本多分类/BRCAdif1_",fea_num,"/result",fea_num,".RData")
  inputpath=paste0("D:/WORK/样本多分类/STAD_latter/",fea_num,"feature/result",fea_num,".RData")
  load(inputpath)
  dif_gene_set=all.result$dif$dif_gene_set[[index]]
  spe_edges_index=all.result$spe_edges[[index]]
  edges=change_index(spe_edges_index,dif_gene_set)
  g<-graph_from_data_frame(edges,directed=F)
  draw_path=paste0("C:/Users/孙佩硕的magic/Desktop/STAD_latter/connection_pics/",subtype,"_",fea_num,".pdf")
  pdf(draw_path)
  plot.igraph(g,layout=layout.auto,
              vertex.label.dist=1,
              main=subtype,
              edge.width=1,
              vertex.shape="circle",
              edge.color="gray57",
              vertex.label.cex=0.3,
              vertex.frame.color=col,
              vertex.size=6,edge.arrow.size=0,vertex.color=col) 
  dev.off()
}
cols=c("chartreuse4","gold2","deepskyblue3","lightsalmon","lightpink2")
fea_nums=c(20,50,100,200)
#BRCA 
#subtypes=c("Normal","ER+","HER2+","TNBC")
#STAD
subtypes=c("Normal","CIN","EBV","MSI","GS")
for(i in 1:length(subtypes)){
  for(j in 1:length(fea_nums)){
    draw_network(subtypes[i],fea_nums[j],i,cols[i])
  }
}


