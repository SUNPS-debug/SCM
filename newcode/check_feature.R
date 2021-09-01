# gene_feature_path="C:/Users/孙佩硕的magic/Desktop/STAD/f_classif/feature200.csv"
# gene_feature=read.csv(gene_feature_path,header = F)
# gene_feature=unlist(gene_feature)
fea_num=20
result_path=paste0("D:/WORK/样本多分类/BRCAdif1_",fea_num,"/result",fea_num,".RData")
load(result_path)
ntype=4
dif_gene_set=all.result$dif$dif_gene_set
spe_edges_loc=all.result$spe_edges
pick_key_gene<-function(dif_gene_set,spe_edges_loc){
  features_gene=list()
  for(i in 1:ntype){
    unique_index=unique(as.numeric(spe_edges_loc[[i]])) 
    features_gene[[i]]=dif_gene_set[[i]][unique_index]
  }
  return(features_gene)
}
dif_genes=pick_key_gene(dif_gene_set,spe_edges_loc)
# for(i in 1:ntype){
#   co_features=dif_genes[[i]]
#   outpath=paste0("C:/Users/孙佩硕的magic/Desktop/STAD_latter/class",i,".txt")
#   write.table(co_features,outpath,quote = F,row.names = F,col.names = F)
# }
dif_genes=unlist(dif_genes)
co_gene_feature=unique(dif_genes)
outpath=paste0("C:/Users/孙佩硕的magic/Desktop/BRCA/co-ex/earlystop/feature",fea_num,".csv")
write.table(co_gene_feature,outpath,row.names = F,col.names = F)