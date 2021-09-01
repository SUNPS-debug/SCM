Enrichment <- function(datafile,path,filename){
  data=datafile
  up=data[data[,2]==1,1]
  all_pathway=all_pathway_canonica_GO_bp_unify.1
  all_symbol=data[,1]
  all_list=unique(all_pathway[,1])
  pathway_enrich_up=NULL
  pvalue=c()
  for(i in 1:length(all_list)){
    gene_pathway=all_pathway[which(all_pathway[,1]==all_list[i]),2]
    up_no=length(intersect(gene_pathway,up))
    hit_no=length(intersect(gene_pathway,all_symbol))
    pval_up=1-phyper(up_no,length(up),length(all_symbol)-length(up),hit_no)
    pathway_enrich_up=rbind(pathway_enrich_up,list(as.character(all_list[i]),pval_up,up_no,hit_no,intersect(gene_pathway,up)))
    pvalue=c(pvalue,pval_up)
  }
  colnames(pathway_enrich_up)=c("pathway","pval_up","up_no","hit_no","gene_list")
  oid=order(as.numeric(pvalue))
  final_pathway_enrich_up=pathway_enrich_up[oid,]
  file1=paste(path,filename,"_up.txt",sep="")
  write.table(final_pathway_enrich_up,file=file1,sep="\t",row.names = FALSE,quote=F)
}

match_gene_list<-function(gene_file,all_symbol){
  up_gene=read.table(gene_file)
  up_gene=up_gene[,1]
  all_symbol=all_symbol[,1]
  allgene_num=length(all_symbol)
  mark=rep(0,allgene_num)
  match_index=match(up_gene,all_symbol)
  mark[match_index]<-1
  data=cbind(all_symbol,mark)
  return(data)
}

typename="class4"
all_symbol=read.table("C:/Users/ËïÅåË¶µÄmagic/Desktop/BRCA/gene_symbol_list_exceptSLC35E2")
gene_file=paste0("C:/Users/ËïÅåË¶µÄmagic/Desktop/BRCA/richment/",typename,".txt")
datafile=match_gene_list(gene_file,all_symbol)
Enrichment(datafile,"C:/Users/ËïÅåË¶µÄmagic/Desktop/BRCA/richment/",typename)