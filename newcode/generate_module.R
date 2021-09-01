library(WGCNA)
library(reshape2)
library(stringr)
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
datacleaning <- function(data) {
  mean_gene=apply(data,1,mean)
  remain_gene=which(mean_gene>=10)
  dataExpr=data[remain_gene,]
  dataExpr = log(dataExpr + 1)
  m.mad <- apply(dataExpr, 1, mad)
  dataExprVar <- dataExpr[which(m.mad >
                                  max(quantile(m.mad, probs = seq(0, 1, 0.1))[2], 0.01)), ]
  
  dataExpr <- as.data.frame(t(dataExprVar))
  print(dim(dataExpr))
  gsg = goodSamplesGenes(dataExpr, verbose = 3)
  if (!gsg$allOK) {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes) > 0)
      printFlush(paste("Removing genes:",
                       paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")))
    
    if (sum(!gsg$goodSamples) > 0)
      printFlush(paste("Removing samples:",
                       paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")))
    
    # Remove the offending genes and samples from the data:
    dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
    print(dim(dataExpr))
  }
  return(dataExpr)
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
mean_corandp<-function(sub_sampling_data){
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
generate_module2<-function(sub_meancor,power,datExpr,classn){
  adj=(sub_meancor)^power
  dissadj = 1-abs(adj)
  geneTree = hclust(as.dist(dissadj), method = "average")
  minModuleSize = 10
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissadj,
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

generate_module<-function(datExpr,classn){
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(datExpr, RsquaredCut = 0.9, 
                          powerVector=powers, 
                          networkType="unsigned", verbose=5)
  power = sft$powerEstimatepower
  adjacency=cor(datExpr,method = "spearman")
  adj=(adjacency)^power
  TOM = TOMsimilarity(abs(adj),TOMType = "unsigned")
  rownames(TOM)=rownames(adj)
  colnames(TOM)=colnames(adj)
  dissTOM = 1-abs(TOM)
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = 10
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = FALSE, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  names(dynamicMods)=rownames(TOM)
  #dynamicColors = labels2colors(dynamicMods)
  ### 合并表达图谱相似的模块
  #计算模块的主成分
  MEList = moduleEigengenes(datExpr, colors = dynamicMods,
                            nPC = 2, softPower = power,)
  MEs = MEList$eigengenes
  varExp=MEList$varExplained
  npcpath=paste0(rootpath,"/",cancertype[classn],"/PCA.csv")
  write.csv(varExp,npcpath)
  
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  MEDissThres = 0.1
  
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


rootpath="/home/rstudio/kitematic/LUSC_Stages_test"
options(stringsAsFactors = FALSE)
# 打开多线程
enableWGCNAThreads()
txtpath=paste0(rootpath,"/LUSC_Stages")
ORpath=paste0(rootpath,"/data")
txtname<<-list.files(txtpath)
ntype<<-length(txtname)
tmp<-strsplit(txtname,split=".",fixed=TRUE)
cancertype<<-unlist(lapply(tmp,head,1)) 
genename_loc<<-paste0(rootpath,"/gene_symbol_list_exceptSLC35E2")
genename<<-read.table(genename_loc)
create_dir(rootpath)
txt_to_csv(txtpath)
data_num=data_num(ORpath)
ORdata=data_num$data
num_list=data_num$numlist
cleandata=list()
for(i in 1:ntype){
  cleandata[[i]]=datacleaning(ORdata[[i]])
}
common_genelist=common_gene(cleandata)
new_common_genelist=common_genelist[11:11629]
genenum<<-length(new_common_genelist)
new_cleandata=filter_cleandata(cleandata,new_common_genelist)
sampling_dir=paste0(rootpath,"/samplingdata/",cancertype)
for(i in 1:ntype){
  if(dir.exists(sampling_dir[i])==F)
    dir.create(sampling_dir[i])
}
sampling_data_all=sampling_all(new_cleandata)
# options(stringsAsFactors = FALSE)
# # 打开多线程
# enableWGCNAThreads()
meancor=list()
meanp=list()
modules=list()
for(j in 1:ntype){
  print(paste0("type:",j))
  # corandp=mean_corandp(sampling_data_all[[j]])
  # meancor[[j]]=corandp$mean_cor
  # meanp[[j]]=corandp$mean_p
  modules[[j]]=generate_module(new_cleandata[[j]],j)
  }

meancor_loc=paste0(rootpath,"/meancor.RData")
meanp_loc=paste0(rootpath,"/meanp.RData")
save(meancor,file=meancor_loc)
save(meanp,file=meanp_loc)
# de=list()
# for(i in 1:ntype){
#   data=new_cleandata[[i]]
#   cor_mat=cor(data,method = "spearman")
#   de[[i]]=cor_mat-meancor[[i]]
#   max_de=max(abs(de[[i]]))
#   print(max_de)
# }
