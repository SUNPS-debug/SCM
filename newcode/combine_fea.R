
traindir="D:/WORK/样本多分类/trainset/ORtrainandtest/train_data"
train_outdir="D:/WORK/样本多分类/trainset/train_data"
combine_feaure<-function(indir,outdir){
  ntype=4
  fea_num=100
  combine_feature=list()
  for(i in 1:10){
    train_file=paste0(indir,i,".csv")
    traindata=read.csv(train_file,check.names = F,row.names = NULL)
    data_len=ncol(traindata)
    train_tar=traindata[,data_len]
    part_fea_train=list()
    cut_point=seq(2, data_len, by= fea_num)
    for(j in 1:ntype){
      part_train=traindata[,cut_point[j]:(cut_point[j+1]-1)]
      #all_fea_train=apply(part_train,1,sum)
      part_fea_train[[j]]=apply(part_train,1,sum)
      # out_file=paste0(filedir,"/test_data",i,j,".csv")
      # write.csv(part_fea_train[[j]],out_file)
    }
    combine_train=do.call(cbind,part_fea_train)
    combine_feature[[i]]=cbind(combine_train,train_tar)
    outpath=paste0(outdir,i,".csv")
    write.csv(combine_feature[[i]],file=outpath,row.names = F)
  }
}
indir="D:/WORK/样本多分类/testset7/testdata/test_data"
outdir="D:/WORK/样本多分类/testset7/test_data"
combine_feaure(indir,outdir)
combine_feaure(traindir,train_outdir)