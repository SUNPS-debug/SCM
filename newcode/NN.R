cvlist=list()
for(cv in 1:ntype){
  datasize=nrow(new_cleandata[[cv]])
  cvlist[[cv]]=CVgroup(10,datasize,cv)
}
model=list()
predictions=list()
conMats=list()
for(pretime in 1:10){
  print(paste0("交叉验证第",pretime,"轮"))
  train_test=trainandtest(new_cleandata,cvlist,pretime)
  traindata=train_test$traindata
  testdata=train_test$testdata
  traintarget=list()
  testtarget=list()
  for(nt in 1:ntype){
    nsam=nrow(traindata[[nt]])
    traintarget[[nt]]=rep(nt,nsam)
    traindata[[nt]]=traindata[[nt]][,features_gene2]
    testdata[[nt]]=testdata[[nt]][,features_gene2]
    testtarget[[nt]]=rep(nt,dim(testdata[[nt]])[1])
  }
  traintar=unlist(traintarget)
  trainda=do.call(rbind,traindata)
  testtar=unlist(testtarget)
  testda=do.call(rbind,testdata)
  train_Targets <- decodeClassLabels(traintar)
  test_Targets <- decodeClassLabels(testtar)
  df=list(inputsTrain=trainda,targetsTrain=train_Targets,inputsTest=testda,targetsTest=test_Targets)
  norm_data<- normTrainingAndTestSet(df)
  print("training")
  model[[pretime]] <- mlp(norm_data$inputsTrain, norm_data$targetsTrain, size=60, learnFunc="BackpropBatch", learnFuncParams=c(10, 0.1), maxit=2000, inputsTest=norm_data$inputsTest, targetsTest=norm_data$targetsTest)
  print("testing")
  predictions[[pretime]] = predict(model[[pretime]],norm_data$inputsTest)
  #生成混淆矩阵，观察预测精度
  print("predict result")
  conMats[[pretime]]=confusionMatrix(norm_data$targetsTest,predictions[[pretime]])
}
pres=do.call(rbind,predictions)