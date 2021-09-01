calculatemacro<-function(remat){
  ntype=nrow(remat)
  p=remat$precision
  macro_p=sum(p)/ntype
  r=remat$recall
  macro_r=sum(r)/ntype
  macro_f1=(2* macro_p*macro_r)/(macro_r+macro_p)
  print("macro-p:")
  print(macro_p)
  print("macro-r:")
  print(macro_r)
  print("macro-f1:")
  print(macro_f1)
}
re=read.table("C:/Users/ËïÅåË¶µÄmagic/Desktop/STAD-200.txt",header= T,check.names = F)
calculatemacro(re)