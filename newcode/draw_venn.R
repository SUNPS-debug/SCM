library(VennDiagram)
setwd("C:/Users/ËïÅåË¶µÄmagic/Desktop/STAD_latter")
combine_fea_sets<-function(fea_select_list,fea_num){
  fea_select_num=length(fea_select_list)
  feature_list=list()
  for(i in 1:fea_select_num){
    print(fea_select_list[i])
    feature_path=paste0(fea_select_list[i],
                        "/feature",fea_num,".csv")
    feature_file=read.csv(feature_path,header = 0)
    feature_list[[i]]=feature_file[,1]
  }
  names(feature_list)=fea_select_list
  return(feature_list)
}
feanum=100
filepath=paste0("Venn plot/venn",feanum,".png")
fea_select_list=c("HSIC-Lasso","Mutual_info","ANOVA","Chi-square test","SCM")
feature_list=combine_fea_sets(fea_select_list,feanum)
venn.plot <- venn.diagram(
  feature_list,
  filename = filepath,
  main="STAD",
  main.cex=2,
  lty = "dotted",
  lwd = 2,
  col = "black",  #"transparent",
  fill = c("LightPink", "lightgoldenrod1", "lightblue2", "tan2", "darkolivegreen3"),
  alpha = 0.60,
  #cat.col="black"
  cat.col = c("black", "black", "black", "black", "black"),
  cat.cex = 1.3,
  cat.fontface = "bold",
  margin = 0.2,
  cex = 1
)



