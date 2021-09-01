library(ggplot2)
setwd("C:/Users/孙佩硕的magic/Desktop/BRCA")
compare_result=read.csv("macro-f1.csv")
methods=compare_result$methods
macro.precision=compare_result$values
number_of_features=compare_result$features_number
png(file = "macro-f1.png",width=800,height=300,res=100)
col=c("darkolivegreen3","cornflowerblue","gold2", "sandybrown","pink2")
ggplot(compare_result, aes(x=factor(number_of_features), y=macro.precision,group=methods, 
                           colour=methods,shape=methods))+ geom_line(size=1)+geom_point(size=4)+
  scale_colour_manual(values=col)+
  xlab("Number of features")+#横坐标名称
  ylab("Macro-avg F1-score")+#纵坐标名称 +#去掉背景灰色
  theme_bw()
  # theme(panel.grid.major=element_line(colour=NA),
  #       panel.background = element_rect(fill = "transparent",colour = NA),
  #       plot.background = element_rect(fill = "transparent",colour = NA),
  #       panel.grid.minor =element_blank())#以上theme中代码用于去除网格线且保留坐标轴边框
  # #       #text = element_text(family = "STXihei"),#设置中文字体的显示
  # #       #legend.position = c(.075,.915),#更改图例的位置，放至图内部的左上角
  # #       #legend.box.background = element_rect(color="black"))#为图例田间边框线

dev.off()


