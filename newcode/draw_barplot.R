library(ggplot2)
classification_path=("C:/Users/ËïÅåË¶µÄmagic/Desktop/STAD_latter/STAD-200feature.csv")
cla_result=read.csv(classification_path,check.names = F)
cols=c("darkolivegreen3","sandybrown","cornflowerblue","gold2","pink2","lightblue2")
p <- ggplot(cla_result, aes(fill=Method, y=Values, x=Method)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = cols)+
  facet_wrap(~Evaluation,strip.position = "bottom") +
  xlab("")+
  theme_classic()+
  ggtitle("STAD-200features")+
  # coord_flip() +
  scale_x_discrete(labels=levels(Evaluation))+ #ÀëÉ¢µÄ±êÇ©
  coord_cartesian(ylim = c(0.5,1))+
  scale_y_continuous(breaks=seq(0.5,1,0.05))+
  theme(axis.text.x=element_text(angle = -90, hjust = 0),
        strip.text = element_text(size = 13),
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size=13),
        plot.title=element_text(size=20)) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())
 
  
ggsave(p,filename = "C:/Users/ËïÅåË¶µÄmagic/Desktop/STAD_latter/200fea_plot.png", width = 20, height = 15, units = "cm")