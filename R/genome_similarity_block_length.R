files = read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/sample_sim_three_sample/simi_binlen.txt", as.is = T, header = F, comment.char = "")
files[which(files[,4]>=0.5),5]=5#0.5
files[which(files[,4]>=0.25&files[,4]<0.5),5]=4#0.25
files[which(files[,4]>=0.125&files[,4]<0.25),5]=3#0.125
files[which(files[,4]>=0.0625&files[,4]<0.125),5]=2#0.0625
files[which(files[,4]>=0.03125&files[,4]<0.0625),5]=1#0.03125
files[which(files[,4]>=0&files[,4]<0.03125),5]=0
files$V5<- factor(files$V5, levels = c("0","0.03125","0.0625","0.125","0.25","0.5"),labels = c('[0,1/32]','[1/32,1/16]','[1/16,1/8]','[1/8,1/4]','[1/4,1/2]','[1/2,1]'))
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

ggplot(files, aes(x=V5, y=log2(V3))) + geom_violin(width = 0.75,colour = "black",fill="#17A3C0")+
  geom_boxplot(outlier.colour = NA,width = 0.075,colour = "black")+ #geom_smooth(method="lm")
  theme(legend.position="right")+
  labs(x="Similarity", y = "Percentage") + 
  scale_y_continuous(limits=c(0, 9), breaks=c(0,2,4,6,7,log2(400)),labels = c(0,4,16,64,128,400))+
  #scale_y_continuous(breaks =seq(0, 7,1))+
  theme_bw()+
  theme(panel.grid =element_blank(),panel.border = element_blank(),
        axis.line.y = element_line(size=1, colour = "black"),
        axis.ticks.y=element_line(size=1, colour = "black"),
        axis.text.y=element_text(size=15,color="black"),
        axis.text.x=element_text(size=15,color="black"),
        axis.title.y=element_text(size=16,color="black"),
        axis.title.x=element_text(size=16,color="black"),
        axis.ticks.x=element_line(size=1, colour = "black"),
        axis.line.x = element_line(size=1, colour = "black"))

ggplot(files, aes(x=V5, y=log2(V3))) + geom_boxplot()+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median,
               position = position_dodge(width = 1))+theme_bw()+
  theme(panel.grid =element_blank(),panel.border = element_blank(),
        axis.line.y = element_line(size=1, colour = "black"),
        axis.ticks.y=element_line(size=1, colour = "black"),
        axis.text.y=element_text(size=15,color="black"),
        axis.text.x=element_text(size=15,color="black"),
        axis.title.y=element_text(size=16,color="black"),
        axis.title.x=element_text(size=16,color="black"),
        axis.ticks.x=element_line(size=1, colour = "black"),
        axis.line.x = element_line(size=1, colour = "black"))+
  ylab("Hap Length") + xlab("")

ggplot(files, aes(x=V5, y=log2(V3))) + geom_point() +stat_smooth(method = 'lm', formula = y~poly(x, 2), se = TRUE, level = 0.95, color = 'green3')

#+ geom_smooth(method = 'loess', span = 0.75, se = TRUE)
