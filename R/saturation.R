library(ggplot2)
files = read.table("/data/user/yangzz/mapping/fieldergenomecompare/saturation_test/combine_stat.txt", as.is = T, header = F, comment.char = "")
filesd <- data.frame(homo=c(),hete=c(),type=c())

num=1
for (i in unique(files$V3)) {
  filesd[num,1] <- sum(files[which(files$V3==i),]$V1)
  filesd[num,2] <- sum(files[which(files$V3==i),]$V2)
  filesd[num,3] <- i
  num= num+1
}
colnames(filesd) <- c("homo","hete","type")
filesd$V4 <- c(1,3,5,6,7,8,9,11,13,14)
filesd$homo <- as.numeric(filesd$homo)
filesd$hete <- as.numeric(filesd$hete)
ggplot(data=filesd,aes(x=V4,y=(homo)/	29569032)) +
  geom_point(colour = "black",size=3)+
  #geom_text(data=DF5,aes(x=CNL,y=CNC,label=gene))+
  theme(legend.position="right")+
  theme_bw()+scale_y_continuous(limits = c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","20%","40%","60%","80%","100%"))+
  geom_line(linetype="longdash")+
  scale_x_continuous(limits = c(1,14),breaks = seq(1,14),labels = c("1x","2x","3x","4x","5x","6x","7x","8x","9x","10x","11x","12x","13x","14x"))+
  theme(panel.grid = element_blank(), panel.background = element_blank(),          
        axis.line.y = element_line(size=0.5, colour = "black"),
        axis.ticks.y=element_line(size=0.5, colour = "black"),
        axis.text.y=element_text(size=15,color="black"),
        axis.text.x=element_text(size=15,color="black"),
        axis.title.y=element_text(size=16,color="black"),
        axis.ticks.x=element_line(size=0.5, colour = "black"),
        axis.title.x=element_text(size=15,color="black"),
        axis.line.x = element_line(size=0.5, colour = "black"),legend.position="none")+
  xlab("Sequencing Coverage")+ylab("Percentage")

