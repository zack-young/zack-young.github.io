library(ggpointdensity)
library(viridis)
library(reshape2)
library(ggrepel)
DF1 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/entropy/CNC_entropy",sep=""),as.is = T, header = F, comment.char = "")
DF1$V4='CNC'
colnames(DF1) <- c("chr1","start1","num1",'type')
DF1$num1 <- sprintf("%.3f",DF1$num1)
DF2 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/entropy/CNL_entropy",sep=""),as.is = T, header = F, comment.char = "")
DF2$V4='CNL'
colnames(DF2) <- c("chr1","start1","num1",'type')
DF2$num1 <- sprintf("%.3f",DF2$num1)
noise <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/noise_CNV.txt",as.is = T, header = F, comment.char = "")
noise$uniq <- paste(noise$V1,noise$V2)
gene_lis <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/gene_all_list_merged.txt",sep=""),as.is = T, header = F, comment.char = "")
colnames(gene_lis) <- c("chr1","start1","type",'num1')
gene_lis$type <- "gene"
gene_lis <- gene_lis[,c(1,2,4,3)]
DF3 <-rbind(DF1,DF2)
DF4 <- dcast(DF3,chr1+start1~type,value.var = 'num1')
DF4$CNC <- as.numeric(DF4$CNC)
DF4$CNL <- as.numeric(DF4$CNL)
DF4$uniq <- paste(DF4$chr1,DF4$start1)
DF3 <-rbind(DF1,DF2,gene_lis)
DF5 <-  dcast(DF3,chr1+start1~type,value.var = 'num1')
DF5 <-DF5[complete.cases(DF5),]
DF5$uniq <- paste(DF5$chr1,DF5$start1)
for (line in noise$uniq) {
  if (length(which(DF4$uniq==line))!=0) {
    DF4 <- DF4[-which(DF4$uniq==line),]
  }
  if (length(which(DF5$uniq==line))!=0) {
    DF5 <- DF5[-which(DF5$uniq==line),]
  }

}

DF5$CNC <- as.numeric(DF5$CNC)
DF5$CNL <- as.numeric(DF5$CNL)

DF6 <- DF5[DF5$CNL-DF5$CNC>2,]
ch_gene <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/ppd1/chosen_gene",as.is = T, header = F, comment.char = "")
DF6 <- DF5[grepl(ch_gene[1,1], DF5$gene),]
for (i in ch_gene$V1[-1]) {
  DF6<-rbind(DF5[grepl(i, DF5$gene),],DF6)
  
}

ggplot(data=DF4,aes(x=CNL,y=CNC)) +
  geom_polygon(data=data.frame(x=c(0,0,2,7,5,3),y=c(2,0,0,5,5,5)),aes(x=x,y=y),fill="red",alpha = 0.2)+
  geom_pointdensity(adjust = 1)+scale_color_viridis()+
  #geom_point(aes(alpha=0.1),colour = "grey")+geom_jitter(width=0.1,height=0.1,colour="grey")+
  #geom_abline(slope=1 ,intercept = -2,)+
  geom_abline(slope=1 ,intercept = 0)+
  #geom_abline(slope=1 ,intercept = 2,)+
  geom_text_repel(data=DF6,aes(x=CNL,y=CNC,label=gene))+
  geom_point(data=DF6,aes(x=CNL,y=CNC,size=1),colour="#F23E30")+ scale_size_continuous(range = 3)+
  theme(legend.position="right")+
  theme_bw()+coord_fixed()+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.line.y = element_line(size=0.5, colour = "black"),
        axis.ticks.y=element_line(size=0.5, colour = "black"),
        axis.text.y=element_text(size=15,color="black"),
        axis.text.x=element_text(size=15,color="black"),
        axis.title.y=element_text(size=16,color="black"),
        axis.ticks.x=element_line(size=0.5, colour = "black"),
        axis.title.x=element_blank(),
        axis.line.x = element_line(size=0.5, colour = "black"),
        legend.title = element_blank())+ scale_x_continuous(limits = c(0,7),expand = c(0,0.05))+
  scale_y_continuous(limits = c(0,5),expand = c(0,0.05))+
  ylab('CNC')+xlab('CNL')
