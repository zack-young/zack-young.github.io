##draw hap types density
library(reshape2)
library(ggplot2)
library(ggridges)
file1 = read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/CNC_CNL_ALL_hap_stat_density.txt", as.is = T, header = F, comment.char = "")
file1 = read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/test.txt", as.is = T, header = F, comment.char = "")
median(file1$V3, na.rm = FALSE)

file1[which(stringr::str_detect(file1$V1,'A')),1] = 'A'
file1[which(stringr::str_detect(file1$V1,'B')),1]= 'B'
file1[which(stringr::str_detect(file1$V1,'D')),1]= 'D'
file1$V8 = file1$V5/file1$V3
file1$V9 = file1$V6/file1$V3
file1$V10 = file1$V7/file1$V3
#aql_all <-melt(file1,measure.vars = c("V1","V2"),id.vars = c( "V1","V2"))

#file2 = read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/sample_type_stat_CNL", as.is = T, header = F, comment.char = "")

#ggplot(file1) +
  #geom_boxplot(aes(x=c(rep('Introduced Cultivars',nrow(file1)),rep('Chinese Landraces',nrow(file1))),
  #                y=c(file1$V8,file1$V9),fill='#8DD3C7'))+
  #geom_bar(aes(x=V2,fill='#747070'))+
  #geom_histogram(data=file1,aes(x=V3,alpha = 1/10,fill=V1),bins = 80)+
  
  #geom_density(data=file1,aes(x=V3,alpha = 1/10,fill=V1),adjust=2,n =300)+ #,,fill=V1
  #geom_density(data=file1,aes(x=V4,alpha = 1/10,fill=V1),n = 300,adjust=2)+ #fill="#C55948",
  #geom_density(data=file1,aes(x=V8,alpha = 1/10),n = 300,fill="#E8D93A")+ #fill="#C55948",
library(car) 
A_CNC <- file1[which(file1$V1=='A'),]$V3
B_CNC <- file1[which(file1$V1=='B'),]$V3
D_CNC <- file1[which(file1$V1=='D'),]$V3
A_CNL <- file1[which(file1$V1=='A'),]$V4
B_CNL <- file1[which(file1$V1=='B'),]$V4
D_CNL <- file1[which(file1$V1=='D'),]$V4
shapiro.test(B_test[1:4999])
qqnorm(A_test)
qqline(A_test, col = 3)
qqnorm(ttest[21:40,2]);
qqPlot(lm(value~group, data = richness_12), simulate = TRUE, main = 'QQ Plot', labels = FALSE)
abline(0,1,lwd=2)
shapiro <- tapply(A_test, D_test, shapiro.test)
t.test(D_CNC,A_CNC,var.equal=F)



ggplot(file1)+
  geom_density_ridges(aes(x = V4, y = V1, fill=V1),alpha=0.8,scale = 4,na.rm=T,bandwidth=1,quantile_lines=TRUE,
                      quantile_fun=function(x,...)mean(x))+
  #geom_line(x=c(10,10),y=c(0,10))+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.border = element_blank(),          
        axis.line.y = element_line(size=0.5, colour = "black"),
        axis.ticks.y=element_line(size=0.5, colour = "black"),
        axis.text.y=element_text(size=15,color="black"),
        axis.text.x=element_text(size=15,color="black"),
        axis.title.y=element_text(size=16,color="black"),
        axis.ticks.x=element_line(size=0.5, colour = "black"),
        axis.title.x=element_blank(),
        axis.line.x = element_line(size=0.5, colour = "black"),
        legend.title = element_blank())+ylab("Proportion")+ #+scale_fill_discrete(labels=c("CNC","CNL"))
  #scale_y_continuous(limits = c(0,0.2))+
  scale_x_continuous(limits = c(0,50))+ylab("") + xlab("")

##
dt_1 <- as.matrix()
dt_1$V1 <- c(24,23,6,2,1,2,3,1,1,2,1)  
barplot(b,names.arg = c("1","2","3","4","5","7","8","9","11","12","27"),
        ylim = c(0,25),xlab = "Sample Counts",ylab="Hap Counts")
b2=c(18,20,6,2,1,1,1,1,1,1,1)
barplot(b2,names.arg = c("1","2","3","4","5","11","12","14","16","21","49"),
        ylim = c(0,25),xlab = "Sample Counts",ylab="Hap Counts")

for (i in c('chr1A','chr2A','chr3A','chr4A','chr5A','chr6A','chr7A','chr1B','chr2B','chr3B','chr4B','chr5B','chr6B','chr7B','chr1D','chr2D','chr3D','chr4D','chr5D','chr6D','chr7D')){#seq(1,(2*nrow(karyotype)),2)){
  file_tmp <- file1[which(file1$V1==i),]
  p3 <- ggplot(file_tmp)+
    
    geom_line(aes(V2,V3), col="#619CFF",size=0.7)+
    geom_line(aes(V2,V4), col="#D85D4A",size=0.7)+
    #geom_line(aes(start1,num1-num2), col="#6A3D9A",size=0.7)+
    #geom_text(data = gene_df,aes(Position*1000000,5,label=Gene),angle=45,hjust = 0,size=3)+
    #geom_text(data = gene_df,aes(Position*1000000,6.8,label=Gene),angle=45,hjust = 0,size=3)+
    #geom_vline(xintercept = gene_df$Position*1000000,linetype=2)+
    #annotate("rect", xmin = qrc_df$V2*1000000, xmax = qrc_df$V3*1000000, ymin = -Inf, ymax = 4, alpha = 0.5,
    #         fill = "#22BBEE")+
    theme_bw()+
    theme(panel.grid =element_blank(),panel.border = element_blank(),
          axis.line.y = element_line(size=0.5, colour = "black"),
          axis.ticks.y=element_line(size=0.5, colour = "black"),
          axis.text.y=element_text(size=15,color="black"),
          axis.title.y=element_text(size=16,color="black"),
          axis.ticks.x=element_line(size=0.5, colour = "black"),
          axis.title.x=element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.text.x=element_blank(),
          plot.margin = unit(c(1,1,1,0), "cm"))+ #top,right,bottom,left
    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
    #geom_line(aes(density, mean), size=1.2, col="red") +
    #geom_vline(xintercept = 30) +
    scale_x_continuous(expand = c(0.02,0))+ #limits=c(0,860000000), breaks = seq(0,860000000,200000000),
    scale_y_continuous(expand = c(0,0))+ 
    #scale_y_continuous(limits=c(4,7), breaks=c(4,5,6,7),labels = c(4,5,6,7)) +
    ylab("hap") + xlab("")  
  # chr_num=chr_num+1
  tmplist[[chr_num]] <- p3
  chr_num=chr_num+1
  
}
