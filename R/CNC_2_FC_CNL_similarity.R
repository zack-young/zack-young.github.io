library(dplyr)
library(reshape2)
files = read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/sample_sim_three_sample/CNC_sample_CNL_FC_sim_count.txt", as.is = T, header = F, comment.char = "")
files$V1 <- files$V1/14075
files$V3 <- files$V3/14075
files <- files[which(files[,1]>=0.25),]

files = read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/CNCsample_2_CNL_F.txt", as.is = T, header = F, comment.char = "")
colnames(files) <- c('sample','line_CNL','line_FC','line_CNC','line_FC_CNL','line_oth')
files$OTH <- rowSums(files[, c(4,5,6)])
files$OTH <- files$line_oth - files$line_CNL -files$line_FC - files$line_CNC-files$line_FC_CNL

#files = read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/CNL_FL.txt", as.is = T, header = F, comment.char = "")
#colnames(files) <- c('sample','value')
#aql_all <-files

dt <- files[order(files$line_FC),]
#aql_all <- melt(dt,measure.vars = c('line_CNL','line_FC','line_FC_CNL','line_CNC','OTH'),id.vars = c( "sample"))
aql_all <- melt(dt,measure.vars = c('line_CNL','line_FC'),id.vars = c( "sample"))
#aql_all <- melt(dt,measure.vars = c("V1","V3"),id.vars = c( "V5"))
aql_all$sample <- factor(aql_all$sample,levels = dt$sample)
aql_all$variable <- factor(aql_all$variable,levels = c("line_CNL","line_FC"))
ggplot(aql_all)+
  geom_bar(aes(x = sample, y = value/14075,alpha = 1/10,fill=variable),stat = "identity",position ='stack')+ 
  theme_bw()+ #scale_fill_manual(values = col_lis)+
  theme(panel.grid = element_blank(), panel.border = element_blank(),          
        axis.line.y = element_line(size=0.5, colour = "black"),
        axis.ticks.y=element_line(size=0.5, colour = "black"),
        axis.text.y=element_text(size=15,color="black"),
        axis.text.x=element_text(size=16,color="black",angle=45,hjust = 1),
        axis.title.y=element_text(size=16,color="black"),
        axis.ticks.x=element_line(size=0.5, colour = "black"),
        axis.title.x=element_text(size=16,color="black"),
        axis.line.x = element_line(size=0.5, colour = "black")
        )+ylab("Proportion")
  #scale_y_continuous(limits = c(0,0.7) ,expand = c(0,0))+
  #scale_y_continuous(limits = c(0,1) ,expand = c(0,0))+
   

ggplot(aql_all,aes(x=variable,y=value/14075)) +
  geom_violin(width = 0.75,colour = "black",fill="#17A3C0",trim = T)+
  geom_boxplot(outlier.colour = NA,width = 0.075,colour = "black")+ 
  theme(legend.position="right")+
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
                           legend.title = element_blank())+
  scale_x_discrete(labels=c("Chinese Landrace","exotic germplasm"))+
  scale_y_continuous(limits = c(0,0.7) ,expand = c(0,0))+
  ylab('Proportion')

