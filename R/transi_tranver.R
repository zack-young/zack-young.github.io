#!/usr/bin/env Rscript
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
if(!is.installed("optparse")){
  warning("Detect package \"optparse\" is not installed in your R enviroment.")
  warning("Trying to install the \"optparse\" package.")
  warning("If failed, please try to install it mannually.")
  install.packages("optparse")
}
if(!is.installed("gplots")){
  warning("Detect package \"gplots\" is not installed in your R enviroment.")
  warning("Trying to install the \"gplots\" package.")
  warning("If failed, please try to install it mannually.")
  install.packages("gplots")
}

## libraries
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
# Arguments
option_list <- list(
  make_option(c("-s", "--sample"), dest = "sample", default = "",
              help = "sample name"),
  make_option(c("-c","--color"), dest = "color", default = "#282658,#32519D,#39BAEC",
              help = "[opt] three colors needed, in order of low,middle,high [default: blue,white,red]"),
  make_option(c("-W","--width"), dest = "figure.width", default = 5,
              help = "[opt] width of figure (inch). [default: 5]"),
  make_option(c("-H","--height"), dest = "figure.height", default = 4,
              help = "[opt] height of figure (inch). [default: 4]")
)

parser <- OptionParser(usage = "transi_transver [options]",
                       option_list=option_list, description = "      (aka mCBinHeatmap)\
                       Description: Plot methylation dynamics of target region for multiple samples [heatmap]\
                       Contact:     Zhu, Ping; pingzhu.work@gmail.com\
                       Last update: 2017-09-16\
                       Example: \
                       mCBinHeatmap.R -i input -m white -o chr1.xxx-xxx.pdf \
                       -Input File Format: \
                       1st line is the header.\
                       Each column contains methylation measurements of a sample. \
                       Example: \
                       Sample1  Sample2 ...  \
                       0.1      0.1     ...  \
                       0.1      0.1     ...  \
                       "
)
## check arguments
arguments <- parse_args(parser)
sample = arguments$sample
#
figure.width <- arguments$figure.width
figure.height <- arguments$figure.height
#
color = arguments$color
if(color == ""){ # default, "FragRegView.Date"
  color = strsplit(c("#282658","#32519D","#39BAEC"), split = ",", fixed = T)[[1]]
}else{
  color <- strsplit(arguments$color, split = ",", fixed = T)[[1]]
}

#
different.color <- color[2]
similar.color <- color[3]
#
dataa<-read.table(paste("/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/",sample,"/combine_total_diff_transi_ver",sep=""),header=F,sep="\t",as.is = F)
datab<-read.table(paste("/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/",sample,"/combine_diff_diff_transi_ver",sep=""),header=F,sep="\t",as.is = F)

t_ac_gt_frame <- dataa[which(dataa[,1]=="AC"|dataa[,1]=="GT"|dataa[,1] == "CA"|dataa[,1] == "TG"),]
t_ac_gt_sum <- sum(t_ac_gt_frame[,2])

d_ac_gt_frame <- datab[which(datab[,1]=="AC"|datab[,1]=="GT"|datab[,1] == "CA"|datab[,1] == "TG"),]
d_ac_gt_sum <- sum(d_ac_gt_frame[,2])

t_ct_ga_frame <- dataa[which(dataa[,1]=="CT"|dataa[,1]=="TC"|dataa[,1] == "GA"|dataa[,1] == "AG"),]
t_ct_ga_sum <- sum(t_ct_ga_frame[,2])

d_ct_ga_frame <- datab[which(datab[,1]=="CT"|datab[,1]=="TC"|datab[,1] == "GA"|datab[,1] == "AG"),]
d_ct_ga_sum <- sum(d_ct_ga_frame[,2])

t_gc_frame <- dataa[which(dataa[,1]=="GC"|dataa[,1]=="CG"),]
t_gc_sum <- sum(t_gc_frame[,2])

d_gc_frame <- datab[which(datab[,1]=="GC"|dataa[,1]=="CG"),]
d_gc_sum <- sum(d_gc_frame[,2])

t_at_frame <- dataa[which(dataa[,1]=="AT"|dataa[,1]=="TA"),]
t_at_sum <- sum(t_at_frame[,2])

d_at_frame <- datab[which(datab[,1]=="AT"|dataa[,1]=="TA"),]
d_at_sum <- sum(d_at_frame[,2])

dt_4 <- data.frame(A=c(t_ac_gt_sum,d_ac_gt_sum),B=rep("A/C+G/T",times=2),C=c("Total","Different"))
dt_5 <- data.frame(A=c(t_ct_ga_sum,d_ct_ga_sum),B=rep("C/T+G/A",times=2),C=c("Total","Different"))
dt_1 <- data.frame(A=c(t_gc_sum,d_gc_sum),B=rep("G/C",times=2),C=c("Total","Different"))
dt_3 <- data.frame(A=c(t_at_sum,d_at_sum),B=rep("A/T",times=2),C=c("Total","Different"))
dt_2 <- rbind(dt_1,dt_4,dt_3,dt_5)

options(scipen=200)
dt_2[,2] = factor(dt_2[,2], levels=c("C/T+G/A","A/C+G/T","G/C","A/T"))
pdf(paste("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/pdf/",sample,"_snp_transi_transver",".pdf",sep=""), height = figure.height <- arguments$figure.height, width = figure.width)
ggplot(dt_2,aes(x=B,y=A,fill=C))+
  geom_bar(stat="identity",position = "dodge")+
  scale_fill_brewer(palette = "Pastel1")+theme_classic() + 
  #scale_y_continuous(breaks=seq(0, 3000, 1000),limits = c(0,3000))+#,labels = c(0,1000,2000,3000))+ 
  labs(title=paste(sample,"SNP transition transversion Density",sep=" "), x="SNP Type", y="Density")+scale_fill_manual(values=c(different.color,similar.color))+
  theme(plot.title = element_text(hjust = 0.5)) ##!!!!!change title position

invisible(dev.off())








