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
  make_option(c("-W","--width"), dest = "figure.width", default = 7,
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
  color = strsplit("#282658,#32519D,#39BAEC", split = ",", fixed = T)[[1]]
  }else{
    color <- strsplit(arguments$color, split = ",", fixed = T)[[1]]
  }

#
different.color <- color[2]
similar.color <- color[3]
datah<-read.table(paste("/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/",sample,"/combine_total_indel_len_count",sep=""),header=F,sep="\t",as.is = F)
datam<-read.table(paste("/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/",sample,"/combine_diff_indel_len_count",sep=""),header=F,sep="\t",as.is = F)

#datah$V3 = datah$V3/11557
datah$V3 = datah$V3
hp10 <-sum(datah[which(datah[,1] >=10),][,2])
hn10 <- sum(datah[which(datah[,1]<=-10),][,2])

#datam$V3 = datam$V3/1526
datam$V3 = datam$V3
mp10 <-sum(datam[which(datam[,1] >=10),][,2])
mn10 <- sum(datam[which(datam[,1]<=-10),][,2])

lowdatam <-cbind(datam[which(datam[,1] <10 & datam[,1] >-10),],c("Different"))
names(lowdatam)<-c("A","B","C")
lowdatah <-cbind(datah[which(datah[,1] <10 & datah[,1] >-10),],c("Total"))
names(lowdatah)<-c("A","B","C")
usedata <- rbind(lowdatam,lowdatah,c(">=10",hp10,"Total"),c("<=-10",hn10,"Total"),c(">=10",mp10,"Different"),c("<=-10",mn10,"Different"))

usedata$A = factor(usedata$A, levels=c("<=-10","-9","-8","-7","-6","-5","-4","-3","-2","-1","0",
                                         "1","2","3","4","5","6","7","8","9",">=10"))
usedata$B = as.numeric(usedata$B)

pdf(paste("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/pdf/",sample,"_InDeL_length",".pdf",sep=""), height = figure.height <- arguments$figure.height, width = figure.width)

ggplot(usedata,aes(x=A,y=B,fill=C))+
  geom_bar(stat="identity",position = "dodge")+
  scale_fill_brewer(palette = "Pastel1")+scale_fill_manual(values=c(different.color,similar.color))+theme_classic() + 
  #scale_y_continuous(breaks=seq(0, 95,5),limits = c(0,95))+#,labels = c(0,200000,400000,600000,800000,1000000,1200000))+ 
  labs(title=paste(sample,"InDel Length Density",sep=" "), x="INDEL Length", y="Density")+
  theme(plot.title = element_text(hjust = 0.5))
invisible(dev.off())
