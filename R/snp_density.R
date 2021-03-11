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
  make_option("--sample", dest = "sample", default = "",
              help = "sample serial number"),
  make_option("--sample_name", dest = "sample_name", default = "",
              help = "sample name"),
  make_option("--sample_number", dest = "sample_number", default = "",
              help = "sample name"),
  make_option(c("-c","--color"), dest = "color", default = "#282658,#32519D,#39BAEC",
              help = "[opt] three colors needed, in order of low,middle,high [default: blue,white,red]"),
  make_option(c("-W","--width"), dest = "figure.width", default = 6,
              help = "[opt] width of figure (inch). [default: 6]"),
  make_option(c("-H","--height"), dest = "figure.height", default = 5,
              help = "[opt] height of figure (inch). [default: 5]")
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
#
figure.width <- arguments$figure.width
figure.height <- arguments$figure.height
#
SAMPLE=arguments$sample
SAMPLE_name=arguments$sample_name
#
SAMPLE_num=arguments$sample_number
#
if(SAMPLE_num == "1"){
  files<-read.table(paste("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/tmp/",SAMPLE,"/combine.",SAMPLE,"_alt_homosnp_density",sep=""),header=F,sep="\t")
  pdf(paste("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/pdf/",SAMPLE_name,"Alt_Homosnp_Density",".pdf",sep=""), height = figure.height, width = figure.width)
  
  hist(log(files[,5],10),breaks = 100,main=paste(SAMPLE_name,"Alt Homosnp Density",sep=""),xlab = "log10(alt homo_snp counts)",prob=T
       #,axes=F
       #,xlim =c(0,2)
       #,ylim = c(0,10)
  )
  dev.off()
  
} else{
  files<-read.table(paste("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/tmp/",SAMPLE,"/combine.",SAMPLE,"_unmatchhomo_snp_density",sep=""),header=F,sep="\t")
  pdf(paste("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/pdf/",SAMPLE_name,"unmatch_homo_snp_density",".pdf",sep=""), height = figure.height, width = figure.width)
  
  hist(log(files[,4],10),breaks = 100,main=paste(SAMPLE_name,"unmatch homo_snp density",sep=""),xlab = "log10(unmatch homo_snp counts)",prob=T
       #,axes=F
       #,xlim =c(0,2)
       ,ylim = c(0,1.5)
  )
  dev.off()
}

#axis(1,seq(from=0, to=4, by=0.5))
#axis(2,seq(from=0, to=2, by=0.5))