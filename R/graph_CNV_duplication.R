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
  make_option("--sample1", dest = "sample1", default = "",
              help = "sample1 serial number"),
  make_option("--sample1_name", dest = "sample1_name", default = "",
              help = "sample1 name"),
  make_option("--sample2", dest = "sample2", default = "",
              help = "sample2 serial number"),
  make_option("--sample2_name", dest = "sample2_name", default = "",
              help = "sample2 name"),
  make_option(c("-c","--color"), dest = "color", default = "#282658,#32519D,#39BAEC",
              help = "[opt] three colors needed, in order of low,middle,high [default: blue,white,red]"),
  make_option(c("-W","--width"), dest = "figure.width", default = 16,
              help = "[opt] width of figure (inch). [default: 5]"),
  make_option(c("-H","--height"), dest = "figure.height", default = 9,
              help = "[opt] height of figure (inch). [default: 4]"),
  make_option(c("-s","--suffix"), dest = "suffix", default = "",
              help = "[opt] suffix of pdf. [default:]")
  
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
color = arguments$color
if(color == ""){ # default, "FragRegView.Date"
  color = strsplit("#282658,#32519D,#39BAEC", split = ",", fixed = T)[[1]]
  }else{
    color <- strsplit(arguments$color, split = ",", fixed = T)[[1]]
  }

#
SAMPLE1=arguments$sample1
SAMPLE1_name=arguments$sample1_name
SAMPLE2=arguments$sample2
SAMPLE2_name=arguments$sample2_name
options(scipen=200)
files = read.table(paste("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plotfile/plotfile_",SAMPLE1,"_",SAMPLE2,'_duplication',sep=""), as.is = T, header = F, comment.char = "")

pdf(paste("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/pdf/",SAMPLE1_name,"_",SAMPLE2_name,arguments$suffix,"_duplicationCNV.pdf",sep=""), height = figure.height, width = figure.width)
plotChrom <- function(xleft, ybottom, height, len, centro, binsize){
  
  #len <- len/binsize
  #centro <- centro/binsize
  
  # r vertical, r horizontal
  rv <- height/2
  rh <- len/70
  rs <- seq(0,pi,len=100)
  
  # left semi-circle
  lx <- c(xleft, xleft + rh - rh*sin(rs), xleft, xleft)
  ly <- c(ybottom, ybottom + rv - rv*cos(rs), ybottom + height, ybottom)
  polygon(lx, ly, border = "white", col = "white")
  
  # right semi-circle
  rx <- c(xleft + len,xleft + len, xleft + len - rh + rh*sin(rs), xleft + len)
  ry <- c(ybottom, ybottom + height, ybottom + rv +rv*cos(rs), ybottom)
  polygon(rx, ry, border = "white", col = "white")
  
  # top tri
  tx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  ty <- c(ybottom + height, ybottom + rv, ybottom + height)
  polygon(tx, ty, border = "white", col = "white")
  
  # bottom tri
  bx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  by <- c(ybottom, ybottom + rv, ybottom)
  polygon(bx, by, border = "white", col = "white")
  
  cbx <- c(xleft + rh - rh*sin(rs), centro - 0.5*rh, centro, centro + 0.5*rh, xleft + len - rh + rh*sin(rs), centro + 0.5*rh, centro, centro - 0.5*rh)
  cby <- c(ybottom + rv - rv*cos(rs), ybottom + height, ybottom + rv, ybottom + height, ybottom + rv + rv*cos(rs), ybottom, ybottom + rv, ybottom)
  polygon(cbx, cby, border = "grey", lwd = 3)
}
plot(x=0, type="n", bty="n", yaxt="n",xaxt="n",
     xlab="", ylab="", 
     xlim=c(-200, 1100), ylim=c(-2, nrow(files)+1),
     xaxs="i", yaxs="i", 
     #main=paste("S",j,sep=""),
     main=paste(SAMPLE1_name,SAMPLE2_name,"snp level and CNV Distribution",sep=" "),
     cex.main = 2.5)
axis(1,at=c(0,100,200,300,400,500,600,700,800),labels = c("0Mbp","","200Mbp","","400Mbp","","600Mbp","","800Mbp"))
#text(x=c(0,200,400,600,800), y=rep(0,5),
#     c("0Mbp","200Mbp","400Mbp","600Mbp","800Mbp"), cex = 1,pos = 1)
#color_pad <- c(colorpanel(low="yellow", high="red", n=3),"green","black","white")
deepblue <- "#282658"
midblue <- "#32519D"
lightblue <- "#39BAEC"
lightred <- "#F6E5D3"
red<-'#F784B6'
purple <- "#A21B5E"
lightpurple <- "#9C7AF0"
lightyellow <- "#FFFEDF"
yellow <- "#F9B731"

#color_pad <- c(deepblue,midblue,lightblue,red,yellow,purple)
color_pad <- c('white','white','white',red,yellow,lightblue)
#color_pad <- c(deepblue,midblue,lightblue,'white','white','white')
num = 22
centromere <- c(215,239,173,337,346,268,300,346,241,300,317,184,254,206,189,286,327,211,357,287,340)
wholen <- c(594102056,
            689851870,
            495453186,
            780798557,
            801256715,
            651852609,
            750843639,
            830829764,
            615552423,
            744588157,
            673617499,
            509857067,
            709773743,
            713149757,
            566080677,
            618079260,
            720988478,
            473592718,
            736706236,
            750620385,
            638686055
)
# rect(xleft = 900, xright = 930, ybottom = 20, ytop = 22,
#      col = color_pad[3], border = color_pad[3]);text(x=c(940), y=c(21),c("Similar"),cex = 1,pos = 4)
# rect(xleft = 900, xright = 930, ybottom = 22, ytop = 24,
#      col = color_pad[1], border = color_pad[1]);text(x=c(940), y=c(23),c("Different"),cex = 1,pos = 4)
# rect(xleft = 900, xright = 930, ybottom = 22, ytop = 24,
#      col = color_pad[3], border = color_pad[3]);text(x=c(940), y=c(23),c("LOW_Diff"),cex = 1,pos = 4)
# rect(xleft = 900, xright = 930, ybottom = 24, ytop = 26,
#      col = color_pad[2], border = color_pad[2]);text(x=c(940), y=c(25),c("MID_Diff"),cex = 1,pos = 4)
# rect(xleft = 900, xright = 930, ybottom = 26, ytop = 28,
#      col = color_pad[1], border = color_pad[1]);text(x=c(940), y=c(27),c("HIGH_Diff"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 28, ytop = 30,
     col = color_pad[4], border = color_pad[4]);text(x=c(940), y=c(29), paste(SAMPLE1,"duplication_CNV",sep=""),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 30, ytop = 32,
     col = color_pad[6], border = color_pad[6]);text(x=c(940), y=c(31),c("Both_duplication_CNV"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 32, ytop = 34,
     col = color_pad[5], border = color_pad[5]);text(x=c(940), y=c(33),paste(SAMPLE2,"duplication_CNV",sep=""),cex = 1,pos = 4)
for(j in 1:nrow(files)){
  if(!is.na(files[j,1])){
    text(x=-100, y=j+0.5, files[j,2], cex = 1)
    if(file.info(files[j,1])$size!=0){
      data1 <- read.table(files[j,1], as.is = T, header = F, comment.char = "")
      #chro <- factor(data1[,5])
      #              if(!is.na(files[j,2])){
      #                chro <- factor(data1[,7])
      #              } else{
      #                chro <- factor(data1[,5])
      #              }
      chro <- factor(data1[,4])
      # chro1 <- factor(chro, levels = c('CNV', 'DP', 'hete', 'homo'), 
      #                 labels = c('1', '2', '3', '4'))
      chro1 <- factor(chro, 
                      #levels = c("level1","level2","level3","unmatchCNV","DP","ownCNV"), 
                      #
                      levels = c("low_com","mid_com","high_com",
                                 paste(SAMPLE1,"duplication_CNV",sep=""),paste(SAMPLE2,"duplication_CNV",sep=""),"duplication_both_CNV"),
                                 #paste(SAMPLE1,"duplication_CNV",sep=""),paste(SAMPLE2,"duplication_CNV",sep=""),"duplication_both_CNV"),
                      labels = c('3','1','1','4',"5","6")
                      #labels = c('3','2','1','4',"5","6")
                      #levels = c("jm22","both","lx99"),
                      #labels = c('4',"5","6")
                      )
      #chro1 <- factor(chro, levels = c('level1', 'level2', 'level3', 'level4',"level5","level6","level7","level8","level9","level10"), 
                      #labels = c('10', '9', '8', '7','6', '5', '4', '3','2', '1')) 
      
      num <- num - 1 
      #assign(paste("sample",num,sep=""),nrow(data1))
      for(i in 1:nrow(data1)){
        rect(xleft = as.numeric(data1[i,2])/1000000,
             #         #ybottom = j+0.1,
             ybottom = j+0.1,
             ytop = j+0.9,
             xright = as.numeric(data1[i,3])/1000000,
             #         #ytop = j+0.9,
             col = color_pad[as.numeric(as.vector(chro1[i]))],
             border = color_pad[as.numeric(as.vector(chro1[i]))] 
        )
      }
      len <- wholen[num]/1000000
      centro <- centromere[num]
      #plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
      #            if(!is.na(files[j,2])){
      #              plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
      #            }
      
      
      # text(x=c(950,950,950,935), y=c(20.5,21.5,22.5,23.5),
      #      c("lx987","nd3097","undefined","CNV DP Hete"), cex = 0.7,pos = 4)
      #            text(x=c(940,940,940,940), y=c(20.5,21.5,22.5,23.5),
      #                 c("level1","level2","level3","level4"), cex = 0.7,pos = 4)
      #text(x=c(940,940,940,940,940,940,940,940,940,940), y=c(20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5),
           #c("level10","level9","level8","level7","level6","level5","level4","level3","level2","level1"),cex = 0.7,pos = 4)
      #rect(xleft = 900, xright = 930, ybottom = 26, ytop = 28,
      #     col = color_pad[4], border = color_pad[4])
      #rect(xleft = 900, xright = 930, ybottom = 28, ytop = 30,
      #     col = color_pad[5], border = color_pad[5])
  
      # rect(xleft = 900, xright = 930, ybottom = 26, ytop = 27,
      #      col = color_pad[7], border = color_pad[7])
      # rect(xleft = 900, xright = 930, ybottom = 27, ytop = 28,
      #      col = color_pad[8], border = color_pad[8])
      # rect(xleft = 900, xright = 930, ybottom = 28, ytop = 29,
      #      col = color_pad[9], border = color_pad[9])
      # rect(xleft = 900, xright = 930, ybottom = 29, ytop = 30,
      #      col = color_pad[10], border = color_pad[10])
  
      # rect(xleft = 900, xright = 930, ybottom = 23, ytop = 24,
      #      col = color_pad[4], border = color_pad[4])
    }
  }
}
num = 22
for(j in 1:nrow(files)){
  if(is.na(files[j,1])){
    j <- j-1
    num <- num - 1
    centro <- centromere[num]
    len <- wholen[num]
    #centro = 200
    #j=3
    len <- wholen[num]/1000000
    plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
    #tx <- c(centro+5,centro,centro-5)
    #ty <- c(j-3.7,j-3,j-3.7)
    #polygon( tx, ty, border = "black", col = "black")

  }
}
# speci_pos <- function(x,y,gene){
#   tx <- c(x+2,x,x-2)
#   ty <- c(y-0.3,y,y-0.3)
#   polygon( tx, ty, border = "black", col = "white")
# }
# #
# pos_oper <- function(x){(x-1)*2+1}
# #
# ya <- pos_oper(21)#1A
# speci_pos(508,ya) 
# text(x=c(508), y=c(ya-0.7),c("Glu-A1"), cex = 0.3,adj = 0.5,col="black")
# speci_pos(4,ya)
# text(x=c(4), y=c(ya-0.7),c("Glu-A3"), cex = 0.3,adj = 0.5,col="black")
# speci_pos(381,ya)
# text(x=c(381), y=c(ya-0.7),c("SnRk2.3-1A"), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(20)#1B
# speci_pos(555,ya) 
# text(x=c(555), y=c(ya-0.7),c('Bx7OE;Bx13'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(5,ya)
# text(x=c(5), y=c(ya-0.7),c('Glu-B3'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(412,ya)
# text(x=c(412), y=c(ya-0.7),c('SnRk2.3-1B'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(19)#1D
# speci_pos(412,ya) 
# text(x=c(412), y=c(ya-0.7),c('Glu-D1'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(18)#2A
# speci_pos(712,ya) 
# text(x=c(712), y=c(ya-0.7),c('ppo-A1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(321,ya)
# text(x=c(321), y=c(ya-0.7),c('zds-A1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(121,ya)
# text(x=c(121), y=c(ya-0.7),c('SUS2-2A'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(508,ya)
# text(x=c(508), y=c(ya-0.7),c('Cwi-A1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(740,ya)
# text(x=c(740), y=c(ya-0.7),c('Flo2-A1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(740,ya)
# text(x=c(740), y=c(ya-0.7),c('Flo2-A1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(37,ya)
# text(x=c(37), y=c(ya-0.7),c('Ppd-A1'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(17)#2B
# speci_pos(171,ya) 
# text(x=c(171), y=c(ya-0.7),c('SUS2-2B'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(449,ya) 
# text(x=c(449), y=c(ya-0.7),c('bas1-B1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(55,ya) 
# text(x=c(55), y=c(ya-0.7),c('Ppd-B1'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(16)#2D
# speci_pos(119,ya) 
# text(x=c(119), y=c(ya-0.7),c('SUS2-2D'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(19,ya) 
# text(x=c(19), y=c(ya-0.7),c('Rht8'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(34,ya) 
# text(x=c(34), y=c(ya-0.7),c('Ppd-D1'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(15)#3A
# speci_pos(448,ya) 
# text(x=c(448), y=c(ya-0.7),c('lyce-3A'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(59,ya) 
# text(x=c(59), y=c(ya-0.7),c('TAR2.1-3A'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(176,ya) 
# text(x=c(176), y=c(ya-0.7),c('GS5-A1/3A'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(532,ya) 
# text(x=c(532), y=c(ya-0.7),c('TGW6-A1'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(14)#3B
# speci_pos(148,ya) 
# text(x=c(148), y=c(ya-0.7),c('Talyc-B1'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(13)#3D
# speci_pos(106,ya) 
# text(x=c(106), y=c(ya-0.7),c('CKX6-D1'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(12)#4A
# speci_pos(688,ya) 
# text(x=c(688), y=c(ya-0.7),c('Wx-B1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(330,ya) 
# text(x=c(330), y=c(ya-0.7),c('TGW6'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(610,ya) 
# text(x=c(610), y=c(ya-0.7),c('CWI-4A'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(582,ya) 
# text(x=c(582), y=c(ya-0.7),c('Rht-A1'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(11)#4B
# speci_pos(612,ya) 
# text(x=c(612), y=c(ya-0.7),c('Tapds-B1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(31,ya) 
# text(x=c(31), y=c(ya-0.7),c('Rht-B1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(657,ya) 
# text(x=c(657), y=c(ya-0.7),c('Vrn2-4B'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(10)#4D
# speci_pos(18,ya) 
# text(x=c(18), y=c(ya-0.7),c('Rht-D1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(509,ya) 
# text(x=c(509), y=c(ya-0.7),c('Vrn2-4D'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(9) #5A
# speci_pos(587,ya) 
# text(x=c(587), y=c(ya-0.7),c('Vrn1-5A'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(650,ya) 
# text(x=c(650), y=c(ya-0.7),c('Q'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(698,ya) 
# text(x=c(698), y=c(ya-0.7),c('Vrn2-5A'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(645,ya) 
# text(x=c(645), y=c(ya-0.7),c('NAC2-5A'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(8) #5B
# speci_pos(698,ya) 
# text(x=c(698), y=c(ya-0.7),c('Vrn2-5A'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(7) #5D
# speci_pos(3,ya) 
# text(x=c(3), y=c(ya-0.7),c('Pina-D1;Pinb-D1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(557,ya) 
# text(x=c(557), y=c(ya-0.7),c('CWI-5D'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(467,ya) 
# text(x=c(467), y=c(ya-0.7),c('Vrn1-5D'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(6) #6A
# speci_pos(413,ya) 
# text(x=c(413), y=c(ya-0.7),c('Rht24'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(461,ya) 
# text(x=c(461), y=c(ya-0.7),c('TPP-6AL1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(237,ya) 
# text(x=c(237), y=c(ya-0.7),c('GW2-6A'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(5) #6B
# speci_pos(291,ya) 
# text(x=c(291), y=c(ya-0.7),c('GW2-6B'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(134,ya) 
# text(x=c(134), y=c(ya-0.7),c('GPC-B1'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(4)#6D
# speci_pos(386,ya) 
# text(x=c(386), y=c(ya-0.7),c('GS1a'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(3)#7A
# speci_pos(729,ya) 
# text(x=c(729), y=c(ya-0.7),c('Psy-A1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(4,ya) 
# text(x=c(4), y=c(ya-0.7),c('6-SFT-A2'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(115,ya) 
# text(x=c(115), y=c(ya-0.7),c('SUS1-7A'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(170,ya) 
# text(x=c(170), y=c(ya-0.7),c('GASR7-A1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(205,ya) 
# text(x=c(205), y=c(ya-0.7),c('TGW-7A'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(519,ya) 
# text(x=c(519), y=c(ya-0.7),c('MOC1-7A'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(585,ya) 
# text(x=c(585), y=c(ya-0.7),c('SAP1-A1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(71,ya) 
# text(x=c(71), y=c(ya-0.7),c('Vrn3-7A'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(2)#7B
# speci_pos(739,ya) 
# text(x=c(739), y=c(ya-0.7),c('Psy-B1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(68,ya) 
# text(x=c(68), y=c(ya-0.7),c('SUS1-7B'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(9,ya) 
# text(x=c(9), y=c(ya-0.7),c('Vrn3-7B'), cex = 0.3,adj = 0.5,col="black")
# #
# ya <- pos_oper(1)#7D
# speci_pos(109,ya) 
# text(x=c(109), y=c(ya-0.7),c('SUS1-7D'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(6,ya) 
# text(x=c(6), y=c(ya-0.7),c('GS-D1'), cex = 0.3,adj = 0.5,col="black")
# speci_pos(68,ya) 
# text(x=c(68), y=c(ya-0.7),c('Vrn3-7D'), cex = 0.3,adj = 0.5,col="black")



dev.off()

