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
  make_option("--sample_ma", dest = "sample_ma", default = "",
              help = "sample1 serial number"),
  make_option("--sample_ma_name", dest = "sample_ma_name", default = "",
              help = "sample1 name"),
  make_option("--sample_fa", dest = "sample_fa", default = "",
              help = "sample2 serial number"),
  make_option("--sample_fa_name", dest = "sample_fa_name", default = "",
              help = "sample2 name"),
  make_option("--sample_son", dest = "sample_son", default = "",
              help = "sample2 serial number"),
  make_option("--sample_son_name", dest = "sample_son_name", default = "",
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
SAMPLE_fa=arguments$sample_fa
SAMPLE_fa_name=arguments$sample_fa_name
SAMPLE_ma=arguments$sample_ma
SAMPLE_ma_name=arguments$sample_ma_name
SAMPLE_son=arguments$sample_son
SAMPLE_son_name=arguments$sample_son_name
options(scipen=200)

files = read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/lx987_5181_3097/plotfile_",SAMPLE_fa,"_",SAMPLE_ma,"_",SAMPLE_son,sep=""), as.is = T, header = F, comment.char = "")

pdf(paste("/data/user/yangzz/mapping/fieldergenomecompare/lx987_5181_3097/pdf/",SAMPLE_fa_name,"_",SAMPLE_ma_name,"_",SAMPLE_son_name,arguments$suffix,".pdf",sep=""), height = figure.height, width = figure.width)
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
     xlim=c(-200, 1100), ylim=c(-2, nrow(files)*2),
     xaxs="i", yaxs="i", 
     #main=paste("S",j,sep=""),
     main=paste(SAMPLE_fa_name,"_",SAMPLE_ma_name,"_",SAMPLE_son_name,sep=" "),
     cex.main = 2.5)
axis(1,at=c(0,100,200,300,400,500,600,700,800),labels = c("0Mbp","","200Mbp","","400Mbp","","600Mbp","","800Mbp"))
#text(x=c(0,200,400,600,800), y=rep(0,5),
#     c("0Mbp","200Mbp","400Mbp","600Mbp","800Mbp"), cex = 1,pos = 1)
#color_pad <- c(colorpanel(low="yellow", high="red", n=3),"green","black","white")

del_fa <- "#31C331"
del_ma <- "#02F791"
del_son <- "#88CD33"
delboth <-  "#005A00"
dup_fa <- "#C1431B"
dup_ma <- "#F3AE83"
dup_son <- '#C8326F'
dupboth <- "#710000"
high <- "#111148"
mid <-"#334899"
low <- "#22BBEE"
unknow <- "#B3B1B4"
fa <- '#B40679'
ma <- "#006363"
both <- "#450063"
fa_cnv <- "#AD66D5"
ma_cnv <- "#5CCCCC"
both_cnv <- "#FFB273"
unknow_cnv <- '#3F3A42'

#color_pad <- c(deepblue,midblue,lightblue,,red,yellow,purple)
#color_pad <- c('white','white','white',red,yellow,purple)
#color_pad <- c(deepblue,midblue,lightblue,'white','white','white')
color_pad_fa <- c(low,mid,high,del_fa,del_son,delboth,dup_fa,dup_son,dupboth)
color_pad_ma <- c(low,mid,high,del_ma,del_son,delboth,dup_ma,dup_son,dupboth)
color_pad_son <- c(fa,ma,both,unknow,fa_cnv,ma_cnv,both_cnv,unknow_cnv)

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
# rect(xleft = 900, xright = 930, ybottom = 18, ytop = 20,
#    col = color_pad[1], border = color_pad[1]);text(x=c(940), y=c(19),c("Similar"),cex = 1,pos = 4)
# rect(xleft = 900, xright = 930, ybottom = 20, ytop = 22,
#    col = color_pad[2], border = color_pad[2]);text(x=c(940), y=c(21),c("Different"),cex = 1,pos = 4)

rect(xleft = 900, xright = 930, ybottom = 8, ytop = 10,
     col = unknow, border = unknow);text(x=c(940), y=c(9),c("Unknown"),cex = 1,pos = 4)

rect(xleft = 900, xright = 930, ybottom = 10, ytop = 12,
     col = both, border = both);text(x=c(940), y=c(11),c("Parent_share"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 12, ytop = 14,
     col = fa, border = fa);text(x=c(940), y=c(13),c(SAMPLE_fa_name),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 14, ytop = 16,
     col = ma, border = ma);text(x=c(940), y=c(15),c(SAMPLE_ma_name),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 16, ytop = 18,
     col = low, border = low);text(x=c(940), y=c(17),c("Low_Diff"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 18, ytop = 20,
   col = mid, border = mid);text(x=c(940), y=c(19),c("Mid_Diff"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 20, ytop = 22,
   col = high, border = high);text(x=c(940), y=c(21),c("High_Diff"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 22, ytop = 24,
     col = del_fa, border = del_fa);text(x=c(940), y=c(23),c(paste(SAMPLE_fa_name,"deletion",sep="")),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 24, ytop = 26,
     col = del_ma, border = del_ma);text(x=c(940), y=c(25),c(paste(SAMPLE_ma_name,"deletion",sep="")),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 26, ytop = 28,
     col = del_son, border = del_son);text(x=c(940), y=c(27),c(paste(SAMPLE_son_name,"deletion",sep="")),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 28, ytop = 30,
     col = delboth, border =  delboth);text(x=c(940), y=c(29),c("Both_deletion"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 30, ytop = 32,
     col = dup_fa, border =  dup_fa);text(x=c(940), y=c(31),c(paste(SAMPLE_fa_name,"duplication",sep="")),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 32, ytop = 34,
     col =  dup_ma, border =  dup_ma);text(x=c(940), y=c(33),c(paste(SAMPLE_ma_name,"duplication",sep="")),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 34, ytop = 36,
     col =  dup_son, border =  dup_son);text(x=c(940), y=c(35),c(paste(SAMPLE_son_name,"duplication",sep="")),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 36, ytop = 38,
     col = dupboth, border = dupboth);text(x=c(940), y=c(37),c("Both_duplication"),cex = 1,pos = 4)
num_line=0
for(j in 1:nrow(files)){
  if(!is.na(files[j,1])){
    data1 <- read.table(files[j,1], as.is = T, header = T, comment.char = "")
    lin = num_line*2 + j
    num_line=num_line+1
    chro <- factor(data1[,4])
    chro_fa <- factor(chro, 
                    levels = c("low","mid","high",
                               paste(SAMPLE_fa,"deletion_CNV",sep=""),paste(SAMPLE_son,"deletion_CNV",sep=""),"deletion_both_CNV",
                               paste(SAMPLE_fa,"duplication_CNV",sep=""),paste(SAMPLE_son,"duplication_CNV",sep=""),"duplication_both_CNV"),
                    labels = c('1','2','3',"4","5",'6',"7","8","9")
    )
    chro <- factor(data1[,5])
    chro_ma <- factor(chro, 
                      levels = c("low","mid","high",
                                 paste(SAMPLE_ma,"deletion_CNV",sep=""),paste(SAMPLE_son,"deletion_CNV",sep=""),"deletion_both_CNV",
                                 paste(SAMPLE_ma,"duplication_CNV",sep=""),paste(SAMPLE_son,"duplication_CNV",sep=""),"duplication_both_CNV"),
                      labels = c('1','2','3',"4","5",'6',"7","8","9")
    )
    chro <- factor(data1[,6])
    chro_son <- factor(chro, 
                    levels = c(SAMPLE_fa,SAMPLE_ma,"both","unknown",
                               paste(SAMPLE_fa,"_CNV",sep=""),paste(SAMPLE_ma,"_CNV",sep=""),"both_CNV",
                               "unknow_CNV"),
                    labels = c('1','2','3',"4","1",'2',"3","4")
                    )
    num <- num - 1 
    for(i in 1:nrow(data1)){
      rect(xleft = as.numeric(data1[i,2])/1000000,
           #         #ybottom = j+0.1,
           ybottom = lin+0.1,ytop = lin+0.9,
           xright = as.numeric(data1[i,3])/1000000,
           #         #ytop = j+0.9,
           col = color_pad_son[as.numeric(as.vector(chro_son[i]))],
           border = color_pad_son[as.numeric(as.vector(chro_son[i]))])
      rect(xleft = as.numeric(data1[i,2])/1000000,
           #         #ybottom = j+0.1,
           ybottom = lin+1.1,ytop = lin+1.9,
           xright = as.numeric(data1[i,3])/1000000,
           #         #ytop = j+0.9,
           col = color_pad_ma[as.numeric(as.vector(chro_ma[i]))],
           border = color_pad_ma[as.numeric(as.vector(chro_ma[i]))])
      rect(xleft = as.numeric(data1[i,2])/1000000,
           #         #ybottom = j+0.1,
           ybottom = lin+2.1,ytop = lin+2.9,
           xright = as.numeric(data1[i,3])/1000000,
           #         #ytop = j+0.9,
           col = color_pad_fa[as.numeric(as.vector(chro_fa[i]))],
           border = color_pad_fa[as.numeric(as.vector(chro_fa[i]))])

    }
    len <- wholen[num]/1000000
    centro <- centromere[num]
    #plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
    #            if(!is.na(files[j,2])){
    #              plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
    #            }
    
    text(x=-100, y=lin+0.5, files[j,2], cex = 1)
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
num = 22
num_line=0
for(j in 1:nrow(files)){
  if(is.na(files[j,1])){
    lin = num_line*2 + j
    num_line=num_line+1
    j <- j-1
    num <- num - 1
    centro <- centromere[num]
    len <- wholen[num]
    #centro = 200
    #j=3
    # len <- wholen[num]/1000000
    # plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
    tx <- c(centro+5,centro,centro-5)
    ty <- c(lin-1.7,lin-1,lin-1.7)
    polygon( tx, ty, border = "black", col = "black")

  }
}
dic = c("A"=0, "B"=1, "D"=2)
pos_file = read.table("/data/user/yangzz/mapping/fieldergenomecompare/gene.txt", as.is = T, header =T, comment.char = "")

speci_pos <- function(x,y,gene){
  tx <- c(x+2,x,x-2)
  ty <- c(y-0.3,y,y-0.3)
  polygon( tx, ty, border = "black", col = "white")
  text(x=c(x), y=c(y-0.7),c(gene), cex = 0.3,adj = 0.5,col="black")
}
#
pos_oper <- function(x){4*x-3}
#
for(j in 1:nrow(pos_file)){

  sub_genome <-  substr(pos_file[j,1], 2,2)
  chr_num <- as.numeric(substr(pos_file[j,1], 1,1))
  pos_num <- 3*(chr_num -1)+  unname(dic[sub_genome])
  ya <- pos_oper(21-pos_num)
  speci_pos(pos_file[j,2],ya,pos_file[j,3])
}
# ya <- pos_oper(21)#1A
# speci_pos(508,ya,"Glu-A1")
# speci_pos(4,ya,"Glu-A3")
# speci_pos(381,ya,"SnRk2.3-1A")
# #
# ya <- pos_oper(20)#1B
# speci_pos(555,ya,'Bx7OE;Bx13')
# speci_pos(5,ya,'Glu-B3')
# speci_pos(412,ya,'SnRk2.3-1B')
# #
# ya <- pos_oper(19)#1D
# speci_pos(412,ya,'Glu-D1')
# #
# ya <- pos_oper(18)#2A
# speci_pos(712,ya,'ppo-A1')
# speci_pos(321,ya,'zds-A1')
# speci_pos(121,ya,'SUS2-2A')
# speci_pos(508,ya,'Cwi-A1')
# speci_pos(740,ya,'Flo2-A1')
# speci_pos(37,ya,'Ppd-A1')
# #
# ya <- pos_oper(17)#2B
# speci_pos(171,ya,'SUS2-2B')
# speci_pos(449,ya,'bas1-B1')
# speci_pos(55,ya,'Ppd-B1')
# #
# ya <- pos_oper(16)#2D
# speci_pos(119,ya,'SUS2-2D')
# speci_pos(19,ya,'Rht8')
# speci_pos(34,ya,'Ppd-D1')
# #
# ya <- pos_oper(15)#3A
# speci_pos(448,ya,'lyce-3A')
# speci_pos(59,ya,'TAR2.1-3A')
# speci_pos(176,ya,'GS5-A1/3A')
# speci_pos(532,ya,'TGW6-A1')
# #
# ya <- pos_oper(14)#3B
# speci_pos(148,ya,'Talyc-B1')
# #
# ya <- pos_oper(13)#3D
# speci_pos(106,ya,'CKX6-D1')
# #
# ya <- pos_oper(12)#4A
# speci_pos(688,ya,'Wx-B1')
# speci_pos(330,ya,'TGW6')
# speci_pos(610,ya,'CWI-4A')
# speci_pos(582,ya,'Rht-A1')
# #
# ya <- pos_oper(11)#4B
# speci_pos(612,ya,'Tapds-B1')
# speci_pos(31,ya,'Rht-B1')
# speci_pos(657,ya,'Vrn2-4B')
# #
# ya <- pos_oper(10)#4D
# speci_pos(18,ya,'Rht-D1')
# speci_pos(509,ya,'Vrn2-4D')
# #
# ya <- pos_oper(9) #5A
# speci_pos(587,ya,'Vrn1-5A')
# speci_pos(650,ya,'Q')
# speci_pos(698,ya,'Vrn2-5A')
# speci_pos(645,ya,'NAC2-5A')
# #
# ya <- pos_oper(8) #5B
# speci_pos(698,ya,'Vrn2-5A')
# #
# ya <- pos_oper(7) #5D
# speci_pos(3,ya,'Pina-D1;Pinb-D1')
# speci_pos(557,ya,'CWI-5D')
# speci_pos(467,ya,'Vrn1-5D')
# #
# ya <- pos_oper(6) #6A
# speci_pos(413,ya,'Rht24')
# speci_pos(461,ya,'TPP-6AL1')
# speci_pos(237,ya,'GW2-6A')
# #
# ya <- pos_oper(5) #6B
# speci_pos(291,ya,'GW2-6B')
# speci_pos(134,ya,'GPC-B1')
# #
# ya <- pos_oper(4)#6D
# speci_pos(386,ya,'GS1a')
# #
# ya <- pos_oper(3)#7A
# speci_pos(729,ya,'Psy-A1')
# speci_pos(4,ya,'6-SFT-A2')
# speci_pos(115,ya,'SUS1-7A')
# speci_pos(170,ya,'GASR7-A1')
# speci_pos(205,ya,'TGW-7A')
# speci_pos(519,ya,'MOC1-7A')
# speci_pos(585,ya,'SAP1-A1')
# speci_pos(71,ya,'Vrn3-7A')
# #
# ya <- pos_oper(2)#7B
# speci_pos(739,ya,'Psy-B1')
# speci_pos(68,ya,'SUS1-7B')
# speci_pos(9,ya,'Vrn3-7B')
# #
# ya <- pos_oper(1)#7D
# speci_pos(109,ya,'SUS1-7D')
# speci_pos(6,ya,'GS-D1')
# speci_pos(68,ya,'Vrn3-7D')




dev.off()

