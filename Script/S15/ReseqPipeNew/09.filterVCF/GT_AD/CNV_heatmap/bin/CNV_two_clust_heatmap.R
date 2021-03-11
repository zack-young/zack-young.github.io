#!/usr/bin/env Rscript

# cgmaptools - mCBinHeatmap.R
# 
# Copyright (C) Ping Zhu
# Contact: Ping Zhu <pingzhu.work@gmail.com>
#   
#   Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.



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
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(grid))

# Arguments
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "chr1A.100k.join",
              help="input file: file name of *.join"),
  make_option(c("-M", "--metadata"), dest = "meta_in", default = "/data2/rawdata2/sample_metadata/metadata_AA.txt",
              help="metadata of samples: only samples in metadata file will appear in the plot"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file name. [default: mCBinHeatmap.SysDate.pdf]"),
  make_option(c("-c", "--cluster"), dest = "cluster", default = TRUE, action="store_true",
              help = "[opt] cluster samples by methylation in regions. [default: TRUE]"),
  make_option(c("-l","--colorLow"), dest = "low.color", default = "blue",
              help = "[opt] color used for the lowest methylation value. [default: blue]"),
  make_option(c("-m","--colorMid"), dest = "mid.color", default = "white",
              help = "[opt] color used for the middle methylation value. [default: white]"),
  make_option(c("-b","--colorHigh"), dest = "high.color", default = "red",
              help = "[opt] color used for the highest methylation value. [default: red]"),
  make_option(c("-n","--colorNumber"), dest = "num.color", default = 10,
              help = "[opt] desired number of color elements in the panel. [default: 10]"),
  make_option(c("-W","--width"), dest = "figure.width", default = 70,
              help = "[opt] width of figure (inch). [default: 70]"),
  make_option(c("-H","--height"), dest = "figure.height", default = 35,
              help = "[opt] height of figure (inch). [default: 35]"),
  make_option(c("-f","--format"), dest = "figure.format", default = "pdf",
              help = "[opt] format of output figure. Alternative: png. [default: pdf]"),
  make_option(c("-R","--resolution"), dest = "figure.resolution", default = 200,
              help = "[opt] Resolution in ppi. Only available for png format. [default: 200]")
)

parser <- OptionParser(usage = "Rscript CNV_heatmap.R [options]",
                       option_list=option_list, description = "\
                       Description: Plot DP levels of target chromo for multiple samples [heatmap]\
                       Contact:     Zhu, Ping; pingzhu.work@gmail.com\
                       Last update: 2017-09-16\
                       Example: \
                       Rscript CNV_heatmap.R -i chr1A.100k.join -m /data2/rawdata2/sample_metadata/metadata_AABBDD.txt \
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
infile <- arguments$infile
if(infile == ""){ # default, STDIN
  infile <- file("stdin")
} else { # user specified
  if( file.access(infile) == -1){ # file not exists
    print_help(parser)
  }
}
outfile <- arguments$outfile
if(outfile == ""){ # default, "FragRegView.Date"
  outfile <- infile
} else { # user specified
  outfile <- gsub(".pdf$|.png$", "", outfile, perl=T)
}
figure.width <- arguments$figure.width
figure.height <- arguments$figure.height
figure.format <- arguments$figure.format
if(! figure.format %in% c("pdf", "png")){ # format not support
  print_help(parser)
} else {
  outfile <- paste(outfile, figure.format, sep = ".")
}
figure.resolution <- arguments$figure.resolution
low.color <- arguments$low.color
mid.color <- arguments$mid.color
if(mid.color == ""){
  mid.color <- NA
}
high.color <- arguments$high.color
num.color <- arguments$num.color

###
meta_in <- arguments$meta_in
#
low.color <- 'blue'
mid.color <- 'white'
high.color <- 'red'
num.color <- 10


meta_in <- "/data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_final_199_headed.txt"
infile <- "/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/chr1B.1M.join_norm"
sample_lis = read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/1B1R_lovrin/1B_1R_time_list", as.is = T, header = F, comment.char = "")

#infile <- "/data/user/yangzz/mapping/fieldergenomecompare/20200416/C12_2020-04-17162314/combine_10M_DP"
#meta_in <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/metadata_cultivar_nogroup.txt"
#infile <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/cultivar_combine_10M_norm"
#
meta_data <- read_delim(meta_in, "\t", escape_double = FALSE, trim_ws = TRUE, comment = "#")
data <- read_delim(infile, "\t", escape_double = FALSE, trim_ws = TRUE)

## filter sample
tmp_index <- colnames(data) %in% meta_data$Sample
data <- data[,tmp_index]

for (num in seq(1:ncol(data))) {
  data[,num][which((data[,num] >= 0.5)&(data[,num] < 1.5)), ] = 1
  data[,num][which(data[,num] < 0.5), ] = 0
  #data[,num][which(data[,num] >= 0.5), ] = 1
  data[,num][which(data[,num] >= 1.5), ] = 2
}

sample_present <- meta_data[meta_data$Sample %in% colnames(data),]
LIST <- sample_present$label
names(LIST) <- sample_present$Sample

colnames(data) = LIST[colnames(data)]

## label color
colorF <- '/data/user/yangzz/mapping/fieldergenomecompare/coloryzz_6cols.txt'
colorFile <- read.table(colorF, header = T, stringsAsFactors = F, sep = "\t", comment.char = "")
colors <- colorFile[,2]
names(colors) <- colorFile[,1]
#colors = c("#5fb951", "#E74D4A", "#D9B460", "#AB6193","#71ACDF", "#F194BE")#, "#424E9F")
#names(colors) = unique(sample_present$Region)

## figure device


## cluster specified
#cluster <- arguments$cluster



## figure margin 

#correlation <- cor(data, method="spearman", use="pairwise.complete.obs")

hc <- hclust(d = dist(t(data)))
dendrogram <- as.dendrogram(hc,hang = -1)

#pdf("/data/user/yangzz/mapping/fieldergenomecompare/pdf/whole_genome_CNV_clust.pdf", height = 10, width = 13)
par(mar=c(3, 0, 0, 0), oma=c(1,1,1,1))
#layout(matrix(1:2, nrow=1), widths=c(1,2))
layout(matrix(c(1,0,2,3), nrow=2), widths=c(1,2,2,2),heights = c(3,1,3,1))
plot(dendrogram, horiz=T,leaflab="none",  dLeaf = 10 ,edgePar=list(lwd=1), yaxt = "none", yaxs="i")#,frame.plot = T)#,xlim = c(14,0))#


segments(x0=seq(0, 18, 2),
         x1=seq(0, 18, 2),lwd=2,
         y0 = -0.5,
         y1 = -1.5,
         xpd = T
)
segments(x0=0,
         x1 = 18,
         y0 = -0.5,
         y1 = -0.5,lwd=2,
         xpd = T
)
text(x=seq(0, 18, 2),
     y=rep(-2, 17),
     seq(0,18, 2),
     adj = c(0.5, 1),
     cex= 1,
     xpd = T
)

#
row_index = c()
for (num in seq(1:nrow(data))) {
  row_index = c(row_index,length(data[num,][, which((data[num,] == 0)|(data[num,] == 2))]))
}

b = which(row_index<=1)
#row_order = seq(1:length(row_index))
row_order <- order(row_index,decreasing = T)
if (length(b)==0) {
  row_order_n <- row_order
} else{
  row_order_n <-row_order[-match(b, row_order)]
}

#data <- data[rev(row_order_n), order.dendrogram(dendrogram)]
#data <- data[row_order_n, order.dendrogram(dendrogram)]
row_order_n <- data[,1]
#data <- data[, order.dendrogram(dendrogram)]
tree_lis <- order.dendrogram(dendrogram)
B_lis <- match(LIST[sample_lis[[1]]],colnames(data))
data <- data[,c(B_lis,tree_lis[-match(B_lis ,order.dendrogram(dendrogram))])]
#

# hc_n <- hclust(d = dist(data),method="average")
# dendrogram_n <- as.dendrogram(hc_n)
# data <- data[order.dendrogram(dendrogram_n), ]

## empty figure panel
#par(mar=c(3, 2, 0, 0), oma=c(1,1,1,1))
#par(mar=c(5, 0, 4, 0),xpd=T)


plot(x=0,  xaxt="n", yaxt="n", type="n", #bty="n",
     xlab="", ylab="", 
     xlim=c(-40, nrow(data)+10), ylim=c(0.5, ncol(data)+1), xaxs="i", yaxs="i",frame.plot = F)
#xlim=c(-10, length(row_order_n)+30)
## methylation to colors

#color.rects <- colorpanel(low=low.color, mid=mid.color, high=high.color, n=num.color)
color.rects <- c('#08306B','#E3E4E6','red')
#methylation.interval <- seq(0, 2, length.out=num.color+1)
row.height <- 0.8
rightPanel.extension <- 1.15
data <- as.data.frame(data)
for(i in 1:ncol(data)){
  #y.i <- (ncol(data):1)[i]
  y.i <- i-0.5
  #x.color <- color.rects[ findInterval(data[,i], methylation.interval, all.inside=T)]
  x.color <- color.rects[data[,i]+1]
  # heatmap body
  rect(xleft = 1:nrow(data) + 0.0 * nrow(data), 
       ybottom = rep(y.i, nrow(data)) , 
       xright = 1:nrow(data) + 0.0 * nrow(data) + 1, 
       ytop = rep(y.i, nrow(data)) + 1,
       col = x.color,
       border = x.color
  )
  
  # methylation level on the right panel, barplot
  # met.mean <- mean(as.numeric(data[,i]), na.rm = T)
  # rect(xleft = nrow(data) + 0.08 * nrow(data) + 1,
  #      ybottom = y.i + 0.1,
  #      xright = met.mean * nrow(data)*(rightPanel.extension - 1.12) + nrow(data) + 0.08 * nrow(data) + 1,
  #      ytop = y.i + 0.1 + row.height,
  #      col = "#66CD00",
  #      border = "#66CD00"
  # )
  # text(x= met.mean * nrow(data)*(rightPanel.extension - 1.12) + nrow(data) + 0.08 * nrow(data) + 2,
  #      y=y.i+0.1+row.height/2,
  #      paste(format(met.mean*100,digits = 3), "%", sep = ""),
  #      cex = 0.7, adj = c(0, 0.5)
  #)
}
segments(x0=seq(1, (length(row_order_n)/10+1)*10, 100),
         x1=seq(1, (length(row_order_n)/10+1)*10, 100),lwd=2,
         y0 = -0.5,
         y1 = -1.5,
         xpd = T
)
segments(x0=1,
         x1 = (length(row_order_n)/10+1)*10,
         y0 = -0.5,
         y1 = -0.5,lwd=2,
         xpd = T
)
text(x=seq(1, (length(row_order_n)/10+1)*10+1, 100),
     y=rep(-1.5, length(seq(1, (length(row_order_n)/10+1)*10, 100))),
     seq(0,(length(row_order_n)/10+1)*10, 100),
     adj = c(0.5, 1),
     cex= 1,
     xpd = T
)
## sample label
sampleLabel.col <- "black"
rect(xleft = rep(-20,ncol(data)), 
     ybottom =0.5:ncol(data), 
     xright = rep(0,ncol(data)), 
     ytop =0.5:ncol(data)+1,
     col = colors[sample_present$Variety[match(colnames(data), sample_present$label)]],
     border =colors[sample_present$Variety[match(colnames(data), sample_present$label)]])

#points(x=rep(-0.5,ncol(data)), y=1:ncol(data)+row.height/2 + 0.1,  pch=16,cex=2.5, xpd=T, col=colors[sample_present$Region[match(colnames(data), sample_present$label)]])

legend(x=nrow(data)+2,y=ncol(data)*3/4, ncol=2,fill = colors, legend = names(colors), border = F, box.col = "white", horiz=F, cex = 1, x.intersp = 0,xpd=T)
plot(x=0,  xaxt="n", yaxt="n", type="n", #bty="n",
     xlab="", ylab="", 
     xlim=c(-40, nrow(data)+10), ylim=c(1, 10), xaxs="i", yaxs="i",frame.plot=F)

plotChrom(1, 9, 0.5,  wholen[2]/1000000 , centromere[2] , 1000000)

dev.off()

  
  #text(x=0, y=1:ncol(data)+row.height/2 + 0.1, names(data), cex=0.8, xpd=T, adj=c(0,0.5), col=colors[sample_present$Region[match(colnames(data), sample_present$label)]])
  
  #text(x=1, y=1:ncol(data)+row.height/2 + 0.1, names(data), cex=0.8, xpd=T, adj=c(1,0.5))
layout(matrix(1))
del_index=c()
dup_index=c()
for (num in seq(1:nrow(data))) {
  del_index = c(del_index,length(data[num,][, which(data[num,] == 0)]))
  dup_index = c(dup_index,length(data[num,][, which(data[num,] == 2)]))
}
del_index <- data.frame(V1=c(1:nrow(data)),V2=del_index,V3="del")
dup_index <- data.frame(V1=c(1:nrow(data)),V2=dup_index,V3="dup")
dt_total <- rbind(del_index,dup_index)
dt_total[which(dt_total[,2] == 0), ][,2] = 1
pdf("/data/user/yangzz/pdf/CNV_heatmap_dist.pdf", height = 8, width = 13)
ggplot(dt_total,aes(V1,log(V2,2),fill=V3))+
  geom_area(position="identity",alpha=1)+
  scale_y_continuous(limits=c(0,8), breaks=seq(0,8,2))+
  scale_fill_manual(values = c("#4889E6", "#FF7F7F"))+theme_classic()
#plot(x=0,xlim=c(0,200),ylim=c(0,200),type='l',col='blue',xaxt='n',xlab = '',yaxt='n', lwd=2)
#polygon(x=c(0,1:nrow(data),nrow(data)),y=c(0,del_index,0))
#spline(x=1:nrow(data),y=del_index,n=100)
#axis(2,cex.axis=2)
#axis(1,cex.axis=2)
#lines(x=1:nrow(data),y=dup_index,type='l',col='red',lwd=2)
#spline(x=1:nrow(data),y=dup_index,n=100)
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
  polygon(lx, ly, border = "white", col = "white",xpd=T)
  
  # right semi-circle
  rx <- c(xleft + len,xleft + len, xleft + len - rh + rh*sin(rs), xleft + len)
  ry <- c(ybottom, ybottom + height, ybottom + rv +rv*cos(rs), ybottom)
  polygon(rx, ry, border = "white", col = "white",xpd=T)
  
  # top tri
  tx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  ty <- c(ybottom + height, ybottom + rv, ybottom + height)
  polygon(tx, ty, border = "white", col = "white",xpd=T)
  
  # bottom tri
  bx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  by <- c(ybottom, ybottom + rv, ybottom)
  polygon(bx, by, border = "white", col = "white",xpd=T)
  
  cbx <- c(xleft + rh - rh*sin(rs), centro - 0.5*rh, centro, centro + 0.5*rh, xleft + len - rh + rh*sin(rs), centro + 0.5*rh, centro, centro - 0.5*rh)
  cby <- c(ybottom + rv - rv*cos(rs), ybottom + height, ybottom + rv, ybottom + height, ybottom + rv + rv*cos(rs), ybottom, ybottom + rv, ybottom)
  polygon(cbx, cby, border = "grey", lwd = 5,xpd=T)
}
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
num=2
centro <- centromere[num]/10
len <- wholen[num]
#centro = 200
#j=3
len <- (wholen[num]/10000000)+1
plotChrom(1+0.0 * nrow(data), -2, 2, len, centro, 10000000)




## methylation value legend
colorLegend.y.min <- nrow(meta_data)-10
colorLegend.y.max <- nrow(meta_data)-30
colorLegend.y <- seq(from=colorLegend.y.max, to=colorLegend.y.min, length.out = num.color+1)
colorLegend.x <- -4
rect(xleft = 88,
     xright = 90,
     ybottom = colorLegend.y[-1], 
     ytop = colorLegend.y[-length(colorLegend.y)],
     col = color.rects,
     border = color.rects,
     xpd = T
)

colorLegend <-seq(from=nrow(meta_data)-50, to=nrow(meta_data)-70, length.out = length(unique(sample_present$Region)))
points(x = rep(89,length(unique(sample_present$Region))),
       y = colorLegend,
       col=colors[1:length(unique(sample_present$Region))],
       cex=8,
       pch=16,
       xpd = T
)
# value label
sta <- (colorLegend.y[2]+colorLegend.y[1])/2
en <- (colorLegend.y[length(colorLegend.y)]+colorLegend.y[length(colorLegend.y)-1])/2
text(x=rep(90.5,3),
     y=seq(from=sta, to=en, length.out = 3),
     seq(0,2,length.out = 3),
     adj = c(0.5, 1),
     cex= 4,
     pos = 4
)

# text label
text(x = rep(90,length(unique(sample_present$Region))),
     y = colorLegend,
     col=colors[1:length(unique(sample_present$Region))],
     names(colors)[1:length(unique(sample_present$Region))],
     cex=3,
     pos = 4,
     xpd = T
)

invisible(dev.off())
infile <- "/data/user/yangzz/mapping/fieldergenomecompare/statistic/1B1R_lovrin/chr1B.lovrin_simmilar.dist"
#meta_in <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/metadata_cultivar_nogroup.txt"

meta_data <- read_delim(meta_in, "\t", escape_double = FALSE, trim_ws = TRUE, comment = "#")
data <- read_delim(infile, "\t", escape_double = FALSE, trim_ws = TRUE)
tmp <- data[ ,-c(1,2,3)]
tmp <- subset(tmp, select = -S212)

tmp[tmp=='low']='2'
tmp[tmp=='deletion_both_CNV']='1'
tmp[tmp=='duplication_both_CNV']='1'
tmp[tmp!='2'&tmp!='1']='0'
par(mar=c(3, 3, 3, 3), oma=c(1,1,1,1))
layout(matrix(1))
plot(x=0,  xaxt="n", yaxt="n", type="n", #bty="n",
     xlab="", ylab="", 
     xlim=c(-10, nrow(tmp)+10), ylim=c(1, ncol(tmp)+1), xaxs="i", yaxs="i",frame.plot = F)
axis(1)
tmp <- as.data.frame(tmp)
color.rects <- c('#E3E4E6','#5C957F','#22BCF0')

hc <- hclust(d = dist(t(tmp)))
dendrogram <- as.dendrogram(hc,hang = -1)
tmp <- tmp[, order.dendrogram(dendrogram)]


for(i in 1:ncol(tmp)){
  #y.i <- (ncol(data):1)[i]
  y.i <- i
  #x.color <- color.rects[ findInterval(data[,i], methylation.interval, all.inside=T)]
  x.color <- color.rects[as.numeric(tmp[,i])+1]
  # heatmap body
  rect(xleft = 1:nrow(tmp) + 0.0 * nrow(tmp), 
       ybottom = rep(y.i, nrow(tmp)) , 
       xright = 1:nrow(tmp) + 0.0 * nrow(tmp) + 1, 
       ytop = rep(y.i, nrow(tmp)) + 1,
       col = x.color,
       border = x.color
  )

}
text(x=rep(-4,ncol(tmp)),y=seq(1,ncol(tmp)),labels = colnames(tmp),cex = 0.5)
