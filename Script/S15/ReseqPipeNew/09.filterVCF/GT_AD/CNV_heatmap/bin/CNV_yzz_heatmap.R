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
meta_in <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/metadata_cultivar_final.txt"
infile <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/chr1B.10M.join_norm"
#
meta_data <- read_delim(meta_in, "\t", escape_double = FALSE, trim_ws = TRUE, comment = "#")
data <- read_delim(infile, "\t", escape_double = FALSE, trim_ws = TRUE)

## filter sample
tmp_index <- colnames(data) %in% meta_data$Sample
data <- data[,tmp_index]

sample_present <- meta_data[meta_data$Sample %in% colnames(data),]
LIST <- sample_present$label
names(LIST) <- sample_present$Sample

colnames(data) = LIST[colnames(data)]

## label color
colors = c("#5fb951", "#E74D4A", "#D9B460", "#AB6193","#71ACDF", "#F194BE", "#424E9F")
names(colors) = unique(sample_present$Region)

## figure device
if (figure.format == "png"){
  png(paste0(outfile, ".png"), type="cairo", height = figure.height, width = figure.width, res = figure.resolution, units = "in")
} else if (figure.format == "pdf"){
  pdf(outfile, height = figure.height, width = figure.width)
}

## cluster specified
cluster <- arguments$cluster
pdf("/data/user/yangzz/pdf/Rplot01.pdf", height = 30, width = 35)

if (cluster){
  ## figure margin 
  par(mar=c(15, 1, 4, 0), oma=c(0,4,2,4))
  layout(matrix(1:2, nrow=1), widths=c(3,5))
  #correlation <- cor(data, method="spearman", use="pairwise.complete.obs")
  hc <- hclust(d = dist(t(data)))
  dendrogram <- as.dendrogram(hc)
  plot(dendrogram, horiz=T, leaflab="none", yaxt = "none", dLeaf = 0, yaxs="i")
  
  # reorder sample index
  data <- data[, order.dendrogram(dendrogram)]
} else {
  ## figure margin 
  par(mar=c(5, 6, 4.1, 2.1))
}

## empty figure panel
par(mar=c(15, 1, 4, 11))
plot(x=0, type="n", bty="n", xaxt="n", yaxt="n", 
     xlab="", ylab="", 
     xlim=c(0, 98), ylim=c(1, ncol(data)+1), xaxs="i", yaxs="i")

## methylation to colors
if ( is.na(mid.color) ){
  color.rects <- colorpanel(low=low.color, high=high.color, n=num.color)
} else {
  color.rects <- colorpanel(low=low.color, mid=mid.color, high=high.color, n=num.color)
}
methylation.interval <- seq(0, 2, length.out=num.color+1)
row.height <- 0.8
rightPanel.extension <- 1.15
data <- as.data.frame(data)
for(i in 1:ncol(data)){
  #y.i <- (ncol(data):1)[i]
  y.i <- i
  x.color <- color.rects[ findInterval(data[,i], methylation.interval, all.inside=T)]
  
  # heatmap body
  rect(xleft = 1:nrow(data) + 0.0 * nrow(data), 
       ybottom = rep(y.i, nrow(data)) + 0.1, 
       xright = 1:nrow(data) + 0.0 * nrow(data) + 1, 
       ytop = rep(y.i, nrow(data)) + 0.1 + row.height,
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

 ## sample label
sampleLabel.col <- "black"
if (cluster){
  points(x=rep(0.5,ncol(data)), y=1:ncol(data)+row.height/2 + 0.1,  pch=16,cex=2.5, xpd=T, col=colors[sample_present$Region[match(colnames(data), sample_present$label)]])
  
  #text(x=0, y=1:ncol(data)+row.height/2 + 0.1, names(data), cex=0.8, xpd=T, adj=c(0,0.5), col=colors[sample_present$Region[match(colnames(data), sample_present$label)]])
} else {
  points(x=rep(0.5,ncol(data)), y=1:ncol(data)+row.height/2 + 0.1, cex=2.5, xpd=T, adj=c(1,0.5))
  
  #text(x=1, y=1:ncol(data)+row.height/2 + 0.1, names(data), cex=0.8, xpd=T, adj=c(1,0.5))
}


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
  polygon(cbx, cby, border = "grey", lwd = 10,xpd=T)
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
segments(x0=seq(1+0.0 * nrow(data),1+0.0 * nrow(data)+85,length.out = 18),
         x1=seq(1+0.0 * nrow(data),1+0.0 * nrow(data)+85,length.out = 18),
         y0 = -3,
         y1 = -4,
         xpd = T
)
segments(x0=1+0.0 * nrow(data),
         x1 = 1+0.0 * nrow(data)+85,
         y0 = -3,
         y1 = -3,
         xpd = T
)
text(x=seq(1+0.0 * nrow(data),1+0.0 * nrow(data)+85,length.out = 18),
     y=rep(-4.5, 17),
     seq(0,85,length.out = 18),
     adj = c(0.5, 1),
     cex= 4,
     xpd = T
)



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

