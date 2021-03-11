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

# Arguments
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "chr1A",
              help="input file"),
  make_option(c("-M", "--metadata"), dest = "meta_in", default = "../metadata_DD.txt",
              help="metadata of samples"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file name. [default: mCBinHeatmap.SysDate.pdf]"),
  make_option(c("-c", "--cluster"), dest = "cluster", default = TRUE, action="store_true",
              help = "[opt] cluster samples by methylation in regions. [default: TRUE]"),
  make_option(c("-l","--colorLow"), dest = "low.color", default = "blue",
              help = "[opt] color used for the lowest methylation value. [default: cyan3]"),
  make_option(c("-m","--colorMid"), dest = "mid.color", default = "#E5E5E5",
              help = "[opt] color used for the middle methylation value. [default: null]"),
  make_option(c("-b","--colorHigh"), dest = "high.color", default = "red",
              help = "[opt] color used for the highest methylation value. [default: coral2]"),
  make_option(c("-n","--colorNumber"), dest = "num.color", default = 10,
              help = "[opt] desired number of color elements in the panel. [default: 10]"),
  make_option(c("-W","--width"), dest = "figure.width", default = 40,
              help = "[opt] width of figure (inch). [default: 7]"),
  make_option(c("-H","--height"), dest = "figure.height", default = 70,
              help = "[opt] height of figure (inch). [default: 7]"),
  make_option(c("-f","--format"), dest = "figure.format", default = "pdf",
              help = "[opt] format of output figure. Alternative: pdf. [default: pdf]"),
  make_option(c("-R","--resolution"), dest = "figure.resolution", default = 200,
              help = "[opt] Resolution in ppi. Only available for png format. [default: 200]"),
  make_option(c("-U","--upper"), dest = "upper", default = 0,
              help = "[opt] upper border on chr to view, in MB. [default: 0]"),
  make_option(c("-L","--lower"), dest = "lower", default = 100,
              help = "[opt] lower border on chr to view, in MB. [default: 0]")
)

parser <- OptionParser(usage = "cgmaptools heatmap [options]",
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
                       Region  Sample1  Sample2 ...  \
                       Region1 0.1      0.1     ...  \
                       Region2 0.1      0.1     ...  \
                       "
)
## check arguments
arguments <- parse_args(parser)
infile <- arguments$infile
if(infile == ""){ # default, STDIN
  infile <- file("stdin")
} else { # user specified
  if( file.access(infile) == -1){ # file not exists
    #print_help(parser)
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
  #outfile <- paste(outfile, figure.format, sep = ".")
}
figure.resolution <- arguments$figure.resolution
low.color <- arguments$low.color
mid.color <- arguments$mid.color
if(mid.color == ""){
  mid.color <- NA
}
high.color <- arguments$high.color
num.color <- arguments$num.color

centr <- data.frame(chr=c("chr1A","chr1B","chr1D","chr2A","chr2B","chr2D","chr3A","chr3B","chr3D","chr4A","chr4B","chr4D","chr5A","chr5B","chr5D","chr6A","chr6B","chr6D","chr7A","chr7B","chr7D"),
                    len=c(215,239,173,337,346,268,300,346,241,300,317,184,254,206,189,286,327,211,357,287,340))

###
meta_in <- arguments$meta_in
meta_data <- read_delim(meta_in, "\t", escape_double = FALSE, trim_ws = TRUE, comment = "#")
data <- read_delim(paste0(infile, ".100k.join"), "\t", escape_double = FALSE, trim_ws = TRUE)
#row.names(data) <- data$loc
#data[,1] <- NULL

upper <- arguments$upper
lower <- arguments$lower
data <- data[(upper*10 + 1):(lower*10 + 1),]

## filter sample
tmp_index <- colnames(data) %in% meta_data$Sample
data <- data[,tmp_index]

sample_present <- meta_data[meta_data$Sample %in% colnames(data),]
LIST <- sample_present$label
names(LIST) <- sample_present$Sample

colnames(data) = LIST[colnames(data)]

## label color
colors = c("#5fb951", "#E74D4A", "#D9B460", "#AB6193","#71ACDF", "#F194BE", "#424E9F")
names(colors) = unique(sample_present$Group)

## figure device
if (figure.format == "png"){
  png(paste0(outfile, ".png"), type="cairo", height = figure.height, width = figure.width, res = figure.resolution, units = "in")
} else if (figure.format == "pdf"){
  pdf(paste(outfile,upper,"MB",lower,"MB","DPheatmap.pdf", sep = "_"), height = figure.height, width = figure.width)
}

## cluster specified
cluster <- arguments$cluster
cluster=T
if (cluster){
  ## figure margin 
  par(mar=c(8, 0, 8, 0), oma=c(0,4,2,3))
  layout(matrix(1:2, nrow=1), widths=c(1,6))

  #correlation <- cor(data, method="spearman", use="pairwise.complete.obs")
  #hc <- hclust(d = as.dist((1-correlation)^2))

  hc <- hclust(d = dist(t(data), method = "euclidean"))
  
  dendrogram <- as.dendrogram(hc)
  plot(dendrogram, horiz=T, leaflab="none", yaxt = "none", dLeaf = 0, yaxs="i")
  
  # reorder sample index
  data <- data[, order.dendrogram(dendrogram)]
} else {
  ## figure margin 
  par(mar=c(3, 6, 4.1, 2.1))
}

## empty figure panel
plot(x=0, type="n", bty="n", xaxt="n", yaxt="n", 
     xlab="", ylab="", 
     xlim=c(0, nrow(data) * 1.15), ylim=c(1, ncol(data)+1), xaxs="i", yaxs="i")

## group info color
colors_info <- c("green", "red")
names(colors_info) <- c("Y", "N")

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
  rect(xleft = 1:nrow(data) + 0.09 * nrow(data), 
       ybottom = rep(y.i, nrow(data)) + 0.1, 
       xright = 1:nrow(data) + 0.09 * nrow(data) + 1, 
       ytop = rep(y.i, nrow(data)) + 0.1 + row.height,
       col = x.color,
       border = x.color
  )
  
  # methylation level on the right panel, barplot
  met.mean <- mean(as.numeric(data[,i]), na.rm = T)
  rect(xleft = nrow(data) + 0.09 * nrow(data) + 1,
       ybottom = y.i + 0.1,
       xright = met.mean * nrow(data)*(rightPanel.extension - 1.12) + nrow(data) + 0.09 * nrow(data) + 1,
       ytop = y.i + 0.1 + row.height,
       col = "#66CD00",
       border = "#66CD00"
  )
  text(x= met.mean * nrow(data)*(rightPanel.extension - 1.12) + nrow(data) + 0.09 * nrow(data) + 2,
       y=y.i+0.1+row.height/2,
       paste(format(met.mean*100,digits = 3), "%", sep = ""),
       cex = 0.7, adj = c(0, 0.5)
  )
}

## sample label
sampleLabel.col <- "black"
if (cluster){
  text(x=0, y=1:ncol(data)+row.height/2 + 0.1, names(data), cex=1, xpd=T, adj=c(0,0.5), col=colors[sample_present$Group[match(colnames(data), sample_present$label)]])
} else {
  text(x=1, y=1:ncol(data)+row.height/2 + 0.1, names(data), cex=0.8, xpd=T, adj=c(1,0.5))
}

## group info by points
points(x=rep(0.09 * nrow(data) - 2, ncol(data)), y=1:ncol(data)+row.height/2 + 0.1,
       col=colors_info[sample_present$Tip_group[match(colnames(data), sample_present$label)]],
       xpd = T,
       cex = 2,
       pch = 16)

## methylation value legend
colorLegend.x.min <- nrow(data)*0.1
colorLegend.x.max <- nrow(data)*0.2
colorLegend.x <- seq(from=colorLegend.x.min, to=colorLegend.x.max, length.out = num.color+1)
colorLegend.y <- -4
rect(xleft = colorLegend.x[-length(colorLegend.x)],
     xright = colorLegend.x[-1],
     ybottom = colorLegend.y, 
     ytop = -2,
     col = color.rects,
     border = color.rects,
     xpd = T
)

# value label
text(x=seq(from=colorLegend.x.min, to=colorLegend.x.max, length.out = 3),
     y=rep(-2, 3),
     seq(0,2,length.out = 3),
     adj = c(0.5, 1),
     cex= 4,
     xpd = T
)

# header
text(x = nrow(data) * 0.6,
     y = ncol(data) + 0.1 + row.height + 2,
     paste(outfile,upper,"MB to",lower,"MB","DPheatmap", sep = " "),
     adj = c(0.5, 0),
     cex = 5,
     xpd = T
)

axis.y <- 0.2
nticks <- ceiling(nrow(data)/100)

# axis
segments(x0 = 1 + 0.09 * nrow(data), 
         y0 = axis.y, 
         x1 = nrow(data) + 0.09 * nrow(data) + 1, 
         y1 = axis.y,
         col = "black",
         xpd = T
)

# axis ticks
segments(x0 = c(seq(from=0, to=(nticks-1) * 100, length.out = nticks)) + 1 + 0.09 * nrow(data),
         x1 = c(seq(from=0, to=(nticks-1) * 100, length.out = nticks)) + 1 + 0.09 * nrow(data),
         y0 = rep(axis.y, nticks),
         y1 = rep(axis.y - 0.1, nticks),
         xpd = T
)

# axis text
text(x = seq(from=0, to=(nticks-1) * 100, length.out = nticks) + 1 + 0.09 * nrow(data),
     y = rep(axis.y - 0.8, nticks),
     c(seq(from=upper, to=upper+(nticks-1) * 10, length.out = nticks)),
     cex = 2,
     xpd = T
)

# centro annotation
if (centr$len[match(infile, centr$chr)] > upper & centr$len[match(infile, centr$chr)] < lower){
  points((centr$len[match(infile, centr$chr)] - upper) * 10 + 1 + 0.09 * nrow(data), axis.y-0.5, pch=2, col = "green", xpd = T, cex = 2.5)
  text(x = (centr$len[match(infile, centr$chr)] - upper) * 10 + 1 + 0.09 * nrow(data),
     y = axis.y-1.2,
     "Centromere",
     cex = 2,
     xpd = T
  )
}
invisible(dev.off())
