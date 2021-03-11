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

# Arguments
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help = "input file"),
  make_option(c("-k", "--karyotype"), dest = "karyotype", default = "",
              help = "[opt] 1st-col:chr-name, 2ed-col:centro-loc, 3rd-col:chr-len, used to plot chrome border and centro"),
  make_option(c("-p", "--prefix"), dest = "prefix", default = "",
              help = "[opt] output file name. [default: DPHeatmap+systime]"),
  make_option(c("-s", "--binsize"), dest = "binsize", default = 1,
              help = "[opt] bin size of input file, in MB [default: 1]"),  
  make_option(c("-c","--color"), dest = "color", default = "blue,white,red",
              help = "[opt] three colors needed, in order of low,middle,high [default: blue,white,red]"),
  make_option(c("-n","--colorNumber"), dest = "num.color", default = 10,
              help = "[opt] desired number of color elements in the panel. [default: 10]"),
  make_option(c("-t","--title"), dest = "figure.title", default = "",
              help = "[opt] plot title. [default: DP heatmap of `sample-name`]"),
  make_option(c("-W","--width"), dest = "figure.width", default = 20,
              help = "[opt] width of figure (inch). [default: 7]"),
  make_option(c("-H","--height"), dest = "figure.height", default = 15,
              help = "[opt] height of figure (inch). [default: 7]"),
  make_option(c("-f","--format"), dest = "figure.format", default = "pdf",
              help = "[opt] format of output figure. Alternative: pdf. [default: png]"),
  make_option(c("-R","--resolution"), dest = "figure.resolution", default = 200,
              help = "[opt] Resolution in ppi. Only available for png format. [default: 200]")
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
prefix = arguments$prefix
if(prefix == ""){ # default, "FragRegView.Date"
  prefix = paste("DPHeatmap", Sys.Date(), sep=".")
} else { # user specified
  prefix <- gsub(".pdf$|.png$", "", prefix, perl=T)
}
#
if(infile == ""){ # default, STDIN
  karyotype = F
} else { # user specified
  if( file.access(infile) == -1){ # file not exists
    print_help(parser)
  }
}
#
figure.width <- arguments$figure.width
figure.height <- arguments$figure.height
figure.format <- arguments$figure.format
#
if(! figure.format %in% c("pdf", "png")){ # format not support
  print_help(parser)
} else {
  outfile <- paste(prefix, figure.format, sep = ".")
}
#
figure.resolution <- arguments$figure.resolution
color <- strsplit(arguments$color, split = ",", fixed = T)[[1]]
high.color <- color[3]
mid.color <- color[2]
low.color <- color[1]
#
binsize <- arguments$binsize
num.color <- arguments$num.color
figure.title <- arguments$figure.title
if(figure.title == ""){ # default, STDIN
  figure.title = prefix
} 
#
# define functions
plotK = F
karyotype <- arguments$karyotype
if(karyotype != ""){ # default, STDIN
  if( file.access(infile) == -1){ # file not exists
    print_help(parser)
  }
  else{
    plotK = T
    karyotype <- read.table(karyotype, sep = "\t", header = T, 
                            stringsAsFactors = F, 
                            colClasses = c("character", "integer", "integer"))
    karyotype[[2]] <- ceiling(karyotype[[2]]/1000000)
    karyotype[[3]] <- ceiling(karyotype[[3]]/1000000)
    plotChrom <- function(xleft, ybottom, height, len, centro, binsize){
      
      len <- len/binsize
      centro <- centro/binsize
      
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
      
      # chrom border
      if (is.na(centro)){i
        cbx <- c(xleft + rh - rh*sin(rs), xleft + len - rh + rh*sin(rs))
        cby <- c(ybottom + rv - rv*cos(rs), ybottom + rv + rv*cos(rs))
        polygon(cbx, cby, border = "grey", lwd = 3)
      }else{
        cbx <- c(xleft + rh - rh*sin(rs), centro - 0.5*rh, centro, centro + 0.5*rh, xleft + len - rh + rh*sin(rs), centro + 0.5*rh, centro, centro - 0.5*rh)
        cby <- c(ybottom + rv - rv*cos(rs), ybottom + height, ybottom + rv, ybottom + height, ybottom + rv + rv*cos(rs), ybottom, ybottom + rv, ybottom)
        polygon(cbx, cby, border = "grey", lwd = 3)
      }
    }
  }
}

## figure device
if (figure.format == "png"){
  png(outfile, height = figure.height, width = figure.width, res = figure.resolution, units = "in")
} else if (figure.format == "pdf"){
  pdf(outfile, height = figure.height, width = figure.width)
}

###

data<- read.table(infile, sep = "\t", header = T)

# reverse chr order
data <- data[, ncol(data):1]

## figure margin 
par(mar=c(3, 6, 4.1, 2.1))

## empty figure panel
plot(x=0, type="n", bty="n", xaxt="n", yaxt="n", 
     xlab="", ylab="", 
     xlim=c(0, nrow(data) * 1.15), ylim=c(-1, ncol(data)+1), xaxs="i", yaxs="i")


color.rects <- colorpanel(low=low.color, mid=mid.color, high=high.color, n=num.color)

methylation.interval <- seq(0, 2, length.out=num.color+1)
row.height <- 0.7
rightPanel.extension <- 1.15

for(i in 1:ncol(data)){
  #y.i <- (ncol(data):1)[i]
  y.i <- i
  x.color <- color.rects[ findInterval(data[,i], methylation.interval, all.inside=T)]
  
  # heatmap body
  rect(xleft = 1:nrow(data) + 0.05 * nrow(data), 
       ybottom = rep(y.i, nrow(data)), 
       xright = 1:nrow(data) + 0.05 * nrow(data) + 1, 
       ytop = rep(y.i, nrow(data)) + row.height,
       col = x.color,
       border = x.color)
       
  # find chromo limit
  # limi = match(NA, data[,i], nomatch = nrow(data))
  
  # draw chromo border
  #rect(xleft = 1 + 0.05 * nrow(data), 
  #     ybottom = y.i, 
  #     xright = limi + 0.05 * nrow(data) + 1, 
  #     ytop =  y.i + 0.1 + row.height,
  #     border = "black",
  #     col = "transparent")
  
  # draw centro and border
  if (plotK == T){
    n <- match(names(data)[i], karyotype[[1]])
    if(!is.na(n)){
      len <- karyotype[[2]][n]
      centro <- karyotype[[3]][n]
      plotChrom(0.05 * nrow(data) + 1, y.i, row.height, len, centro, binsize)
    }
  }
}

## sample label
sampleLabel.col <- "black"
text(x= 0.02 * nrow(data), y=1:ncol(data)+row.height/2 + 0.1, names(data), cex=2, xpd=T, adj=c(1,0.5))

## methylation value legend
colorLegend.x.min <- nrow(data)*0.1
colorLegend.x.max <- nrow(data)*0.2
colorLegend.x <- seq(from=colorLegend.x.min, to=colorLegend.x.max, length.out = num.color+1)
colorLegend.y <- 0
rect(xleft = colorLegend.x[-length(colorLegend.x)],
     xright = colorLegend.x[-1],
     ybottom = -0.5, 
     ytop = 0,
     col = color.rects,
     border = color.rects,
     xpd = T
)

# value label
text(x=seq(from=colorLegend.x.min, to=colorLegend.x.max, length.out = 3),
     y=rep(-0.5, 3),
     seq(0,2,length.out = 3),
     adj = c(0.5, 1),
     cex= 2,
     xpd = T
)
# header
text(x = nrow(data) * 0.6,
     y = ncol(data) + 0.1 + row.height + 1,
     figure.title,
     adj = c(0.5, 0),
     cex = 3,
     xpd = T
)

axis.y <- 0.5
# axis
segments(x0 = 1 + 0.05 * nrow(data), 
         y0 = axis.y, 
         x1 = 1 + 1.05 * nrow(data), 
         y1 = axis.y
)

# axis ticks
segments(x0 = seq(from=0, to=800/binsize, length.out = 5) + 1 + 0.05 * nrow(data),
         x1 = seq(from=0, to=800/binsize, length.out = 5) + 1 + 0.05 * nrow(data),
         y0 = rep(axis.y, 5),
         y1 = rep(axis.y - 0.1, 5),
         xpd = T
)

# axis text
text(x = seq(from=0, to=800/binsize, length.out = 5) + 1 + 0.05 * nrow(data),
     y = rep(axis.y - 0.4, 5),
     seq(from=0, to=800, length.out = 5),
     cex = 2,
     xpd = T
)


invisible(dev.off())

# grid
#vps <- baseViewports()
#pushViewport(vps$inner, vps$figure, vps$plot)
#grid.roundrect(x = unit(0.05 * nrow(data), "native"),
#               y = unit(1.1, "native"),
#               height =  unit(0.8, "native"),
#               width = unit(nrow(data) + 0.05 * nrow(data),"native"),
#               fill="transparent")
