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
  make_option(c("-H","--height"), dest = "figure.height", default = 50,
              help = "[opt] height of figure (inch). [default: 50]"),
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

infile <- '/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/combine.1M.join_norm'
meta_in  <- '/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/metadata_cultivar_nogroup.txt'
colorF <- '/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/sub_genome_tree/coloryzz_6cols.txt'

colorFile <- read.table(colorF, header = T, stringsAsFactors = F, sep = "\t", comment.char = "")
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
colors <- colorFile[,2]
names(colors) <- colorFile[,1]

## figure device
if (figure.format == "png"){
  png(paste0(outfile, ".png"), type="cairo", height = figure.height, width = figure.width, res = figure.resolution, units = "in")
} else if (figure.format == "pdf"){
  pdf(outfile, height = figure.height, width = figure.width)
}



  ## figure margin 
par(mar=c(3,5 , 3, 5), oma=c(0,0,0,0))
#layout(matrix(1:2, nrow=1), widths=c(1,6))
#correlation <- cor(data, method="spearman", use="pairwise.complete.obs")

dt <-dist(t(data),method = "euclidean")
NJ_use <- njs(dt)
DF_CNV <- full_join(as_tibble(NJ_use), sample_present, by = 'label')

hc <- hclust(dt)
dendrogram <- as.dendrogram(hc)
phylo <- as.phylo(hc)
col.vector <- vector(mode="character",length=nrow(phylo$edge))
n.tips <- length(phylo$tip.label)
col.vector[phylo$edge[,2]>n.tips] <- "black"
edge.data <- as.data.frame(phylo$edge)
for(i in seq_along(phylo$tip.label)){
  edge.row <- as.numeric(rownames(edge.data[edge.data$V2==i,]))
  col.vector[edge.row] <- colors[DF_CNV$Group[match(phylo$tip.label[i], DF_CNV$label)]]
}
tip.color <- colors[DF_CNV$Group]
names(tip.color) <- as.character(DF_CNV$label)
plot(phylo, type = 'u',show.tip.label = T,
      edge.color = col.vector,cex = 0.2)



col.vector <- vector(mode="character",length=nrow(NJ_use$edge))
n.tips <- length(NJ_use$tip.label)
col.vector[NJ_use$edge[,2]>n.tips] <- "black"
edge.data <- as.data.frame(NJ_use$edge)
for(i in seq_along(NJ_use$tip.label)){
  edge.row <- as.numeric(rownames(edge.data[edge.data$V2==i,]))
  col.vector[edge.row] <- colors[DF_CNV$Group[match(NJ_use$tip.label[i], DF_CNV$label)]]
}
tip.color <- colors[DF_CNV$Group]
names(tip.color) <- as.character(DF_CNV$label)
plot(NJ_use, type = "u", show.tip.label = T, edge.color = col.vector, cex=0.2)
legend("bottom", fill = colors, legend = names(colors), border = F, box.col = "white", horiz=T, cex = 1, x.intersp = 0,xpd=T)

# reorder sample index
data <- data[, order.dendrogram(dendrogram)]

## empty figure panel



invisible(dev.off())
