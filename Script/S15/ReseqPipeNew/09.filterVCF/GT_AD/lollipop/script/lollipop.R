#!/usr/bin/env Rscript

# cgmaptools - mCLollipop
#
# Copyright (C) Weilong Guo
# Contact: Weilong Guo <guoweilong@126.com>
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

# Guo, Weilong; guoweilong@126.com; 2015-05-07

is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
if(!is.installed("optparse")){
  warning("Detect package \"optparse\" is not installed in your R enviroment.")
  warning("Trying to install the \"optparse\" package.")
  warning("If failed, please try to install it mannually.")
  install.packages("optparse")
}

# Argument
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "infile3.txt",
              help="input file, use STDIN if ommited, multiple-chr is not suggested"),
  make_option(c("-a", "--annotation"), dest = "annofile", default = "/data2/rawdata2/database/genome_func_anno/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3.gz",
              help="[opt] annotation file name, refFlat format"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "test.pdf",
              help = "[opt] output file"),
  make_option(c("-f", "--format"), dest = "format", default = "pdf",
              help = "[opt] the format for output figure: pdf (default), png, eps"),
  make_option(c("-l", "--left"), dest = "left", default = "65000",
              help = "[opt] Left-most position, use the 1st position if omitted"),
  make_option(c("-r", "--right"), dest = "right", default = "105000",
              help = "[opt] Right-most position, use the last position of input if omitted"),
  make_option(c("-c", "--chr"), dest = "chr", default = "chr1A",
              help = "[opt] chromosome name, use the chr in 1st line of input file if omitted"),
  make_option(c("-t", "--title"), dest = "title", default = "SnpFrep",
              help = "[opt] text shown on title"),
  make_option(c("-p", "--filetag"), dest = "filetag", default = "SnpFrep",
              help = "[opt] text shown on title"),
  make_option(c("-w", "--width"), dest = "fig_width", default = 15,
              help = "[opt] width (in inch). Default: 8."),
  make_option(c("-h", "--height"), dest = "fig_height", default = 10,
              help = "[opt] height (in inch). Default: 8."),
  make_option(c("-n", "--ntrans"), dest = "ntrans", default = 10,
              help = "[opt] transcript number to show. Default: 10")
)
#
parser <- OptionParser(usage = "cgmaptools lollipop [options] file", add_help_option = F,
                       option_list=option_list, description = "      (aka mCLollipop) \
                       Description: Plot local mC level for multiple samples \
                       Contact:     Guo, Weilong; guoweilong@126.com\
                       Last Update: 2018-04-10 \
                       Example: \
                       mCLollipop [-i input] -o gene.png \
                       -Input Format (-i)\
                       Can be output by \"cgmaptools mergelist tomatrix\". Use STDIN if omitted.\
                       The 1st line (header line) is required.\
                       Example: \
                       chr     pos     ratio   severity\
                       Chr1    111403  0.30    nan\
                       Chr1    111406  0.66    0.40\
                       #   chr   left       right      region-description \
                       -annotation file (-a), refFlat Format:\
                       To show the structure of genes/transcripts. One-line in annotation, one-track in figure. \
                       Example: \
                       GeneA   TransA  chr2  +	     1000      2000       1100    1950     3     1100,1500,1700,  1200,1580,1950,\
                       #   GeneID  TrandID ChrID Strand TransLeft TransRight CDSLeft CDSRight nExon ExonLefts        ExonRights\
                       "
)
#
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}
#
gffRead <- function(gffFile, nrows = -1) {
  
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}
#
arguments <- parse_args(parser)
infile <- arguments$infile
outfile = arguments$outfile
annofile = arguments$annofile
plot_title <- arguments$title
chr=arguments$chr
#
if(infile == "") { # default, STDIN
  infile = file("stdin")
  #print_help(parser)
  #stop(sprintf("input file is not specified"))
}  else {
  if(file.access(infile) == -1) {
    stop(sprintf("Specified file ( %s ) does not exist", infile))
  }
}
#
# Read the input file
mainDF <- read.table(infile, header=T, stringsAsFactors = F, sep="\t")
n_group <- ncol(mainDF)-3
group_names <- colnames(mainDF)[-c(1,2,3)]
#
pos <- mainDF[,2]
if( arguments$left=="" || arguments$right=="" ) {
  LeftMost = min(pos)
  RightMost = max(pos)
} else {
  LeftMost = as.integer(arguments$left)
  RightMost = as.integer(arguments$right)
}
#
if(LeftMost==RightMost){
  LeftMost = LeftMost-5000
  RightMost = RightMost+5000
}
#
# Read the annotation file
# tabix iwgsc_refseqv1.0_HighConf_2017Mar13.gff3.gz chr7A:1-100000
tmpannoDF <- system2("/home/wangzh/bin/htslib/bin/tabix", args = c(annofile, paste0(chr, ":", LeftMost, "-", RightMost)), stdout = T)
annoDF <- data.frame(do.call(rbind, strsplit(tmpannoDF, "\t", fixed=TRUE)), stringsAsFactors=F)
annoDF[,4] <- as.numeric(annoDF[,4])
annoDF[,5] <- as.numeric(annoDF[,5])
colnames(annoDF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
annoDF$ID <- getAttributeField(annoDF$attributes, "ID")
annoDF$Par <- getAttributeField(annoDF$attributes, "Parent")
N_anno <- sum(annoDF$feature == "mRNA", na.rm = TRUE)

if( outfile == "") {
  outfile = paste(infile, figure.format, sep=".")
}
#
figure.format = arguments$format
figure.width = as.integer(arguments$fig_width)
figure.height = as.integer(arguments$fig_height)
if (figure.format == "pdf") {
  pdf(outfile, width = figure.width, height = figure.height )
} else { # png
  png(outfile, width = figure.width*72, height = figure.height*72 )
}
#
# ========================
# plot region definition
# plot1: main plot; plot2: header and legend
layout(mat = matrix(c(2,1,3),ncol = 1), heights = c(1,5,1))
# xlim = c(LeftMost, RightMost)
# ylim = c(1,10)
#
color = c("#00D9FFFF", "#00FF19FF", "#FF0000FF", "#FFB300FF","#0026FFFF", "gray")
names(color) <- c("missense_variant", "synonymous_variant", "frameshift_variant", "stop_gained/stop_lost", "splice_region_variant", "others")
#
# ========================
# plot1: main plot
#
par(mar=c(0,6,0,3))
xlength = RightMost - LeftMost
# Set the xlim and ylim
ntrans <- arguments$ntrans
if(ntrans < N_anno) {N_anno <- ntrans}
plot(NULL, NULL, 
     xlim = c(LeftMost, RightMost), 
     ylim = c(1 - N_anno*1.5 - 1, 6 + (n_group-1)*6),
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n", frame.plot = FALSE,
     cex.axis = 1.5, cex = 1.5)
# -------------------------
# Draw the DNA strands
#abline( h = 1, lwd = 2 )
# plot axis
# xaxs does not work as expected
axis(1, pos = 1, lwd = 4, cex.axis = 2, xaxs = "i")
for (j in seq_len(n_group)){
  abline( h = 1 + (j-1)*6, lwd = 4 )
  axis(2, at = 1:5 + (j-1)*6, labels = c("0%", "25%", "50%", "75%", "100%"), las=2, cex.axis = 2, lwd = 4)
  text(x = LeftMost, y = j*6, group_names[j], pos = 1, cex=2.5)
  # Draw the loll
  for( i in seq_len(nrow(mainDF))){
    # lines( c(pos[i], pos[i]), c(0, mC[i,k])+k*yscale , lwd=5)
    segments( mainDF[i,2], 1 + (j-1)*6, mainDF[i,2], 1 + mainDF[i,3+j]*4 + (j-1)*6, lwd = 3, col = color[mainDF[i,3]] )
    points( mainDF[i,2], 1 + mainDF[i,3+j]*4 + (j-1)*6, cex = 2, pch = 16, col = color[mainDF[i,3]] )
  }
}
#
# ==========================
# Draw the gene annotation
if(N_anno > 0){
  Ttrans <- annoDF[annoDF$feature == "mRNA", "ID"]
  for (i in c(1:N_anno)){
    trans <- Ttrans[i]
    annoDFtmp <- rbind(annoDF[annoDF$feature == "mRNA" & annoDF$ID == trans,], 
                       annoDF[annoDF$Par == trans, ])
    annoDFtmp <- annoDFtmp[!is.na(annoDFtmp$seqname),]
    annoDFtmp <- annoDFtmp[order(annoDFtmp$feature, decreasing = T),]

    y_center = 0.5 - N_anno*1.5 + i*1.5 - 1.5
    # gene name on upmost
    for (j in seq_len(nrow(annoDFtmp)) ) {
      annotype = annoDFtmp[j,3]
      Left = annoDFtmp[j,4]
      Right = annoDFtmp[j,5]
      x_center = Left + (Right-Left)/2
      
      if (annotype == "exon"){
        rect(Left, y_center-0.2, Right, y_center+0.2, col = "white")
      } else if (annotype == "CDS"){
        rect(Left, y_center-0.2, Right, y_center+0.2, col = "black")
      } else if (annotype == "mRNA"){
        strand = annoDFtmp[j,7]
        StrandSign=(c(1,-1)[c("+","-")==strand])
        x0 = c(Left, Right)[c("+","-")==strand]
        x2 = c(Right, Left)[c("+","-")==strand]
        x1 = x2 + (RightMost-LeftMost)/30*StrandSign
        segments(Left, y_center, Right, y_center, col = "black")
        arrows(x0, y_center, x1, y_center, length = 0.08 )
        text(Left, y_center+0.4, labels = annoDFtmp[j,10], cex = 1.5, pos = 4)
      }
    }
  }
}
#
# ========================
# plot2: header and legend
#
par( mar = c(0,3,3,3) )
plot(x=0, type="n", bty="n", xaxt="n", yaxt="n", 
     xlab="", ylab="", 
     xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i")

# Text ex:  "chr2: 90,936-91,653"
# text( 0.5, 1, cex = 3,
#      paste(chr, " : ", 
#            format(LeftMost, big.mark = ",", scientific = F), " - ", 
#            format(RightMost, big.mark = ",", scientific = F), sep = "") )
text(0.5, 1, cex = 3, plot_title, pos = 1)

# Draw the color legends
legend("bottomleft", fill = color, legend = names(color), border = F, box.col = "white", horiz=T, cex = 2, x.intersp = 0.5)

#
# ========================
# plot3: tag and paras
#
#filenametag <- paste0("SnpHub SnpFreq ", Sys.time())
filenametag <- arguments$filetag
#
par( mar = c(3,3,0,3) )
plot(NULL, NULL, type="n", bty="n", xaxt="n", yaxt="n", 
     xlab="", ylab="", 
     xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i")
text(x = 1, y = 0.5, paste0(plot_title,"\n",filenametag), pos = 2, cex=2)


dev.off()
#
