#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2018-08-08

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mixtools"))
option_list <- list(
  make_option(c("-p", "--pathway"), dest = "pathway", default = "",
              help="output path"),  
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="input file"),
  make_option(c("-o", "--outgraph"), dest = "outgraph", default = "",
              help = "[opt] output graph"),
  make_option(c("-e", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file"),
  make_option(c("-t", "--title"), dest = "title", default = "", 
              help = "[opt] text shown on title"),
  make_option(c("-f", "--format"), dest = "format", default = "pdf",
              help = "[opt] the format for output figure: pdf (default), png, eps"),
  make_option(c("-W", "--width"), dest = "fig_width", default = 6,
              help = "[opt] width (in inch). Default: %default."),
  make_option(c("-H", "--height"), dest = "fig_height", default = 6,
              help = "[opt] height (in inch). Default: %default.")
)
#
parser <- OptionParser(usage = "%prog [options] file", 
   option_list=option_list, description = "Author: Guo, Weilong; guoweilong@126.com; 2018-08-08 \
Description: Draw the Boxplot chart for the five quantile number\
./DP_hist.R -i  /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/DP/s1_DP/chr1A.DP -p /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/ -H 6 -W 12 -f pdf -t Density plot of S1_1A deepth -o S1_1A deepth.pdf -e S1_1A.txt \
Input examples: \
  sample min q1 median q3 max\
  A	    2    4   6    8   12\
  B	    5    6   7    8   11"
)
#
arguments <- parse_args(parser)
infile <- arguments$infile
figure.format = arguments$format
#
outgraph = paste(arguments$pathway, arguments$outgraph, sep="")
outfile = paste(arguments$pathway, arguments$outfile, sep="")

if( outgraph == "") {
  outgraph = paste( "output", "png", sep=".")
}
#
#infile = "Yellow_module_David_ZhuList_simple.txt"
if(infile == "") {
  #T = read.table( stdin(), header=TRUE, check.names = FALSE, sep = "\t", comment.char = "")
  T = read.table( file("stdin"), header=FALSE, check.names = FALSE, sep = "\t", comment.char = "")
} else {
  T = read.table(infile, header=FALSE, check.names = FALSE, sep = "\t", comment.char = "")
}
mydata=T[complete.cases(T),]
#png(outfile, width = arguments$png_width, height = arguments$png_height)
figure.width = as.integer(arguments$fig_width)
figure.height = as.integer(arguments$fig_height)
options(scipen = 200)
a  <- mydata[,3] # remove NA
  #i = 7
if (figure.format == "pdf") {
  pdf(outgraph, width = figure.width, height = figure.height )
} else if (figure.format == "eps") {
  postscript(outgraph, width = figure.width, height = figure.height)
} else { # png
  png(outgraph, width = figure.width*72, height = figure.height*72 )
}
A <- density(a)
mixmdla <- normalmixEM(a,k=2)
R <- qnorm(0.995,mixmdla$mu[1], mixmdla$sigma[1])
L <- qnorm(0.005,mixmdla$mu[1], mixmdla$sigma[1])
#while (R > 11 | L < 8.5)
#while (R > 9 | L < 5)
#{
#   mixmdla <- normalmixEM(a,k=2)
#   R <- qnorm(0.995,mixmdla$mu[1], mixmdla$sigma[1])
#   L <- qnorm(0.005,mixmdla$mu[1], mixmdla$sigma[1])
#}
plot.new()
plot.window(xlim=c(min(A$x),12), ylim=c(0,max(A$y)+1))
plot(mixmdla, which =2, 
     breaks=100,
     xlab2 = "SNP Average Deepth",ylab2 = "Density",main2 = arguments$title,
     add = TRUE)
axis(1); axis(2); title(main=arguments$title, 
                        xlab="SNP Average Deepth", ylab="Density",
                        cex.main = 3,cex.lab=1.5)
print(outfile)
print(outgraph)
abline(v=R)
abline(v=L)
text(x = R+0.5, y = 0.8 , labels = paste("Confidence level(99.5%)=",format(R , digits = 3),sep=""))
text(x = R+0.5, y = 0.6 , labels = paste("Confidence level(0.5%)=",format(L , digits = 3),sep=""))
dev.off()
cat( L,R, file=outfile, sep = "\t")
#

