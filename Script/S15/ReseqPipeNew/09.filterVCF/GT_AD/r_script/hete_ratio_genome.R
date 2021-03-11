#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2018-08-08

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mixtools"))
option_list <- list(
  make_option(c("-p", "--pathway"), dest = "pathway", default = "",
              help="output path"),
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="input file"),
  make_option(c("-O", "--outFDR"), dest = "outFDR", default = "",
              help = "[opt] output FDR"),
  make_option(c("-o", "--outgraph"), dest = "outgraph", default = "",
              help = "[opt] output graph"),
  make_option(c("-e", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file"),
  make_option(c("-t", "--title"), dest = "title", default = "",
              help = "[opt] text shown on title"),
  make_option(c("-T", "--Title"), dest = "Title", default = "",
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
  A         2    4   6    8   12\
  B         5    6   7    8   11"
)
#
arguments <- parse_args(parser)
infile <- arguments$infile
figure.format = arguments$format
#
outFDR = paste(arguments$pathway, arguments$outFDR, sep="")
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
options(scipen = 200)
Adata <- mydata[which(mydata[,7]<0.04),]
#Adata <- mydata
A <- density(Adata[,7])
LENA <- length(Adata[,7])
SDA <- sd(Adata[,7])
MEA <- mean(Adata[,7])
DTA <- rnorm(LENA, mean = MEA, sd = SDA)
A2 <- density(DTA)
#
figure.width = as.integer(arguments$fig_width)
figure.height = as.integer(arguments$fig_height)
  #i = 7
if (figure.format == "pdf") {
  pdf(outFDR, width = figure.width, height = figure.height )
} else if (figure.format == "eps") {
  postscript(outFDR, width = figure.width, height = figure.height)
} else { # png
  png(outFDR, width = figure.width*72, height = figure.height*72 )
}
#
data_list = list() 
loop = min(which(A2$x>0.00))
loop2 = min(which(A$x>0.00))
del = abs(loop - loop2)
for (j in seq(min(which(A2$x>0.00)),(length(A$x))-2,1)) {
  b = 0
  b2 = 0
  loop = loop + 1
  loop2 = loop2+1
  for (i in seq(1,length(A2$x)-loop,1)){
    #i = 400
    a = ((max(A$x)-min(A$x))/length(A$x))*A$y[loop2+i]
    b = b + a
    a2 = ((max(A2$x)-min(A2$x))/length(A2$x))*A2$y[loop+i]
    b2 = b2 + a2
    #print(i)
  }
  data_list[[j]] = c(b2/b,A2$x[j])
}
tea <- do.call(rbind, data_list)
#df <- data.frame(ID = matrix(unlist(data_list), byrow=T))
ad <- min(tea[which(tea[,1]>=0.01)]) # FDR = 0.01
ad2 <- tea[tea[,1] >= ad & tea[,1] < ad + 0.00000001,][2]
plot(x = tea[,2],y = tea[,1],
     #ylim = c(0,1),
     type = "l",xlab = "Hete/Homo",ylab = "FDR"
     ,xlim = c(0,max(tea[,2]+0.02))
#     ,ylim = c(0,1)
     ,main = arguments$title)
abline(v=ad2,lwd=1,col="blue")
text(x = ad2,y = 0.2 , pos = 2,
     labels = paste("FDR=0.01 ratio=",format(ad2 , digits = 3),
                    sep=""))
dev.off()
cat( ad2, file=outfile)
#
if (figure.format == "pdf") {
  pdf(outgraph, width = figure.width, height = figure.height )
} else if (figure.format == "eps") {
  postscript(outgraph, width = figure.width, height = figure.height)
} else { # png
  png(outgraph, width = figure.width*72, height = figure.height*72 )
}
plot(A,freq = F ,
     breaks=100,
     xlab = "Hete/Homo",ylab = "Density"
     ,main = arguments$Title,
     col = "red"
     #xlim =c(0,1)
     #,add = TRUE)
)
#hist(DTA,freq = F,breaks=100,add = TRUE)
lines(A2,lty=2,lwd=2, col = "blue")
legend('topright',pch=c(-1,-1),lty=c(1,1),col=c("blue","red"),
       legend = c("simulations", "real data"))
R <- qnorm(0.95,MEA, SDA)
dev.off()
# 画出 density function 
#plot(density(x),breaks=c(0,100),freq=T)
#

