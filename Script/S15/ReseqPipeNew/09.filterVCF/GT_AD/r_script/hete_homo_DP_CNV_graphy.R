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
  A         2    4   6    8   12\
  B         5    6   7    8   11"
)
#
arguments <- parse_args(parser)
infile <- arguments$infile
figure.format = arguments$format
#
outgraph = paste(arguments$pathway, arguments$outgraph, sep="")

if( outgraph == "") {
  outgraph = paste( "output", "png", sep=".")
}
#
#infile = "Yellow_module_David_ZhuList_simple.txt"
if(infile == "") {
  #T = read.table( stdin(), header=TRUE, check.names = FALSE, sep = "\t", comment.char = "")
  mydata = read.table( file("stdin"), header=FALSE,check.names = FALSE, sep = "\t", comment.char = "")
} else {
  mydata = read.table(infile, header=FALSE,check.names = FALSE, sep = "\t", comment.char = "")
}
#print(class(as.character(mydata[1,1])))
#print(files[1,1])
#png(outfile, width = arguments$png_width, height = arguments$png_height)
figure.width = as.integer(arguments$fig_width)
figure.height = as.integer(arguments$fig_height)
options(scipen = 200)
  #i = 7
if (figure.format == "pdf") {
  pdf(outgraph, width = figure.width, height = figure.height )
} else if (figure.format == "eps") {
  postscript(outgraph, width = figure.width, height = figure.height)
} else { # png
  png(outgraph, width = figure.width*72, height = figure.height*72 )
}
options(scipen=200)
plot(x=0, type="n", bty="n", yaxt="n",
        xlab="", ylab="", 
#        xlim=c(-200, 1000), ylim=c(0, nrow(files)+1),
        xlim=c(-200, 1000), ylim=c(0, nrow(mydata)+1),
        xaxs="i", yaxs="i", main=arguments$title,
cex.main = 3)

color_pad <- c("black","gray", "green", "blue2")

for(j in 1:nrow(mydata)){
    if(!is.na(mydata[j,1])){
            data1 <- read.table(as.character(mydata[j,1]), as.is = T, header = F, comment.char = "")
            chro <- factor(data1[,5])
            chro1 <- factor(chro, levels = c('CNV', 'DP', 'hete', 'homo'), 
                            labels = c('1', '2', '3', '4'))
            for(i in 1:nrow(data1)){
              rect(xleft = as.numeric(data1[i,3])/1000000,
            #         #ybottom = j+0.1,
              ybottom = j+0.1,
              ytop = j+0.9,
              xright = as.numeric(data1[i,4])/1000000,
            #         #ytop = j+0.9,
              col = color_pad[as.numeric(as.vector(chro1[i]))],
              border = color_pad[as.numeric(as.vector(chro1[i]))] 
              )
            }
            #text(x=-100, y=j+0.5, files[j,2], cex = 0.7)
             text(x=-100, y=j+0.5, mydata[j,2], cex = 1)
             text(x=c(950,950,950,950), y=c(20.5,21.5,22.5,23.5), 
                  c("CNV","DP","Hete","Homo"), cex = 0.7)
             rect(xleft = 900, xright = 930, ybottom = 20, ytop = 21,
                  col = color_pad[1], border = color_pad[1])
             rect(xleft = 900, xright = 930, ybottom = 21, ytop = 22,
                  col = color_pad[2], border = color_pad[2])
             rect(xleft = 900, xright = 930, ybottom = 22, ytop = 23,
                  col = color_pad[3], border = color_pad[3])
             rect(xleft = 900, xright = 930, ybottom = 23, ytop = 24,
                  col = color_pad[4], border = color_pad[4])

    }
}



dev.off()

