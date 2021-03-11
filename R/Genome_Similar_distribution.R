#!/usr/bin/env Rscript
## libraries
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "LA_abs_gene_rank_count.log.txt,HA_abs_gene_rank_count.log.txt,core_PAV_gene_count.norm.txt",
              help = "input file"),
  make_option(c("-k", "--karyotype"), dest = "karyotype", default = "/data2/rawdata2/readDepth/DP_bySampel/CS_karyotype.txt",
              help = "[opt] 1st-col:chr-name, 2ed-col:centro-loc, 3rd-col:chr-len, used to plot chrome border and centro"),
  make_option(c("-p", "--prefix"), dest = "prefix", default = "test",
              help = "[opt] output file name. [default: DPHeatmap+systime]"),
  make_option(c("-t","--title"), dest = "figure.title", default = "",
              help = "[opt] plot title. [default: DP heatmap of `sample-name`]"),
  make_option(c("-W","--width"), dest = "figure.width", default = 20,
              help = "[opt] width of figure (inch). [default: 7]"),
  make_option(c("-H","--height"), dest = "figure.height", default = 15,
              help = "[opt] height of figure (inch). [default: 7]")
)

parser <- OptionParser(usage = "ABD subgenome Similar distribution [options]",
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
#arguments <- parse_args(parser)
#infile <- strsplit((arguments$infile), split=",")[[1]]
#karyotype_use <- arguments$karyotype
karyotype_use <-"/data2/rawdata2/readDepth/DP_bySampel/CS_karyotype.txt"
infile <- strsplit((infile), split=",")#[[1]]
nfiles <- length(infile)
prefix = arguments$prefix
if(prefix == ""){ # default, "FragRegView.Date"
  prefix = paste("DPHeatmap", Sys.Date(), sep=".")
} else { # user specified
  prefix <- gsub(".pdf$|.png$", "", prefix, perl=T)
}
#
#
figure.width <- arguments$figure.width
figure.height <- arguments$figure.height
#
figure.title <- arguments$figure.title
if(figure.title == ""){ # default, STDIN
  figure.title = prefix
} 
#


karyotype <- read.table(karyotype_use, header = T, stringsAsFactors = F)

A <-  karyotype[grep(".*?(A).*?",karyotype$Chr),]  # 提取name中包含字符 a 的行
B <-  karyotype[grep(".*?(B).*?",karyotype$Chr),]
D <-  karyotype[grep(".*?(D).*?",karyotype$Chr),]
chr_max_length <- max(karyotype$Length)/1000000

color <- c("#005DCC", "#009E00", "red")
track_height <- 0.2
chr_height <- track_height
total_height <- chr_height*nrow(karyotype)/3
x_blank <- chr_max_length*0.05

pdf(paste0(prefix, ".pdf"), height = figure.height, width = figure.width)
par(mar=c(2, 2, 0, 2))
## empty figure panel
plot(x=0, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="", 
     xlim=c(0, (980+x_blank)*3),
     #xlim=c(0, x_blank+chr_max_length),
     ylim=c(-0.2, total_height+0.2), xaxs="i", yaxs="i")

for (i in 1:nrow(A)){
  chrID <- A[i,1]
  anchor_y <- total_height - i*chr_height
  centro <- x_blank+A[i,3]/1000000
  DF2 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/similar_distri/",chrID,".CN_combine_hmmed_level_count_compensent",sep=""),as.is = T, header = F, comment.char = "")
  # max value
  y_pad <- anchor_y
  lis = DF2$V3#log(DF2$V3+1,10)
  print(max(lis))
  polygon(x=x_blank+c(0,DF2$V1/1000000,DF2$V1[length(DF2$V1)]/1000000), y=c(y_pad,y_pad+lis,y_pad), col = color[1], border = color[1])
  # border
  # rect(xleft = x_blank, 
  #      ybottom = anchor_y, 
  #      #xright = x_blank+chr_max_length, 
  #      xright =850+x_blank,
  #      ytop = anchor_y + chr_height)
  #points(karyotype[i,3]/1000000+x_blank, anchor_y-0.3, pch=17, cex=1.3)
  tx <- c(centro+10,centro,centro-10)
  ty <- c(y_pad-0.02,y_pad,y_pad-0.02)
  polygon( tx, ty, border = "black", col = "black")
  text(x_blank+30, y_pad, substr(chrID, 4,5), cex=1.5, pos = 2,xpd=T)
}
rect(xleft = x_blank,
     ybottom = 0,
     #xright = x_blank+chr_max_length,
     xright =850+x_blank,
     ytop = total_height)
ticky_x <- c(0,200,400,600,850)
segments(x0 = ticky_x+x_blank,
         x1 = ticky_x+x_blank,
         y0 = 0,
         y1 = -0.02,
         xpd = T)
text(ticky_x+x_blank, -0.02,paste(ticky_x,'Mb',sep = "") , cex=0.75, pos = 1,xpd=T)


x_blankB <- x_blank+ 900+130
for (i in 1:nrow(B)){
  chrID <- B[i,1]
  anchor_y <- total_height - i*chr_height
  centro <- x_blankB+B[i,3]/1000000
  DF2 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/similar_distri/",chrID,".CN_combine_hmmed_level_count_compensent",sep=""),as.is = T, header = F, comment.char = "")
  # max value
  y_pad <- anchor_y
  lis = DF2$V3#log(DF2$V3+1,10)
  print(max(lis))
  polygon(x=x_blankB+c(0,DF2$V1/1000000,DF2$V1[length(DF2$V1)]/1000000), y=c(y_pad,y_pad+lis,y_pad), col = color[1], border = color[1])
  # border
  #points(karyotype[i,3]/1000000+x_blank, anchor_y-0.3, pch=17, cex=1.3)
  tx <- c(centro+10,centro,centro-10)
  ty <- c(y_pad-0.02,y_pad,y_pad-0.02)
  polygon( tx, ty, border = "black", col = "black")
  text(x_blankB+30, y_pad, substr(chrID, 4,5), cex=1.5, pos = 2,xpd=T)
}
rect(xleft = x_blankB, 
     ybottom = 0, 
     #xright = x_blank+chr_max_length, 
     xright =850+x_blankB,
     ytop = total_height  )

segments(x0 = ticky_x+x_blankB,
         x1 = ticky_x+x_blankB,
         y0 = 0,
         y1 = -0.02,
         xpd = T)
text(ticky_x+x_blankB, -0.02,paste(ticky_x,'Mb',sep = "") , cex=0.75, pos = 1,xpd=T)
## sample label
sampleLabel.col <- "black"

x_blankD <- x_blankB+ 900+130
for (i in 1:nrow(D)){
  chrID <- D[i,1]
  anchor_y <- total_height - i*chr_height
  centro <- x_blankD+D[i,3]/1000000
  DF2 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/similar_distri/",chrID,".combine_hmm_level_count_compensent",sep=""),as.is = T, header = F, comment.char = "")#.CN_combine_hmmed_level_count_compensent
  DF2$V3[which(DF2$V3>0.3)] <-0.3
  # max value
  y_pad <- anchor_y
  lis = DF2$V3#log(DF2$V3+1,10)
  print(max(lis))
  polygon(x=x_blankD+c(0,DF2$V1/1000000,DF2$V1[length(DF2$V1)]/1000000), y=c(y_pad,y_pad+lis,y_pad), col = color[1], border = color[1])
  # border
  #points(karyotype[i,3]/1000000+x_blank, anchor_y-0.3, pch=17, cex=1.3)
  tx <- c(centro+10,centro,centro-10)
  ty <- c(y_pad-0.02,y_pad,y_pad-0.02)
  polygon( tx, ty, border = "black", col = "black")
  text(x_blankD+30, y_pad, substr(chrID, 4,5), cex=1.5, pos = 2,xpd=T)
}

rect(xleft = x_blankD, 
     ybottom = 0, 
     #xright = x_blank+chr_max_length, 
     xright =900+x_blankD,
     ytop = total_height)
segments(x0 = ticky_x+x_blankD,
         x1 = ticky_x+x_blankD,
         y0 = 0,
         y1 = -0.02,
         xpd = T)
text(ticky_x+x_blankD, -0.02,paste(ticky_x,'Mb',sep = "") , cex=0.75, pos = 1,xpd=T)
## sample label
sampleLabel.col <- "black"




# header
if(F){
text(x = nrow(data) * 0.6,
     y = ncol(data) + 0.1 + row.height + 1,
     figure.title,
     adj = c(0.5, 0),
     cex = 3,
     xpd = T
)
}
invisible(dev.off())
