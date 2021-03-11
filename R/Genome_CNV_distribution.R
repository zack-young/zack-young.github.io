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
infile <- strsplit((arguments$infile), split=",")[[1]]
karyotype_use <- arguments$karyotype

a <- '/data2/rawdata2/PAV/Tibet_population/Core_PAV_HA_vs_LA'
infile <- c('/data2/rawdata2/PAV/Tibet_population/Core_PAV_HA_vs_LA/LA_abs_gene_rank_count.log.txt',
            '/data2/rawdata2/PAV/Tibet_population/Core_PAV_HA_vs_LA/HA_abs_gene_rank_count.log.txt',
            '/data2/rawdata2/PAV/Tibet_population/Core_PAV_HA_vs_LA/core_PAV_gene_count.norm.txt')

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
chr_max_length <- max(karyotype$Length)/1000000

color <- c("#4889E6", "#C15065", "red")
track_height <- 1
gap <- 0.1
chr_height <- track_height*2 + gap*2
total_height <- chr_height*nrow(karyotype)
x_blank <- chr_max_length*0.05

pdf(paste0(prefix, ".pdf"), height = figure.height, width = figure.width)
par(mar=c(2, 2, 0, 2))
pdf("/data/user/yangzz/mapping/fieldergenomecompare/pdf/201_sample_CNV_distri.pdf", height = 13, width = 10)

## empty figure panel
plot(x=0, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="", 
     xlim=c(0, 900+x_blank),
     #xlim=c(0, x_blank+chr_max_length),
     ylim=c(-1, total_height+1), xaxs="i", yaxs="i")

for (i in 1:nrow(karyotype)){
  chrID <- karyotype[i,1]
  anchor_y <- total_height - i*chr_height
  
  #for (j in 1:nfiles){
  #  DF2 <- DF1[[j]][DF1[[j]]$V1==chrID,]
  #  #DF2$V4 <- DF2$V4/
  #  y_pad <- (j-1)*(track_height+track_gap)+anchor_y
  #  polygon(x=x_blank+c(DF2$V2,DF2$V2[length(DF2$V2)],0), y=c(DF2$V4+y_pad,y_pad,y_pad), col = color[j], border = color[j])
  #}
  DF2 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/CNV_count/",chrID,".combine_mask_CNV_deletion_compensent",sep=""),as.is = T, header = F, comment.char = "")
  #/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/field_CNV_count/
  
  # max value
  y_pad <- track_height + gap + anchor_y
  #lis = log(2-DF2$V3,10)
  lis = DF2$V3
  print(max(lis))
  polygon(x=x_blank+c(0,DF2$V1/1000000,DF2$V1[length(DF2$V1)]/1000000), y=c(y_pad,y_pad-lis,y_pad), col = color[1], border = color[1])
  
  DF2 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/CNV_count/",chrID,".combine_mask_CNV_duplication_compensent",sep=""),as.is = T, header = F, comment.char = "")
  lis = DF2$V3
  #lis = log(DF2$V3,10)
  print(max(lis))
  # max value
  
  polygon(x=x_blank+c(0,DF2$V1/1000000,DF2$V1[length(DF2$V1)]/1000000), y=c(y_pad,lis+y_pad,y_pad), col = color[3], border = color[3])
  lines(x=x_blank+c(0,DF2$V1[length(DF2$V1)]/1000000),y=c(y_pad,y_pad))
#  j=3
#  DF2 <- DF1[[j]][DF1[[j]]$V1==chrID,]
#  y_pad <- (j-1)*(track_height+track_gap)+anchor_y
#  polygon(x=x_blank+c(DF2$V2,DF2$V2[length(DF2$V2)],0), y=c(DF2$V4+y_pad,y_pad,y_pad), col = color[j], border = color[j])
  
  # border
  #rect(xleft = x_blank, 
  #     ybottom = anchor_y, 
  #     #xright = x_blank+chr_max_length, 
  #     xright =900+x_blank,
  #     ytop = anchor_y + chr_height)
  #points(karyotype[i,3]/1000000+x_blank, anchor_y-0.3, pch=17, cex=1.3)
  
  text(x_blank, y_pad, chrID, cex=1.5, pos = 2,xpd=T)
}
rect(xleft = x_blank,
    ybottom = 0,
    #xright = x_blank+chr_max_length,
    xright =900+x_blank,
    ytop = total_height)
ticky_x <- c(0, 150,300,450, 600,750,900)
segments(x0 = ticky_x+x_blank,
         x1 = ticky_x+x_blank,
         y0 = 0,
         y1 = -0.5,
         xpd = T)
text(ticky_x+x_blank, -0.1,paste(ticky_x,'Mb',sep = "") , cex=1, pos = 1,xpd=T)
dev.off()

## sample label
sampleLabel.col <- "black"
#ext(x= 0.02 * nrow(data), y=1:ncol(data)+row.height/2 + 0.1, names(data), cex=2, xpd=T, adj=c(1,0.5))

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
