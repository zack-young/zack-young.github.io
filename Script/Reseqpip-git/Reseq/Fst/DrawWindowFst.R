#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2017-11-14
#
#
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-c", "--chr"), dest = "chr", default = "",
              help = "The chr name"),
  make_option(c("-t", "--threshold"), dest = "thr", default = 0,
              help = "The threshold value to draw the horizonal line")
)
#
parser <- OptionParser(
    usage = "%prog -c chr1",
    option_list=option_list,
    description = 
      "Author: Guo, Weilong; guoweilong@126.com; 2018-07-12 \
      Description: Draw the Fst distribution"
)
#
arguments <- parse_args(parser)
#
CHR <- arguments$chr
THR <- arguments$thr
#
pdf( paste0( CHR, ".windowed.weir.fst.pdf"), width = 8, height = 4)
par( mai = c(0.5,1,0.5,0.2) )
TAB = read.table( paste0( CHR, ".windowed.weir.fst"), header = TRUE)
Data = TAB[,c(2,5)]
plot(Data[,1]/1000000, Data[,2],
     pch=16, cex = 0.2, 
     col = "grey2",
     main = CHR, cex.main = 1.5,
     ylim = c(-0.1, 1),
     xlab = "(M bp)",
     ylab=expression(F[ST]) , cex.lab = 1.5
)
DataRmNA = Data[!is.na(Data[,2]),]
lines(smooth.spline(DataRmNA[,1]/1000000, DataRmNA[,2], spar=0.35), col='blue', lwd=1)
abline( h=THR, col='red', lwd=1.5)
dev.off()
#

