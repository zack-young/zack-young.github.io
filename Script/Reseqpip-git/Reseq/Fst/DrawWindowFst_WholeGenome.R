#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2017-11-14
#
#
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-c", "--chr"), dest = "chr", default = "",
              help = "The chr name"),
  make_option(c("-w", "--window_size"), dest = "wSize", default = 200000,
              help = "The window size [default: %default]"),
  make_option(c("-t", "--threshold"), dest = "thr", default = 0,
              help = "The threshold for drawing the horizontal line  [default: %default]")
)
#
parser <- OptionParser(
    usage = "%prog -c chr1",
    option_list=option_list,
    description = 
      "Author: Guo, Weilong; guoweilong@126.com; 2018-08-12 \
      Description: Draw the Fst distribution for whole genome"
)
#
arguments <- parse_args(parser)
#
CHR <- arguments$chr
wSize <- as.integer(arguments$wSize)
THR <- as.double(arguments$thr)
#CHR = "chrA.HA_vs_LA.200k"
#
pdf( paste0( CHR, ".windowed.weir.fst.pdf"), width = 20, height = 4)
par( mai = c(0.5,1,0.5,0.2) )
TAB = read.table( paste0( CHR, ".windowed.weir.fst"), header = TRUE)
Data = TAB[,c(1,2,3,5)]
ChrMaxLen = tapply(Data[,3], Data[,1], max)
ChrLenCumsum = cumsum( c(0, t(ChrMaxLen) ))
MaxFst = max(Data[,4]); MaxFst
ChrMidPos = (ChrLenCumsum[-1] + ChrLenCumsum[-length(ChrLenCumsum)])/2
plot(NULL, NULL,
     xlim = c(0, ChrLenCumsum[length(ChrLenCumsum)]/1000000),
    main = CHR, cex.main = 1.5,
    xlab = "(M bp)",
    ylim = c(0, 1),
    ylab=expression(F[ST]) , cex.lab = 1.2, 
    xaxt="n"
)
abline(v=0/1000000, col = "black")
abline(h=THR, col = "red", lwd = 1.5)
axis(1, at=ChrMidPos/1000000,
     labels = paste0(names(ChrMaxLen),"\n", as.integer(ChrMaxLen/1000000), "M"),
     tick = FALSE,
     las=1)

for( i in c(1:length(ChrMaxLen)) ) {
  chrID = names(ChrMaxLen)[i]
  X = Data[Data[,1]==chrID, 2] + ChrLenCumsum[i]
  Y = Data[Data[,1]==chrID, 4]
  points(X/1000000, Y,
       pch=16, cex = 0.2, 
       col = "grey2"
  )
  X1 = X[!is.na(Y)]
  Y1 = Y[!is.na(Y)]
  lines(smooth.spline(X1/1000000, Y1, spar=0.35), 
        col='blue', lwd=2)
  abline(v=ChrLenCumsum[i+1]/1000000, col = "black")
}

dev.off()


