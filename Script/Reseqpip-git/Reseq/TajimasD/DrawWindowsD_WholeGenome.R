#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2018-08-109
#
#
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-c", "--chr"), dest = "chr", default = "",
              help = "The chr name"),
  make_option(c("-w", "--window_size"), dest = "wSize", default = 200000,
              help = "The window size [default: %default]")
)
#
parser <- OptionParser(
    usage = "%prog -c chr1",
    option_list=option_list,
    description = 
      "Author: Guo, Weilong; guoweilong@126.com; 2018-08-09 \
      Description: Draw the distribution for Tajima's D"
)
#
arguments <- parse_args(parser)
#
CHR <- arguments$chr
wSize <- as.integer(arguments$wSize)
#CHR = "chrA.HA_vs_LA.200k"
#
pdf( paste0( CHR, ".Tajima.D.pdf"), width = 20, height = 4)
par( mai = c(0.5,1,0.5,0.2) )
TAB = read.table( paste0( CHR, ".Tajima.D"), header = TRUE)
Data = TAB[,c(1,2,3,4)]
ChrMaxLen = tapply(Data[,2], Data[,1], max)+wSize
ChrLenCumsum = cumsum( c(0, t(ChrMaxLen) ))
MaxTajimaD = max(Data[,4]); MaxTajimaD
ChrMidPos = (ChrLenCumsum[-1] + ChrLenCumsum[-length(ChrLenCumsum)])/2
plot(NULL, NULL,
     xlim = c(0, ChrLenCumsum[length(ChrLenCumsum)]/1000000),
    main = CHR, cex.main = 1.5,
    xlab = "(M bp)",
    ylim = c(-3, 5),
    ylab = "Tajima\'s D", cex.lab = 1.2, 
    xaxt="n"
)
abline(v=0/1000000, col = "red")
axis(1, at=ChrMidPos/1000000,
     labels = paste0(names(ChrMaxLen),"\n", as.integer(ChrMaxLen/1000000), "M"),
     tick = FALSE,
     las=1)
#
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
  abline(v=ChrLenCumsum[i+1]/1000000, col = "red")
}
#
dev.off()
#

