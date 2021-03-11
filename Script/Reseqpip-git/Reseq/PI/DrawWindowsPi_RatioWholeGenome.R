#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2018-07-24
#
#
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-c", "--chr"), dest = "chr", default = "",
              help = "The chr name"),
  make_option(c("-s", "--size"), dest = "size", default = 50,
              help = "The window size (k bp)")
)
#
parser <- OptionParser(
    usage = "%prog -c chr1",
    option_list=option_list,
    description = 
      "Author: Guo, Weilong; guoweilong@126.com; 2018-08-04 \
      Description: Draw the ratios for two Pi distribution across whole genome"
)
#
arguments <- parse_args(parser)
#
CHR <- arguments$chr
size <- as.integer(arguments$size)
#
pdf( paste0( CHR, ".", size, "k.windowed.pi.pdf"), width = 24, height = 4)
TAB = read.table( paste0( CHR, ".", size, "k.windowed.pi2"), header = TRUE)
Data = TAB

ChrMaxLen = tapply(Data[,3], Data[,1], max)
ChrLenCumsum = cumsum( c(0, t(ChrMaxLen) ))
ChrMidPos = (ChrLenCumsum[-1] + ChrLenCumsum[-length(ChrLenCumsum)])/2
plot(NULL, NULL,
    xlim = c(0, ChrLenCumsum[length(ChrLenCumsum)]/1000000),
    main = CHR, cex.main = 1.5,
    xlab = "(M bp)",
    ylim = c(0, 60),
    ylab=expression(pi[1]/pi[2]) , cex.lab = 1.5,
    xaxt="n"
)
abline(v=0/1000000, col = "grey2")
axis(1, at=ChrMidPos/1000000,
     labels = paste0(names(ChrMaxLen),"\n", as.integer(ChrMaxLen/1000000), "M"),
     tick = FALSE,
     las=1)


for( i in c(1:length(ChrMaxLen)) ) {
  chrID = names(ChrMaxLen)[i]
  X = Data[Data[,1]==chrID, 2] + ChrLenCumsum[i]
  Y1 = Data[Data[,1]==chrID, 4]/Data[Data[,1]==chrID, 5]
  Y2 = Data[Data[,1]==chrID, 5]/Data[Data[,1]==chrID, 4]
  lines(X/1000000, Y1, col = "red2", lwd = 2)
  lines(X/1000000, Y2, col = "blue4", lwd = 2)
  abline(v=ChrLenCumsum[i+1]/1000000, col = "grey2")
}

legend("topright", legend = c( expression(pi[A]/pi[B]), expression(pi[B]/pi[A]) ), 
       col = c("red2", "blue4"), lty = 1, lwd = 2)
dev.off()


