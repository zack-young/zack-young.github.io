#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2018-01-07
#
#
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-c", "--chr"), dest = "chr", default = "",
              help = "The chr name"),
  make_option(c("-s", "--size1"), dest = "size1", default = 50,
              help = "The small window size"),
  make_option(c("-S", "--size2"), dest = "size2", default = 1000,
              help = "The larger window size")
#  make_option(c("-p", "--prefix"), dest = "prefix", default = "",
#              help = "The prefix for filename"),
#  make_option(c("-c", "--color"), dest = "color", default = "red2",
#              help = "The color for drawing the line")
)
#
parser <- OptionParser(
    usage = "%prog -c chr1A -s 50 -S 1000",
    option_list=option_list,
    description = 
      "Author: Guo, Weilong; guoweilong@126.com; 2018-08-09 \
      Description: Draw the distribution of Tajima's D for single chromosome"
)
#
arguments <- parse_args(parser)
#
#PRE <- arguments$prefix
CHR <- arguments$chr
size1 <- as.integer(arguments$size1)
size2 <- as.integer(arguments$size2)
#

#
pdf( paste0( CHR, ".Tajima.D.pdf"), width = 16, height = 4)
TAB = read.table( paste0( CHR, ".", size1, "k.Tajima.D"), header = TRUE)
Data = TAB[,c(2,4)]
Data = Data[!is.na(Data[,2]),]
plot(Data[,1]/1000000, Data[,2],
     pch=16, cex = 0.2,
     col = "grey2",
     main = CHR, cex.main = 1.5,
     ylim = c(-3, 5),
     xlab = "(M bp)",
     ylab = "Tajima's D", cex.lab = 1.5
)
TAB_L = read.table( paste0( CHR, ".", size2, "k.Tajima.D"), header = TRUE)
Data_L = TAB_L[,c(2,4)]
Data_L = Data_L[!is.na(Data_L[,2]),]
lines(Data_L[,1]/1000000, Data_L[,2], col = "red2", lwd = 2)
dev.off()



