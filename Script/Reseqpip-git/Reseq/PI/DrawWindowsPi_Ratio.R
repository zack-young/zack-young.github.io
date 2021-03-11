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
      Description: Draw the ratios for two Pi distribution"
)
#
arguments <- parse_args(parser)
#
CHR <- arguments$chr
size <- as.integer(arguments$size)
#
pdf( paste0( CHR, ".windowed.pi.pdf"), width = 16, height = 4)
TAB = read.table( paste0( CHR, ".", size, "k.windowed.pi2"), header = FALSE)
Data = TAB[,c(2,4,5)]
X=Data[,1]/1000000
Y1=Data[,2]/Data[,3]
Y2=Data[,3]/Data[,2]
plot(X, Y1,
     type="l", lwd = 2,
     cex = 0.2, 
     col = "red2",
     ylim = c(0, 60),
     main = CHR, cex.main = 1.5,
     xlab = "(M bp)",
     ylab=expression(pi[1]/pi[2]), 
     cex.lab = 1.5
)
lines(X, Y2, col = "blue4", lwd = 2)
legend("topright", legend = c( expression(pi[A]/pi[B]), expression(pi[B]/pi[A]) ), 
       col = c("red2", "blue4"), lty = 1, lwd = 2)
dev.off()

# ==
# write the ratio to file

#write.table( cbind( CHR=TAB[,1], START=TAB[,2], END=TAB[,3], 
#                   PIa_to_PIb = format(Y1, digit=4),  
#                   PIb_to_PIa = format(Y2, digit=4) ),
#             file =  paste0( CHR, ".", size, "k.windowed.piratio.xls"),
#             quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
# )


