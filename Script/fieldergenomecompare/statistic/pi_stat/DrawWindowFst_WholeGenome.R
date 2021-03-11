#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help = "The chr name"),
  make_option(c("-c", "--color"), dest = "colorS", default = 1,
              help = "The chr name"),
  make_option(c("-w", "--window_size"), dest = "wSize", default = 200000,
              help = "The window size [default: %default]"),
  make_option(c("-t", "--threshold"), dest = "thr", default = 0.5,
              help = "The threshold for drawing the horizontal line  [default: %default]"),
  make_option(c("-g", "--genefile"), dest = "genefile", default = "/data2/rawdata2/GenomeScaning/Fst/scripts/DEDO_genelist.txt",
              help = "The threshold for drawing the horizontal line  [default: %default]"),
  make_option(c("-n", "--datacol"), dest = "datacol", default = "5",
              help = "The threshold for drawing the horizontal line  [default: %default]"),
  make_option(c("-a", "--astat"), dest = "astat", default = "Fst",
              help = "The threshold for drawing the horizontal line  [default: %default]"),
  make_option(c("-Y", "--yuplimP"), dest = "yuplimP", default = 1,
              help = "The threshold for drawing the horizontal line  [default: %default]"),
  make_option(c("-y", "--ylolimP"), dest = "ylolimP", default = 0,
              help = "The threshold for drawing the horizontal line  [default: %default]"),
  make_option(c("-T", "--titleP"), dest = "titleP", default = "Fst-scanning",
              help = "The threshold for drawing the horizontal line  [default: %default]"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "outfile",
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
infile <- arguments$infile
colorS <- as.integer(arguments$colorS)
THR <- as.double(strsplit(as.character(arguments$thr), split = ",")[[1]])
datacol <- as.integer(arguments$datacol)
genefile <- arguments$genefil
astat <- arguments$astat
yuplimP <- as.double(arguments$yuplimP)
ylolimP <- as.double(arguments$ylolimP)
titleP <- arguments$titleP
outfile <- arguments$outfile

if (!is.na(genefile)) genelist <- read.table(genefile)
#CHR = "chrA.HA_vs_LA.200k"
#
# define colors
# red: #ED1C24, #EA841F
# yellow: #FFDE17, #00A14B
# pur: #7F3F98, #21409A
# pur: "#8C4E9E", "#1B8EB0"
pdf( outfile, width = 10, height = 4)
par( mai = c(0.8,0.5,0.5,0.2) ,mgp=c(2,1,0))
TAB = read.table(infile, header = TRUE)
Data = TAB[,c(1,2,3,datacol)]
ChrMaxLen = tapply(Data[,3], Data[,1], max)
ChrLenCumsum = cumsum( c(0, t(ChrMaxLen) ))
MaxFst = max(Data[,4])
ChrMidPos = (ChrLenCumsum[-1] + ChrLenCumsum[-length(ChrLenCumsum)])/2




  mainP <- titleP

plot(NULL, NULL,
    xlim = c(0, ChrLenCumsum[length(ChrLenCumsum)]/1000000),
    main = mainP, cex.main = 2,
    xlab = "Chromosome",
    ylim = c(ylolimP, yuplimP),
    ylab=NA , cex.lab = 2, 
    xaxt="n",
    yaxt="n",
    frame.plot = FALSE,
)
#
axis(2,pos = 0, tck=-0.005, lwd = 2, las=2, cex.axis=1.3)
axis(1, at=ChrMidPos/1000000, pos = ylolimP, labels = names(ChrMaxLen), cex.axis=1.3)
segments(0, ylolimP, ChrLenCumsum[8]/1000000,ylolimP,lwd = 2)
for (i in 1:length(THR)){
  segments(0, THR[i], ChrLenCumsum[8]/1000000, THR[i], lwd = 2, lty=2, col = "darkgray")
}
segments(0, ylolimP, 0, yuplimP, lwd = 2)
title(ylab = astat, mgp = c(0.3, 1, 0), cex.lab = 2)

for( i in c(1:length(ChrMaxLen)) ) {
  if (colorS == 1){
    if (i%%2 == 1) color <- "#ED1C24" else color <- "#EA841F"
  } else if (colorS == 2) {
    if (i%%2 == 1) color <- "#46C7F4" else color <- "#00A14B"
  } else {
    if (i%%2 == 1) color <- "#8C4E9E" else color <- "#4043FF"
  }
  chrID = names(ChrMaxLen)[i]
  X = Data[Data[,1]==chrID, 2] + ChrLenCumsum[i]
  Y = Data[Data[,1]==chrID, 4]
  points(X/1000000, Y,
       pch=16, cex = 0.5, 
       col = color
  )
  X1 = X[!is.na(Y)]
  Y1 = Y[!is.na(Y)]
  lines(smooth.spline(X1/1000000, Y1, spar=0.35), 
        col="darkgray", lwd=2)
  #
  if (!is.na(genefile)){
    genes <- genelist[genelist$V1 == chrID,]
    #print("there is gene") 
    if (nrow(genes) != 0){
      for (j in 1:nrow(genes)){
        gX <- (genes[j,2]+ChrLenCumsum[i])/1000000
        arrows(gX, yuplimP, gX, yuplimP-(yuplimP-ylolimP)*0.1, col="black", length = 0.1, lwd=1.5)
        if (j==1) POS=2 else POS=4 # special case only
        text(gX, yuplimP*0.95, genes[j,3], pos = POS, offset=2, cex=2)
      }
    }
  }
}

dev.off()
