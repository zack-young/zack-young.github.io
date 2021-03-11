#!/usr/bin/env Rscript

#
#SAMPLE1=arguments$sample1
#SAMPLE1_name=arguments$sample1_name
#SAMPLE2=arguments$sample2
#SAMPLE2_name=arguments$sample2_name
options(scipen=200)
files = read.table("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plotfile_6811_parent", as.is = T, header = F, comment.char = "")

#pdf("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/pdf/parent_define.pdf", height = 9, width = 12)
plotChrom <- function(xleft, ybottom, height, len, centro, binsize){
  
  #len <- len/binsize
  #centro <- centro/binsize
  
  # r vertical, r horizontal
  rv <- height/2
  rh <- len/70
  rs <- seq(0,pi,len=100)
  
  # left semi-circle
  lx <- c(xleft, xleft + rh - rh*sin(rs), xleft, xleft)
  ly <- c(ybottom, ybottom + rv - rv*cos(rs), ybottom + height, ybottom)
  polygon(lx, ly, border = "white", col = "white")
  
  # right semi-circle
  rx <- c(xleft + len,xleft + len, xleft + len - rh + rh*sin(rs), xleft + len)
  ry <- c(ybottom, ybottom + height, ybottom + rv +rv*cos(rs), ybottom)
  polygon(rx, ry, border = "white", col = "white")
  
  # top tri
  tx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  ty <- c(ybottom + height, ybottom + rv, ybottom + height)
  polygon(tx, ty, border = "white", col = "white")
  
  # bottom tri
  bx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  by <- c(ybottom, ybottom + rv, ybottom)
  polygon(bx, by, border = "white", col = "white")
  
  cbx <- c(xleft + rh - rh*sin(rs), centro - 0.5*rh, centro, centro + 0.5*rh, xleft + len - rh + rh*sin(rs), centro + 0.5*rh, centro, centro - 0.5*rh)
  cby <- c(ybottom + rv - rv*cos(rs), ybottom + height, ybottom + rv, ybottom + height, ybottom + rv + rv*cos(rs), ybottom, ybottom + rv, ybottom)
  polygon(cbx, cby, border = "grey", lwd = 3)
}

plot(x=0, type="n", bty="n", yaxt="n",xaxt="n",
     xlab="", ylab="", 
     xlim=c(-200, 1100), ylim=c(-2, nrow(files)+1),
     xaxs="i", yaxs="i", 
     #main=paste("S",j,sep=""),
     main="JiMai44 parent define",
     cex.main = 2.5)
text(x=c(0,200,400,600,800), y=rep(0,5),
     c("0Mbp","200Mbp","400Mbp","600Mbp","800Mbp"), cex = 1,pos = 1)
#color_pad <- c(colorpanel(low="yellow", high="red", n=3),"green","black","white")
fa_snp <- "#3E3A92"
fa_CNV <- "#39BAEC"
ma_snp <- "#F6594C"
ma_CNV <- "#F8C0BB"
CNV_un <- "#22221F"
SNP_un <- "#CDCDC9"
wait <- "#39BAEC"

color_pad <- c(fa_snp,fa_CNV,ma_snp,ma_CNV,CNV_un,SNP_un ,wait)
#color_pad <- c('white','white','white',red,yellow,purple)
num = 22
centromere <- c(215,239,173,337,346,268,300,346,241,300,317,184,254,206,189,286,327,211,357,287,340)
wholen <- c(594102056,
            689851870,
            495453186,
            780798557,
            801256715,
            651852609,
            750843639,
            830829764,
            615552423,
            744588157,
            673617499,
            509857067,
            709773743,
            713149757,
            566080677,
            618079260,
            720988478,
            473592718,
            736706236,
            750620385,
            638686055
)
rect(xleft = 900, xright = 930, ybottom = 20, ytop = 22,
     col = color_pad[1], border = color_pad[1]);text(x=c(940), y=c(21),c("JiNan17_SNP"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 22, ytop = 24,
     col = color_pad[3], border = color_pad[3]);text(x=c(940), y=c(23),c("954072_SNP"),cex = 1,pos = 4)

rect(xleft = 900, xright = 930, ybottom = 26, ytop = 28,
     col = color_pad[2], border = color_pad[2]);text(x=c(940), y=c(27),c("JiNan17_CNV"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 28, ytop = 30,
     col = color_pad[4], border = color_pad[4]);text(x=c(940), y=c(29),c("954072_CNV"),cex = 1,pos = 4)
rect(xleft = 900, xright = 930, ybottom = 32, ytop = 34,
     col = color_pad[6], border = color_pad[6]);text(x=c(940), y=c(33),c("Undefined"),cex = 1,pos = 4)


for(j in 1:nrow(files)){
  if(!is.na(files[j,1])){
    data1 <- read.table(files[j,1], as.is = T, header = T, comment.char = "")
    #chro <- factor(data1[,5])
    #              if(!is.na(files[j,2])){
    #                chro <- factor(data1[,7])
    #              } else{
    #                chro <- factor(data1[,5])
    #              }
    chro <- factor(data1[,4])
    # chro1 <- factor(chro, levels = c('CNV', 'DP', 'hete', 'homo'), 
    #                 labels = c('1', '2', '3', '4'))
    chro1 <- factor(chro, 
                    #levels = c("level1","level2","level3","unmatchCNV","DP","ownCNV"), 
                    #
                    levels = c("JiNan17_snp","JiNan17_CNV","954072_snp","954072_CNV","undefined_CNV","undefined_snp","wait_define"),
                    labels = c('1','2','3','4',"6","6","6")
                    #levels = c("jm22","both","lx99"),
                    #labels = c('4',"5","6")
    )
    #chro1 <- factor(chro, levels = c('level1', 'level2', 'level3', 'level4',"level5","level6","level7","level8","level9","level10"), 
    #labels = c('10', '9', '8', '7','6', '5', '4', '3','2', '1')) 
    
    num <- num - 1 
    #assign(paste("sample",num,sep=""),nrow(data1))
    for(i in 1:nrow(data1)){
      rect(xleft = as.numeric(data1[i,2])/1000000,
           #         #ybottom = j+0.1,
           ybottom = j+0.1,
           ytop = j+0.9,
           xright = as.numeric(data1[i,3])/1000000,
           #         #ytop = j+0.9,
           col = color_pad[as.numeric(as.vector(chro1[i]))],
           border = color_pad[as.numeric(as.vector(chro1[i]))] 
      )
    }
    len <- wholen[num]/1000000
    centro <- centromere[num]
    #plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
    #            if(!is.na(files[j,2])){
    #              plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
    #            }
    
    text(x=-100, y=j+0.5, files[j,2], cex = 1)
    # text(x=c(950,950,950,935), y=c(20.5,21.5,22.5,23.5),
    #      c("lx987","nd3097","undefined","CNV DP Hete"), cex = 0.7,pos = 4)
    #            text(x=c(940,940,940,940), y=c(20.5,21.5,22.5,23.5),
    #                 c("level1","level2","level3","level4"), cex = 0.7,pos = 4)
    #text(x=c(940,940,940,940,940,940,940,940,940,940), y=c(20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5),
    #c("level10","level9","level8","level7","level6","level5","level4","level3","level2","level1"),cex = 0.7,pos = 4)
    #rect(xleft = 900, xright = 930, ybottom = 26, ytop = 28,
    #     col = color_pad[4], border = color_pad[4])
    #rect(xleft = 900, xright = 930, ybottom = 28, ytop = 30,
    #     col = color_pad[5], border = color_pad[5])
    
    # rect(xleft = 900, xright = 930, ybottom = 26, ytop = 27,
    #      col = color_pad[7], border = color_pad[7])
    # rect(xleft = 900, xright = 930, ybottom = 27, ytop = 28,
    #      col = color_pad[8], border = color_pad[8])
    # rect(xleft = 900, xright = 930, ybottom = 28, ytop = 29,
    #      col = color_pad[9], border = color_pad[9])
    # rect(xleft = 900, xright = 930, ybottom = 29, ytop = 30,
    #      col = color_pad[10], border = color_pad[10])
    
    # rect(xleft = 900, xright = 930, ybottom = 23, ytop = 24,
    #      col = color_pad[4], border = color_pad[4])
    
  }
}
num = 22
for(j in 1:nrow(files)){
  if(is.na(files[j,1])){
    j <- j-1
    num <- num - 1 
    centro <- centromere[num]
    len <- wholen[num]
    #centro = 200
    #j=3
    len <- wholen[num]/1000000
    plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
    #tx <- c(centro+5,centro,centro-5)
    #ty <- c(j-3.7,j-3,j-3.7)
    #polygon( tx, ty, border = "black", col = "black")
    
  } 
}
speci_pos <- function(x,y,gene){
  tx <- c(x+2,x,x-2)
  ty <- c(y-0.7,y,y-0.7)
  polygon( tx, ty, border = red, col = red)
}

pos_oper <- function(x){(x-1)*2+1}
ya <- pos_oper(21)
speci_pos(508,ya) #1A
text(x=c(508+2.2), y=c(ya-0.5),c("Glu-A1"), cex = 0.5,pos = 4,col=red)
speci_pos(4,ya)
text(x=c(4+2.2), y=c(ya-0.5),c("Glu-A3"), cex = 0.5,pos = 4,col=red)
ya <- pos_oper(20)
speci_pos(555,ya) #1B
text(x=c(555+2.2), y=c(ya-0.5),c('Bx7OE;Bx13'), cex = 0.5,pos = 4,col=red)
speci_pos(5,ya)
text(x=c(5+2.2), y=c(ya-0.5),c('Glu-B3'), cex = 0.5,pos = 4,col=red)
ya <- pos_oper(19)
speci_pos(412,ya) #1D
text(x=c(412+2.2), y=c(ya-0.5),c('Glu-D1'), cex = 0.5,pos = 4,col=red)
ya <- pos_oper(18)
speci_pos(712,ya) #2A
text(x=c(712+2.2), y=c(ya-0.5),c('ppo-A1'), cex = 0.5,pos = 4,col=red)
speci_pos(321,ya)
text(x=c(321+2.2), y=c(ya-0.5),c('zds-A1'), cex = 0.5,pos = 4,col=red)
ya <- pos_oper(15)
speci_pos(448,ya) #3A
text(x=c(448+2.2), y=c(ya-0.5),c('lyce-3A'), cex = 0.5,pos = 4,col=red)
ya <- pos_oper(14)
speci_pos(148,ya) #3B
text(x=c(148+2.2), y=c(ya-0.5),c('Talyc-B1'), cex = 0.5,pos = 4,col=red)
ya <- pos_oper(12)
speci_pos(688,ya) #4A
text(x=c(688+2.2), y=c(ya-0.5),c('Wx-B1'), cex = 0.5,pos = 4,col=red)
ya <- pos_oper(11)
speci_pos(612,ya) #4B
text(x=c(612+2.2), y=c(ya-0.5),c('Tapds-B1'), cex = 0.5,pos = 4,col=red)
ya <- pos_oper(7)
speci_pos(3,ya) #5D
text(x=c(3+2.2), y=c(ya-0.5),c('Pina-D1;Pinb-D1'), cex = 0.5,pos = 4,col=red)
ya <- pos_oper(3)
speci_pos(729,ya) #7A
text(x=c(729+2.2), y=c(ya-0.5),c('Psy-A1'), cex = 0.5,pos = 4,col=red)
ya <- pos_oper(2)
speci_pos(739,ya) #7B
text(x=c(739+2.2), y=c(ya-0.5),c('Psy-B1'), cex = 0.5,pos = 4,col=red)


#dev.off()

