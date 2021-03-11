files = read.table("/data/user/yangzz/mapping/fieldergenomecompare/10+_genome/chr6A_hap_data.txt", as.is = T, header = F, comment.char = "")
meta <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/1.diff_dev/metadata_cultivar_10+.txt", as.is = T, header = F, comment.char = "")
bedfile <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/1.diff_dev/chr6A.1M.bed", as.is = T, header = F, comment.char = "")
hap_dt <- as.data.frame(setNames(replicate(8+2,numeric(0), simplify = F), c("CHR","POS",meta$V1)))
hap_dt[,"CHR"] <- files$V1
hap_dt[,"POS"] <- files$V2
for (item in meta$V1){
  for (i in 1:nrow(files)){
    hap_dt[i,item] <- files[i,which(files[i,]==item)+1]
  }
}
hap_dt <- hap_dt[c(1,2,3,5,8,9,4,6,7,10)]
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

del1 <- "#60D6A9"
del2 <- "#9FEE00"
delboth <-  "#007046"
dup1 <- "#FF8B73"
dup2 <- "#FFF073"
dupboth <- "#AA0000"
high <- "#111148"
mid <-"#334899"
low <- "#22BBEE"
color_pad <- c("#008081","#FF7F2B","#DE87CC",del1,del2,"white")
centromere <- c(215,239,173,337,346,268,300,346,241,300,317,184,254,206,189,286,327,211,357,287,340)
wholen <- c(594102056,689851870,495453186,
            780798557,801256715,651852609,
            750843639,830829764,615552423,
            744588157,673617499,509857067,
            709773743,713149757,566080677,
            618079260,720988478,473592718,
            736706236,750620385,638686055)
plot(x=0, type="n", bty="n", yaxt="n",xaxt="n",
     xlab="", ylab="", 
     xlim=c(-200, 1100), ylim=c(-5, ncol(hap_dt)+5),
     xaxs="i", yaxs="i", 
     #main=paste("S",j,sep=""),
     #main=paste(SAMPLE1_name,"VS",SAMPLE2_name,sep=" "),
     cex.main = 2.5)
axis(1,at=c(0,100,200,300,400,500,600,700,800),labels = c("0Mbp","","200Mbp","","400Mbp","","600Mbp","","800Mbp"))

for(j in 3:ncol(hap_dt)){
    chro <- factor(hap_dt[,j])
    chro1 <- factor(chro, 
                    levels = c(paste(rep("type",5),seq(1,5),sep = ""),"other"),
                    labels = c('1','2','3',"4","5",'6'))
    for(k in 1:nrow(hap_dt)){
      rect(xleft = as.numeric(hap_dt[k,2])/1000000,
           #         #ybottom = j+0.1,
           ybottom = j+0.1,
           ytop = j+0.9,
           xright = as.numeric(hap_dt[k,2])/1000000+1,
           #         #ytop = j+0.9,
           col = color_pad[as.numeric(as.vector(chro1[k]))],
           border = color_pad[as.numeric(as.vector(chro1[k]))] 
      )
    }
    #len <- wholen[num]/1000000
    #centro <- centromere[num]
    #plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
    #            if(!is.na(files[j,2])){
    #              plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
    #            }
    
    text(x=-100, y=j+0.5, colnames(hap_dt)[j], cex = 1)
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

num = 22
tx <- c(centro+5,centro,centro-5)
ty <- c(j-0.7,j,j-0.7)
polygon( tx, ty, border = "black", col = "black")

for(j in 3:ncol(hap_dt)){
    #j <- j-1
    #num <- num - 1
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
