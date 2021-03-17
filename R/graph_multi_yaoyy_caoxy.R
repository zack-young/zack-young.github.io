library(gplots)
files = read.table("/data/user/yangzz/mapping/fieldergenomecompare/202009_11_yaoyy/plotfile_CP07_CP09_CP11", as.is = T, header = F, comment.char = "")
qrc_lis = read.table("/data/user/yangzz/mapping/fieldergenomecompare/202009_11_yaoyy/statistic/gene_list.txt", as.is = T, header = F, comment.char = "")

options(scipen=200)
#pdf(paste("C:/Users/lenovo/Desktop/2019_3_6/T",j,"_both_ratio_new.pdf",sep=""), height = 10, width = 15.69)
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
  polygon(lx, ly, border = NA, col = "white")
  
  # right semi-circle
  rx <- c(xleft + len,xleft + len, xleft + len - rh + rh*sin(rs), xleft + len)
  ry <- c(ybottom, ybottom + height, ybottom + rv +rv*cos(rs), ybottom)
  polygon(rx, ry, border = NA, col = "white")
  
  # top tri
  tx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  ty <- c(ybottom + height, ybottom + rv, ybottom + height)
  polygon(tx, ty, border = NA, col = "white")
  
  # bottom tri
  bx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  by <- c(ybottom, ybottom + rv, ybottom)
  polygon(bx, by, border = NA, col = "white")
  
  cbx <- c(xleft + rh - rh*sin(rs), centro - 0.5*rh, centro, centro + 0.5*rh, xleft + len - rh + rh*sin(rs), centro + 0.5*rh, centro, centro - 0.5*rh)
  cby <- c(ybottom + rv - rv*cos(rs), ybottom + height, ybottom + rv, ybottom + height, ybottom + rv + rv*cos(rs), ybottom, ybottom + rv, ybottom)
  polygon(cbx, cby, border = "grey", lwd = 1)
}

low1 <- "#F9B731"
low2 <- "#3689EA"
delboth <-  "#007046"
dupboth <- "#AA0000"
undefined <- "#7B7B7B"
#color_pad <- c(low,mid,high,del1,del2,delboth,dup1,dup2,dupboth)
color_pad <- c(low1,low2,low1,low2,delboth,dupboth,undefined)
color_pad1 <- c(rgb(1, 0, 0, 0.5),rgb(0, 1, 0, 0.5),rgb(0, 0, 1, 0.5))
plot(x=0, type="n", bty="n", yaxt="n",xaxt="n",
     xlab="", ylab="", 
     xlim=c(-200, 1100), ylim=c(-5, nrow(files)+1),
     xaxs="i", yaxs="i", 
     #main=paste("S",j,sep=""),
     #main="ji200040919 jimai229 954072 Unmatch Level",
     cex.main = 3)
text(x=c(0,200,400,600,800), y=rep(-1,5),
     c("0Mbp","200","400","600","800"), cex = 1.5,pos = 1)
segments(x0=c(0,200,400,600,800),
         x1=c(0,200,400,600,800),lwd=2,
         y0 = 0,
         y1 = -1,
         xpd = T
)
segments(x0=0,
         x1 = 800,
         y0 = 0,
         y1 = 0,lwd=2,
         xpd = T
)
#text(x=c(940,940,940), y=c(21,24,27),
#     c("Similar","LOW_Diff","HIGH_Diff"), cex = 0.8,pos = 4)
#text(x=rep(940,6), y=c(30,33,36,21,24,27),
#     c("Both Deletion","Both Duplication","Undefined","JiMai22","YanNong15","SN224"), cex = 0.8,pos = 4)
text(x=rep(860,8), y=seq(18,39,3),#x=seq(-190,-50,20), y=rep(30.5,8),
     c("Undefiend","ji200040919","954072","Both Deletion","Both Duplication","HMW","LMW","Gliadin"),cex = 0.8,pos=4)#,srt=90)

rect(xleft = rep(850,5), xright = rep(860,5),
     ybottom = seq(17.5,30.5,3),
     ytop = seq(18.5,31.5,3),
     col = c(undefined,low1,low2,delboth,dupboth), border = NA)

x=rep(855,3)
y=seq(32.5,38.5,3)
polygon( c(x[1]+5,x[1],x[1]-5), c(y[1]+1,y[1],y[1]+1), border = NA, col = color_pad1[1])
polygon( c(x[2]+5,x[2],x[2]-5), c(y[2]+1,y[2],y[2]+1), border = NA, col = color_pad1[2])
polygon( c(x[3]+5,x[3],x[3]-5), c(y[3]+1,y[3],y[3]+1), border = NA, col = color_pad1[3])

#color_pad <- c(blue,red)

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
num = 22
for(j in 1:nrow(files)){
  if(!is.na(files[j,1])){
    text(x=-50, y=j+1, substr(files[j,2],4,5), cex = 1.5)#,srt = 90)
    if(file.exists(files[j,1]) == 0) {
      text(x=-100, y=j+0.5, files[j,2], cex = 1.5)
      #plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
    } else{
      data1 <- read.table(files[j,1], as.is = T, header = F, comment.char = "")
      chro <- factor(data1[,4])
      #SAMPLE1 <- strsplit(strsplit(files[j,1],split = '/')[[1]][8],split = '_')[[1]][1]
      #SAMPLE2 <- strsplit(strsplit(files[j,1],split = '/')[[1]][8],split = '_')[[1]][2]
      SAMPLE1 <- strsplit(strsplit(files[j,1],split = '/')[[1]][9],split = '_')[[1]][2]
      #      chro1 <- factor(chro, 
      #                      levels = c("low1","low2","low3",'deletion_both_CNV','duplication_both_CNV','undefined'),
      #                      labels = c('1','2','3','4','5','6')
      #      )
      chro1 <- factor(chro, 
                      levels = c("low1","low2","low3","low4",'deletion_both_CNV','duplication_both_CNV',"undefined"),
                      #levels = c("low","mid","high",
                      #           paste(SAMPLE1,"deletion_CNV",sep=""),paste(SAMPLE2,"deletion_CNV",sep=""),"deletion_both_CNV",
                      #           paste(SAMPLE1,"duplication_CNV",sep=""),paste(SAMPLE2,"duplication_CNV",sep=""),"duplication_both_CNV"),
                      #labels = c('3','1','1','4',"5","6")
                      labels = c('1','2','3',"4","5",'6',"7")
                      #levels = c("jm22","both","lx99"),
                      #labels = c('4',"5","6")
      )
      if(!is.na(files[j,2])){
        #points(x=-10,y=j+0.45,pch=2,cex=0.7)
        num <- num - 1
        centro <- centromere[num]
        len <- wholen[num]
        #centro = 200
        #j=3
        len <- wholen[num]/1000000
        for(i in 1:nrow(data1)){
          rect(xleft = as.numeric(data1[i,2])/1000000,
               #         #ybottom = j+0.1,
               ybottom = (j-0.3),
               ytop = (j+0.9),
               xright = as.numeric(data1[i,3])/1000000,
               #         #ytop = j+0.9,
               col = color_pad[as.numeric(as.vector(chro1[i]))],
               border = NA#color_pad[as.numeric(as.vector(chro1[i]))] 
          )
        }
        plotChrom(0, j-0.3, 1.25, len, centro, 1000000)
        qrc_dt <- qrc_lis[which(qrc_lis[,1]==files[j,2]),]
        if (nrow(qrc_dt)!=0){
          chro_gene <- factor(qrc_dt[,4])
          chro_gene2 <- factor(chro_gene, 
                          levels = c("HMW","LMW","GLIA"),
                          labels = c('1','2','3')

          )
          for(a in 1:nrow(qrc_dt)){
            x <- as.numeric(qrc_dt[a,2])/1000000
            tx <- c(x+1,x,x-1)
            ty <- c(j+2.9,j+2.3,j+2.9)
            polygon( tx, ty, border = NA, col = color_pad1[as.numeric(as.vector(chro_gene2[a]))])
            #text(x=(as.numeric(qrc_dt[a,2])+as.numeric(qrc_dt[a,3]))/2,y=j+2.5,a)
          }}
      }else{
        #points(x=-10,y=j+0.75,pch=0,cex=0.7)
        for(i in 1:nrow(data1)){
          if (files[j,3]==1){
            ybow=j+0.4
          }else{
            ybow=j-0.15
          }
          rect(xleft = as.numeric(data1[i,2])/1000000,
               #         #ybottom = j+0.1,
               ybottom = ybow,
               ytop = ybow+0.3,
               xright = as.numeric(data1[i,3])/1000000,
               #         #ytop = j+0.9,
               col = color_pad[as.numeric(as.vector(chro1[i]))],
               border = NA#color_pad[as.numeric(as.vector(chro1[i]))] 
          )
          
        }
      }
      #text(x=900, y=j+0.5, lis[SAMPLE1], cex = 0.75,adj=0)
      #text(x=-100, y=j+2, files[j,2], cex = 1)
      
    }
  }
}
# num = 22
# for(j in 1:nrow(files)){
#   if(is.na(files[j,1])){
#     j <- j-1
#     num <- num - 1 
#     centro <- centromere[num]
#     len <- wholen[num]
#     #centro = 200
#     #j=3
#     len <- wholen[num]/1000000
#     #plotChrom(0, j+0.1, 0.8, len, centro, 1000000)
#     tx <- c(centro+5,centro,centro-5)
#     ty <- c(j-1.7,j-1,j-1.7)
#     polygon( tx, ty, border = "black", col = "black")
#     
#   } 
# }
num = 22
for(j in 1:nrow(files)){
  if(is.na(files[j,1])){
    j <- j-4
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

pos_oper <- function(x){(x-1)*3+1}
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



