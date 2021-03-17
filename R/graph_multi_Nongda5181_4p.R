library(gplots)
files = read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/NongDa5181/plotfile_nd5181", as.is = T, header = F, comment.char = "")
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
deepblue <- "#282658"
midblue <- "#32519D"
lightblue <- "#39BAEC"
lightred <- "#F6E5D3"
red<-'#F784B6'
purple <- "#C471F2"
lightpurple <- "#9C7AF0"
lightyellow <- "#FFFEDF"
yellow <- "#F9B731"
del1 <- "#60D6A9"
del2 <- "#9FEE00"
delboth <-  "#007046"
dup1 <- "#FF8B73"
dup2 <- "#FFF073"
dupboth <- "#AA0000"
high <- "#111148"
mid <-"#334899"
low <- "#22BBEE"
S67 <- "#A5D6F2"
S6554 <- "#F49DA9"
YM03_undefined <- "#BEB8BB"
S67_S6554_shared <- "#5AE21A"
NZ05 <- "#266534"
S140_undefined <- "#BEB8BB"
YM03_NZ05_shared <- "#5AE21A"
YM03 <- "#DD80E1"

qrc_lis <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/QRC_list.txt",as.is = T, header = F, comment.char = "",sep = '\t')
#color_pad <- c(low,mid,high,del1,del2,delboth,dup1,dup2,dupboth)
color_pad <- c(S67,S6554,YM03_undefined,S67_S6554_shared,NZ05,YM03,S140_undefined,YM03_NZ05_shared)

plot(x=0, type="n", bty="n", yaxt="n",xaxt="n",
     xlab="", ylab="", 
     xlim=c(-200, 1100), ylim=c(-5, nrow(files)+1),
     xaxs="i", yaxs="i", 
     #main=paste("S",j,sep=""),
     #main="ji200040919 jimai229 954072 Unmatch Level",
     cex.main = 3)
text(x=c(0,200,400,600,800), y=rep(-1.5,5),
     c("0Mbp","200","400","600","800"), cex = 1.5,pos = 1)
segments(x0=c(0,200,400,600,800),
         x1=c(0,200,400,600,800),lwd=2,
         y0 = -1,
         y1 = -2,
         xpd = T
)
segments(x0=0,
         x1 = 800,
         y0 = -1,
         y1 = -1,lwd=2,
         xpd = T
)
#text(x=c(940,940,940), y=c(21,24,27),
#     c("Similar","LOW_Diff","HIGH_Diff"), cex = 0.8,pos = 4)
#text(x=rep(940,6), y=c(30,33,36,21,24,27),
#     c("Both Deletion","Both Duplication","Undefined","JiMai22","YanNong15","SN224"), cex = 0.8,pos = 4)
text(x=rep(920,4), y=seq(18,27,3),#x=seq(-190,-50,20), y=rep(30.5,8),
     c("LunXuan987","NongDa3097","Unknown","Shared"),cex = 0.8,pos=4)#,srt=90)

rect(xleft = rep(920,4), xright = rep(930,4),
     ybottom = seq(17.5,27.5,3),
     ytop = seq(18.5,28.5,3),
     col = color_pad[5:8], border = color_pad[5:8])
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
lis <- c("NongDa3097","NongDa5181")
names(lis) <-c("nd3097","nd5181")
for(j in 1:nrow(files)){
#for(j in 1:12){
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
                      levels = c("S67","S6554","YM03_unknown","S67_S6554_both","NZ05","YM03","S140_unknown","YM03_NZ05_both"),
                      #levels = c("low","mid","high",
                      #           paste(SAMPLE1,"deletion_CNV",sep=""),paste(SAMPLE2,"deletion_CNV",sep=""),"deletion_both_CNV",
                      #           paste(SAMPLE1,"duplication_CNV",sep=""),paste(SAMPLE2,"duplication_CNV",sep=""),"duplication_both_CNV"),
                      #labels = c('3','1','1','4',"5","6")
                      labels = c('1','2','3',"4","5",'6',"7","8")
                      #levels = c("jm22","both","lx99"),
                      #labels = c('4',"5","6")
      )
      if(!is.na(files[j,2])){
        #points(x=-10,y=j+0.45,pch=2,cex=0.7)
        num <- num - 1
        centro <- centromere[num]
        #len <- wholen[num]
        #centro = 200
        #j=3
        len <- wholen[num]/1000000
        for(i in 1:nrow(data1)){
          rect(xleft = as.numeric(data1[i,2])/1000000,
               #         #ybottom = j+0.1,
               ybottom = (j-1),
               ytop = (j+0.2),
               xright = as.numeric(data1[i,3])/1000000,
               #         #ytop = j+0.9,
               col = color_pad[as.numeric(as.vector(chro1[i]))],
               border = NA#color_pad[as.numeric(as.vector(chro1[i]))] 
          )
        }
        plotChrom(0, j-1, 1.2, len, centro, 1000000)
        qrc_dt <- qrc_lis[which(qrc_lis[,1]==files[j,2]),]
        if (nrow(qrc_dt)!=0){
        for(a in 1:nrow(qrc_dt)){
          rect(xleft = as.numeric(qrc_dt[a,2]),
            ybottom = (j+1.9),
            ytop = (j+2.5),
            xright = as.numeric(qrc_dt[a,3]),
             #         #ytop = j+0.9,
            col = "red",
            border = NA)
          #text(x=(as.numeric(qrc_dt[a,2])+as.numeric(qrc_dt[a,3]))/2,y=j+2.5,a)
        }
        
          }
      }else{
        #points(x=-10,y=j+0.75,pch=0,cex=0.7)
        for(i in 1:nrow(data1)){
          if (files[j,3]==1){
            ybow=j-0.4
          }else{
            ybow=j-0.75
            
          }
          rect(xleft = as.numeric(data1[i,2])/1000000,
               #         #ybottom = j+0.1,
               ybottom = ybow,
               ytop = ybow+0.6,
               xright = as.numeric(data1[i,3])/1000000,
               #         #ytop = j+0.9,
               col = color_pad[as.numeric(as.vector(chro1[i]))],
               border = NA#color_pad[as.numeric(as.vector(chro1[i]))] 
          )
          
        }
        if (files[j,3]==1){
          text(x=len, y=j, "SGR with Jingdong6", cex = 0.5 ,pos=4)
          text(x=len, y=j+1.2, "QTL-rich clusters", cex = 0.5 ,pos=4)
        }
        if (files[j,3]==2){
          rect(xleft = 0,
               #         #ybottom = j+0.1,
               ybottom = j+0.6-2,
               ytop = j+2.5-2,
               xright = len,
               #         #ytop = j+0.9,
               col = NA,
               border = "black"#color_pad[as.numeric(as.vector(chro1[i]))] 
          )
          segments(x0=0,x1 = len,y0 = j+1.25-2,y1 = j+1.25-2)
          segments(x0=0,x1 = len,y0 = j+1.9-2,y1 = j+1.9-2)
          text(x=len, y=j-0.45, "SGR with NongDa3338", cex = 0.5 ,pos=4)
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


