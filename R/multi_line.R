#/bin/Rscript
dt_CNV <- read.table(paste("/data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/fielder/fielder_CNV_counts/chr","1","A",".1M.CNV_duplication_count_compensent",sep=""),as.is = T, header = F, comment.char = "")
dt_CNV_1 <- read.table(paste("/data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/fielder/fielder_CNV_counts/chr","1","A",".1M.CNV_deletion_count_compensent",sep=""),as.is = T, header = F, comment.char = "")

df1 <- read.table("CS_abs_gene_sum.txt", header = F, stringsAsFactors = F)
df2 <- read.table("Zang1817_abs_gene_sum.txt", header = F, stringsAsFactors = F)

df1$V1 <- df1$V1/max(df1$V1)
df2$V1 <- df2$V1/max(df2$V1)

sarm <- 20
larm <- 30
tarm <- sarm + larm

pdf("RD_quantile.pdf", width = 10, height = 6)
par(mgp = c(2,0,0))
plot(x=NULL,y=NULL,xlim=c(-10,800*3.5-200),ylim=c(0,32), frame.plot = F, yaxt = "n", xaxt = "n", xlab = "Physical Chromosome Position", ylab = "Relative CNV Density", cex.lab=1.5, main = "Relative CNV Density in Fielder and Cultivars", cex.main = 1.5)

# rect(xleft = seq(0,13.5,1.5), ybottom = rep(0,10)+0.03, xright = seq(1.5,15,1.5), ytop = df1$V1[11:20]+0.03, col = rgb(51, 126, 187, maxColorValue=255), border = F)
centromere <- c(215,239,173,337,346,268,300,346,241,300,317,184,254,206,189,286,327,211,357,287,340)

chro_list=c("A","B","D")
for (i in seq(0,6)){
  for (a in seq(0,2)){
    dt_CNV <- read.table(paste("/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/field_CNV_count/chr",i+1,chro_list[a+1],".combine_mask_CNV_deletion_compensent",sep=""),as.is = T, header = F, comment.char = "")
    dt_CNV_1 <- read.table(paste("/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/field_CNV_count/chr",i+1,chro_list[a+1],".combine_mask_CNV_duplication_compensent",sep=""),as.is = T, header = F, comment.char = "")
    list_num=c(6,5,4,3,2,1,0)
    
    #lines(x=dt_CNV[,1]/1000000+a*950,y=dt_CNV[,3]+0.05+i*4 ,col = "#EE8805")
    #lines(x=dt_CNV_1[,1]/1000000+a*950,y=dt_CNV_1[,3]+i*4 ,col = "#3CB7EF")
    polygon(x=c(a*950,dt_CNV[,1]/1000000+a*950,dt_CNV[nrow(dt_CNV),1]/1000000+a*950),y=c(1+0.05+list_num[i+1]*5, dt_CNV[,3]+0.05+list_num[i+1]*5, 1+0.05+list_num[i+1]*5),col = "#EE8805",border = "#EE8805")
    polygon(x=c(a*950, dt_CNV_1[,1]/1000000+a*950, dt_CNV_1[nrow(dt_CNV_1),1]/1000000+a*950),y=c(1+0.05+list_num[i+1]*5, dt_CNV_1[,3]+0.05+list_num[i+1]*5, 1+0.05+list_num[i+1]*5),col = "#3CB7EF",border = "#3CB7EF")
    #centrosome
    polygon( c(centromere[(i)*3+a+1]+a*950-10,centromere[(i)*3+a+1]+a*950,centromere[(i)*3+a+1]+a*950+10), c(list_num[i+1]*5-0.2, list_num[i+1]*5, list_num[i+1]*5-0.2), border = "black", col = "black")
    
    # yaxis
    segments(x0 = 0+a*950, 
             x1 = 0+a*950, 
             y0 = 0+list_num[i+1]*5,
             y1 = 2+list_num[i+1]*5,
             col = "black",
             xpd = T)
    # yaxis
    segments(x0 = 0+a*950,
             x1 = 800+a*950,
             y0 = 0+list_num[i+1]*5,
             y1 = 0+list_num[i+1]*5,
             xpd = T)
    # yaxis ticks
    ticky_y <- c(0, 1, 2)
    ticky_x <- c(0, 200,400,600, 800)
    segments(x0 = 0+a*950,
             x1 = -10+a*950,
             y0 = ticky_y+list_num[i+1]*5,
             y1 = ticky_y+list_num[i+1]*5,
             xpd = T)
    # xaxis ticks
    segments(x0 = ticky_x+a*950,
             x1 = ticky_x+a*950,
             y0 = 0+list_num[i+1]*5,
             y1 = -0.1+list_num[i+1]*5,
             xpd = T)
    # axis texts
    text(x = -0.1+a*950,
         y = ticky_y+list_num[i+1]*5,
         labels = c(0, 1, 2),
         cex = 0.7,
         xpd = T,
         pos = 2)
    text(x = ticky_x+a*950,
         y = 0+list_num[i+1]*5,
         labels = c("0Mbp","200Mbp","400Mbp","600Mbp","800Mbp"),
         cex = 0.7,
         xpd = T,
         pos = 1)
    text(x = 400+a*950,
         y = 0+list_num[i+1]*5-0.7,
         labels = paste("chr",i+1,chro_list[a+1],sep=""),
         cex = 1,
         xpd = T,
         pos = 1)
  }
}
speci_pos <- function(x,y,gene){
  list_num=c(6,5,4,3,2,1,0)
  segments(x0 = x+(y[2]-1)*950, 
           x1 = x+(y[2]-1)*950, 
           y0 = list_num[y[1]]*5,
           y1 = 2+list_num[y[1]]*5,
           col = "black",
           xpd = T)
  text(x=x+(y[2]-1)*950,
       y=list_num[y[1]]*5+2.8,
       c(gene), 
       cex = 0.5,
       adj = 0.5,
       col="black",
       srt = 45)
}
#
#pos_oper <- function(x){(x-1)*2+1}
#
ya <- c(1,1)#1A
speci_pos(508,ya,"Glu-A1")
speci_pos(4,ya,"Glu-A3")
speci_pos(381,ya,"SnRk2.3-1A")
#
ya <- c(1,2)#1B
speci_pos(555,ya,'Bx7OE;Bx13')
speci_pos(5,ya,'Glu-B3')
speci_pos(412,ya,'SnRk2.3-1B')
#
ya <- c(1,3)#1D
speci_pos(412,ya,'Glu-D1')
#
ya <- c(2,1)#2A
speci_pos(712,ya,'ppo-A1')
speci_pos(321,ya,'zds-A1')
speci_pos(121,ya,'SUS2-2A')
speci_pos(508,ya,'Cwi-A1')
speci_pos(740,ya,'Flo2-A1')
speci_pos(37,ya,'Ppd-A1')
#
ya <- c(2,2)#2B
speci_pos(171,ya,'SUS2-2B')
speci_pos(449,ya,'bas1-B1')
speci_pos(55,ya,'Ppd-B1')
#
ya <- c(2,3)#2D
speci_pos(119,ya,'SUS2-2D')
speci_pos(19,ya,'Rht8')
speci_pos(34,ya,'Ppd-D1')
#
ya <- c(3,1)#3A
speci_pos(448,ya,'lyce-3A')
speci_pos(59,ya,'TAR2.1-3A')
speci_pos(176,ya,'GS5-A1/3A')
speci_pos(532,ya,'TGW6-A1')
#
ya <- c(3,2)#3B
speci_pos(148,ya,'Talyc-B1')
#
ya <- c(3,3)#3D
speci_pos(106,ya,'CKX6-D1')
#
ya <- c(4,1)#4A
speci_pos(688,ya,'Wx-B1')
speci_pos(330,ya,'TGW6')
speci_pos(610,ya,'CWI-4A')
speci_pos(582,ya,'Rht-A1')
#
ya <- c(4,2)#4B
speci_pos(612,ya,'Tapds-B1')
speci_pos(31,ya,'Rht-B1')
speci_pos(657,ya,'Vrn2-4B')
#
ya <- c(4,3)#4D
speci_pos(18,ya,'Rht-D1')
speci_pos(509,ya,'Vrn2-4D')
#
ya <- c(5,1) #5A
speci_pos(587,ya,'Vrn1-5A')
speci_pos(650,ya,'Q')
speci_pos(698,ya,'Vrn2-5A')
speci_pos(645,ya,'NAC2-5A')
#
ya <- c(5,2) #5B
speci_pos(698,ya,'Vrn2-5A')
#
ya <- c(5,3) #5D
speci_pos(3,ya,'Pina-D1;Pinb-D1')
speci_pos(557,ya,'CWI-5D')
speci_pos(467,ya,'Vrn1-5D')
#
ya <- c(6,1) #6A
speci_pos(413,ya,'Rht24')
speci_pos(461,ya,'TPP-6AL1')
speci_pos(237,ya,'GW2-6A')
#
ya <- c(6,2) #6B
speci_pos(291,ya,'GW2-6B')
speci_pos(134,ya,'GPC-B1')
#
ya <- c(6,3)#6D
speci_pos(386,ya,'GS1a')
#
ya <- c(7,1)#7A
speci_pos(729,ya,'Psy-A1')
speci_pos(4,ya,'6-SFT-A2')
speci_pos(115,ya,'SUS1-7A')
speci_pos(170,ya,'GASR7-A1')
speci_pos(205,ya,'TGW-7A')
speci_pos(519,ya,'MOC1-7A')
speci_pos(585,ya,'SAP1-A1')
speci_pos(71,ya,'Vrn3-7A')
#
ya <- c(7,2)#7B
speci_pos(739,ya,'Psy-B1')
speci_pos(68,ya,'SUS1-7B')
speci_pos(9,ya,'Vrn3-7B')
#
ya <- c(7,3)#7D
speci_pos(109,ya,'SUS1-7D')
speci_pos(6,ya,'GS-D1')
speci_pos(68,ya,'Vrn3-7D')
# rect(xleft = seq(0,13.5,1.5), ybottom = rep(0,10)-0.03, xright = seq(1.5,15,1.5), ytop = -df2$V1[11:20]-0.03, col = rgb(228, 30, 38, maxColorValue=255), border = F)



dev.off()