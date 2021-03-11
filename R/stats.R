library(mixtools)
library(dplyr)
require(qcc)
library(ggplot2)

for(i in c('lx99','jm22')){
  for( j in c('snp','indel')){
    a <- paste("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/",i,"_snp_indel_analysis/",j,"/combine_",j,"_miss_DP",sep="")

    mydata2<-read.table(a,header=F,sep="\t")
    b <- paste("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/",i,"_snp_indel_analysis/",j,"/combine_",j,"_homo_DP",sep="")
    mydata<-read.table(b,header=F,sep="\t")
    c <- paste("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/",i,"_snp_indel_analysis/",j,"/combine_",j,"_hete_DP",sep="")
    mydata1<-read.table(c,header=F,sep="\t")
    x <- mydata1$V2 / sum(mydata1$V2)
    y <- mydata$V2 / sum(mydata$V2)
    z <- mydata2$V2 / sum(mydata2$V2)
    data2<-cbind(mydata1,x)
    data1<-cbind(mydata,y)
    data3<-cbind(mydata2,z)
    
    par(mar=c(5, 5, 5, 8),xpd=F)
    plot(x = data2[,1],y = data2[,3],'l',
         main="hete",xlab = "genotype quality",ylab = 'Density'
         ,ylim=c(0,0.2)
         ,xlim=c(0,100)
         ,col=('blue')
         #,add=TRUE
    ) 
    par(new=TRUE)
    plot(x = data1[,1],y = data1[,3],'l',
         main="hete",xlab = "genotype quality",ylab = 'Density'
         ,ylim=c(0,0.2)
         ,xlim=c(0,100)
         ,col=('red')
         #,add=TRUE
    ) 
    par(new=TRUE)
    plot(x = data3[,1],y = data3[,3],'l',
         main="hete",xlab = "genotype quality",ylab = 'Density'
         ,ylim=c(0,0.2)
         ,xlim=c(0,100)
         ,col=('green')
         #,add=TRUE
    )
    axis(1,seq(0,100,5))
    par(xpd=T)
    legend(106,0.1, bty='n',lty=c(1,1,1),c("Hete", "Homo","Miss"),col=c('blue','red',"green")
           #, fill = c("yellow", "orange")
    )
    axis(1,seq(0,20,5))
  }

}


files<-read.table("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/jm22_snp_indel_analysis/indel/combine_indel_altHOMO",
                  header=T,sep="\t")

alt=as.numeric(as.character(files[which(files$ALL_HOMO != 'ALL_HOMO'),]$ALT_HOMO)) #files[which(files$ALT_HOMO == 0),]
all=as.numeric(as.character(files[which(files$ALL_HOMO != 'ALL_HOMO'),]$ALL_HOMO))
#cleardata <- mydata[which(mydata[,1] != "."), ]

#newdata <- mydata[order(as.numeric(as.character(mydata$V1))),]
#a <- cleardata$V1
hist(log((alt/all),10),breaks = 100,main="lx99 snp alt homo ratio",xlab = "log10(alt_homo/all_homo)"
     ,xlim = c(-3,0)
     #,ylim = c(0,200)
)
axis(1,seq(-2.75,-0.25,0.5))
sumb = 0
suma = 0
for(i in 1:nrow(newdata)){
  a <- newdata[i,1]*newdata[i,2]
  b <- newdata[i,2]
  suma = suma + a
  sumb = sumb + b
}
plot(x = newdata$V1,
     'l',
     #lty = 3,
     y = newdata$V2,
     #xlim = c(0,25),
     xlab = "Gene Quality", 
     #xlab = "SNP Deepth", 
     ylab = "Counts", main = "F1")

Rbind <- suma/sumb + 1.96*sqrt(suma/sumb)
Lbind <- suma/sumb - 1.96*sqrt(suma/sumb)
abline(v=Rbind)
abline(v=Lbind)
axis(1,seq(0,100,5))


exactPoiCI <- function (X, conf.level=0.95) {
  alpha = 1 - conf.level
  upper <- 0.5 * qchisq(1-alpha/2, 2*X+2)
  lower <- 0.5 * qchisq(alpha/2, 2*X)
  return(c(lower, upper))
}

axis(1,0:100)
usedata=rpois(226959767, 7.256369)
exp(confint(glm(usedata ~ 1, family=poisson)))
confint(x=usedata,sigma=7.256369,alpha=7.256369)
poisson.test(25,T=1,alternative = "t",conf.level = 0.95)
ratio <- result[estimate,value]
7.256369 + 1.96*sqrt(7.256369)
sum(mydata[seq(1,3),2])

(2078556+5509903)/226959767
dpois(x=15, lambda=7.256369)
poisson.test(7,347)

files<-read.table("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/tmp/NZ10_YM05/combine.NZ10_YM05_unmatchhomo_snp_density",header=F,sep="\t")

par(xpd=T)
hist(log(files[,4],10),breaks = 100,main="FuGou583 S742 unmatch homo_snp density",xlab = "log10(unmatch homo_snp counts)",prob=T
     ,axes=F
     #,xlim =c(0,2)
     #,ylim = c(0,10)
)
axis(1,seq(from=0, to=4, by=0.5))
axis(2,seq(from=0, to=2, by=0.5))
abline(v=0.5,lty=2)
lines(density(files[,1]),col="red")



max(files[,4],files[,5])
files1 = files %>%
  mutate(V7 = pmax(V4,V5))
ratio <- files1[,6]/files1[,7]
chr <- files[,1]
p1 <- ggplot(data = files)+
  geom_line(aes(y=ratio, x=files[,2]
                #,color=color
                ))
  #ylim(0,0.75)+ 
  #scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.5, 0.6,0.7, 0.75))+
  #geom_line(aes(y=files[,7], x=files[,3],color=color))+
  #xlab(c)
p1 + facet_grid(files[,1] ~ .
                #,scales = 'free_y'
                #,space = 'free'
)
#files1[which(files1[,6]/files1[,7] >1), ]
ratio <- files[,8]/files[,9]
a <- files[which(files$V1=="chr4B"),]
hist(files[,8]/files[,9],xaxt="n",
     breaks = 300,main="lx99 jm22 ",xlab = "Homo SNP Unmatch Ratio per 1Mbp"
     #,xlim = c(6.5,9)
     ,ylim = c(0,200)
     )
axis(1,seq(0,1,0.1))     
files<-read.table("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/lx99_jm22/combine_different",header=F,sep="\t")

filt_data<-files[which(log(files$V4,10)<=3),]
hist(log(filt_data$V4,10))
mixmdla <- normalmixEM(log(filt_data$V4+1,10),k=2)
plot.new()
plot.window(xlim=c(0,3), ylim=c(0,3))
plot(mixmdla, which =2, 
     breaks=100,
     #xlab2 = "Unmatch SNP per 1Mb",ylab2 = "Density",
     xlab2 = "SNP Number Subtraction per 1Mb ",ylab2 = "Density",
     main2 = "LX99 JM22",add=TRUE
     #add = TRUE
)
axis(1); axis(2); title(xlab="log10(diff_homoSNP_density)", ylab="Density")


files[which(files$V9==0),]
for (i in unique(files[,1])) {
  assign(paste(i,"lx99",sep=""),sum(files[which(files$V1==i),4]))
  assign(paste(i,"lx99d",sep=""),sum(files[which(files$V1==i),5]))
  assign(paste(i,"jm22",sep=""),sum(files[which(files$V1==i),6]))
  assign(paste(i,"jm22d",sep=""),sum(files[which(files$V1==i),7]))
  assign(paste(i,"diff",sep=""),sum(files[which(files$V1==i),8]))
}

student<-t(data.frame(
                    #ID=c(unique(files[,1])),
                    lx99_snp=c(chr1Alx99,chr1Blx99,chr1Dlx99,chr2Alx99,chr2Blx99,chr2Dlx99,chr3Alx99,
                                                chr3Blx99,chr3Dlx99,chr4Alx99,chr4Blx99,chr4Dlx99,chr5Alx99,chr5Blx99,
                                                chr5Dlx99,chr6Alx99,chr6Blx99,chr6Dlx99,chr7Alx99,chr7Blx99,chr7Dlx99),
                    jm22_snp=c(chr1Ajm22,chr1Bjm22,chr1Djm22,chr2Ajm22,chr2Bjm22,chr2Djm22,chr3Ajm22,
                             chr3Bjm22,chr3Djm22,chr4Ajm22,chr4Bjm22,chr4Djm22,chr5Ajm22,chr5Bjm22,
                             chr5Djm22,chr6Ajm22,chr6Bjm22,chr6Djm22,chr7Ajm22,chr7Bjm22,chr7Djm22),
                    diff=c(chr1Adiff,chr1Bdiff,chr1Ddiff,chr2Adiff,chr2Bdiff,chr2Ddiff,chr3Adiff,
                                chr3Bdiff,chr3Ddiff,chr4Adiff,chr4Bdiff,chr4Ddiff,chr5Adiff,chr5Bdiff,
                                chr5Ddiff,chr6Adiff,chr6Bdiff,chr6Ddiff,chr7Adiff,chr7Bdiff,chr7Ddiff)))
names(student) <- c("1A", "1B", "1C",
                    "2A", "2B", "2C",
                    "3A", "3B", "3C",
                    "4A", "4B", "4C",
                    "5A", "5B", "5C",
                    "6A", "6B", "6C",
                    "7A", "7B", "7C")
barplot(student,
        #a = c(student[2,1],student[3,1],student[4,1]),b = c(student[2,2],student[3,2],student[4,2]),
        # #names.arg = mydata[,1]，
        #xlim = c(0,20),
        xlab= "Chromosome",
        ylab = "Counts",
        main="lx99 and jm22", 
        #col = "blue",border = "black",
        beside = T,
        legend.text = TRUE,
        names.arg = c("1A", "1B", "1D",
                      "2A", "2B", "2D",
                      "3A", "3B", "3D",
                      "4A", "4B", "4D",
                      "5A", "5B", "5D",
                      "6A", "6B", "6D",
                      "7A", "7B", "7D")
)

ggplot(data=student,aes(x=c(1,2,3),y=c(student$Name,student$Gender,student$Birthdate)))+
  geom_bar(position='stack')

c(unique(files[,1]))
sum(files[which(files$V1=="chr1A"), 4])

     
high <- mydata[which(mydata[,5] >4300), ]
low <- mydata[which(mydata[,5] <1000), ]
mid <- mydata[1000 < which(mydata[,5] < 4300), ]

hist(high[,6]/high[,5],breaks = 100,main="high",xlab = "hete/homo",xlim = c(0,1),ylim = c(0,200))
hist(low[,6]/low[,5],breaks = 100,main="low",xlab = "hete/homo",xlim = c(0,1),ylim = c(0,200))
hist(mid[,6]/mid[,5],breaks = 100,main="mid",xlab = "hete/homo",xlim = c(0,1),ylim = c(0,200))

hist(mydata[,6]/mydata[,5],breaks = 100,main="main",xlab = "hete/homo",xlim = c(0,1),ylim = c(0,200))
hist(mydata[,5],breaks = 100,main="F1",xlab = "snp_counts")

mydata<-read.table(
  "/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/CNV/combine.norm",header=F,sep="\t")
hist(mydata[,1],breaks = 100,main="high",xlab = "hete/homo"
     #,xlim = c(0,1)
     #,ylim = c(0,200)
     )
hist(mydata[,2],breaks = 100,main="high",xlab = "hete/homo"
     #,xlim = c(0,1)
     #,ylim = c(0,200)
)
hist(mydata[,1]-mydata[,2],breaks = 100,main="high",xlab = "hete/homo"
     #,xlim = c(0,1)
     #,ylim = c(0,200)
)
p1 <- ggplot(data = mydata)+
  geom_histogram(aes(x=mydata[,1]-mydata[,2])
                 ,binwidth=0.005
                #,color=color
  )+
xlim(-0.1,0.1) 
#+scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.5, 0.6,0.7, 0.75))+
#geom_line(aes(y=files[,7], x=files[,3],color=color))+
#xlab(c)
p1 + facet_grid(mydata[,3] ~ .
                #,scales = 'free_y'
                #,space = 'free'
)



#text(x=seq(0,1,0.1), y=rep(0,10),
#     seq(0,1,0.1), cex = 1,pos = 1)

rows <- nrow (mydata)
plot(x = mydata[,1],y = mydata[,2],
     #log = "xy" ֱ?Ӷ?xy??ȡlog
     #xlim = c(0,2),ylim = c(0,2),
     pch=".", cex = 1.2, #pch = 20,
     #xlim=c(0,1000),  ylim=c(0,1000),
     #xlab = 'fs/fm',#log = "xy",
     ylab = 'ms/fm',col = 'red',main = '')
barplot(mydata[,1],
        names.arg = mydata[,1],
        #xlim = c(0,20),
        xlab= "unmatch snp/snp number subtraction",
        ylab = "counts",
        main="lx99 and jm22", 
        col = "blue",border = "blue"
)
plot(0,xlim = c(0,4), ylim = c (0,4))
abline(a = 1.3, b = 0.5)
abline(a = -1.9, b = 1.7)
abline(a = -0.3, b = 1.1)
