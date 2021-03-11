library(mixtools)
options(scipen = 200)
j = "4"
b <- paste("D:/idea/data/DP/s",j,"_DP/combine_data",sep="")
c <- paste("chr",j,sep="")
#pdf(a, height = 4, width = 7)
#a  <- mydata[complete.cases(mydata),] # remove NA
# max(hist(a[,7])$density)
# lines(A, lty = 2, lwd = 2)
# plot(A)
#den <- as.data.frame(table(a[,7])) looking for density of specific data
for (i in seq(1,7,1)) {
  #i = 6
  a <- paste("D:/idea/data/S",j,"/chr",i,"A.ratio",sep="")
  mydata<-read.table(a,header=F,sep="\t")
  b <- paste("D:/idea/data/S",j,"/chr",i,"B.ratio",sep="")
  mydatb<-read.table(b,header=F,sep="\t")
  d <- paste("D:/idea/data/S",j,"/chr",i,"D.ratio",sep="")
  mydatd<-read.table(d,header=F,sep="\t")
  # dataa <- a[which(a$V1==paste("chr",i,"A",sep="")),]
  # datab <- a[which(a$V1==paste("chr",i,"B",sep="")),]
  # datad <- a[which(a$V1==paste("chr",i,"D",sep="")),]
  A <- density(mydata[,7])
  B <- density(mydatb[,7])
  D <- density(mydatd[,7])
  #plot(density(x),breaks=c(0,100),freq=T)
  #
  Namea <- paste("C:/Users/lenovo/Desktop/2019_1_15/S",j,"_",i,"A_hete_ratio.pdf",sep="")
  pdf(Namea, height = 9, width = 16)
  mixmdla <- normalmixEM(mydata[,7],k=2)
  plot.new()
  plot.window(xlim=c(0,
                     #0.2
                     max(A$x)
                     )
              , ylim=c(0,32))
  plot(mixmdla, which =2, 
       breaks=100,
       xlab2 = "SNP Average Deepth",ylab2 = "Density",main2 = paste("Density plot of S",j,"_",i,"A deepth",sep=""),
       add = TRUE)
  axis(1); axis(2); title(main=paste("Density plot of S",j,"_",i,"A_hete_ratio",sep=""), 
                          xlab="Hete/Homo", ylab="Density",
                          cex.main = 3,cex.lab=1.5)
  R <- qnorm(0.995,mixmdla$mu[1], mixmdla$sigma[1])
#  L <- qnorm(0.005,mixmdla$mu[1], mixmdla$sigma[1])
  abline(v=R)
#  abline(v=L)
  text(x = R+0.05, y = 20 , labels = paste("Confidence level(99.5%)=",format(R , digits = 3),sep=""))
#  text(x = R+0.01, y = 18 , labels = paste("Confidence level(0.5%)=",format(L , digits = 3),sep=""))
  dev.off()
  #
  Nameb <- paste("C:/Users/lenovo/Desktop/2019_1_15/S",j,"_",i,"B_hete_ratio.pdf",sep="")
  pdf(Nameb, height = 9, width = 16)
  mixmdlb <- normalmixEM(mydatb[,7],k=2)
  plot.new()
  plot.window(xlim=c(0,
                     #0.2
                     max(B$x)
                     ), ylim=c(0,32))
  plot(mixmdlb, which =2, 
       breaks=100,
       xlab2 = "SNP Average Deepth",ylab2 = "Density",main2 = paste("Density plot of S",j,"_",i,"B deepth",sep=""),
       add = TRUE)
  axis(1); axis(2); title(main=paste("Density plot of S",j,"_",i,"B_hete_ratio",sep=""), 
                          xlab="Hete/Homo", ylab="Density",
                          cex.main = 3,cex.lab=1.5)
  R <- qnorm(0.995,mixmdlb$mu[1], mixmdlb$sigma[1])
  abline(v=R)
  text(x = R+0.05, y = 20 , labels = paste("Confidence level(99.5%)=",format(R , digits = 3),sep=""))
  dev.off()
  #
  Named <- paste("C:/Users/lenovo/Desktop/2019_1_15/S",j,"_",i,"D_hete_ratio.pdf",sep="")
  pdf(Named, height = 9, width = 16)
  mixmdld <- normalmixEM(mydatd[,7],k=2)
  plot.new()
  plot.window(xlim=c(0,
                     #0.2
                     max(D$x)
                     ), ylim=c(0,32))
  plot(mixmdld, which =2, 
       breaks=100,
       xlab2 = "SNP Average Deepth",ylab2 = "Density",main2 = paste("Density plot of S",j,"_",i,"D deepth",sep=""),
       add = TRUE)
  axis(1); axis(2); title(main=paste("Density plot of S",j,"_",i,"D_hete_ratio",sep=""), 
                          xlab="Hete/Homo", ylab="Density",
                          cex.main = 3,cex.lab=1.5)
  R <- qnorm(0.995,mixmdld$mu[1], mixmdld$sigma[1])
  abline(v=R)
  text(x = R+0.05, y = 20 , labels = paste("Confidence level(99.5%)=",format(R , digits = 3),sep=""))
  dev.off()

}

