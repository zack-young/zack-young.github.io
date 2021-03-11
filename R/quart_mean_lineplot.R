library(mixtools)
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(grid))
library(reshape2)
present_sample <- read.table("/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/cultivar_list",as.is = T, header = F, comment.char = "")


tmp_index <- colnames(dt) %in% present_sample[[1]] # match return the position of x in table
dist_mat <- dt[,tmp_index ]
dist_mat <- dt
#order(unname(mat[1,]))
# 
# floor(length(use_num)*0.95)
# ceiling(length(use_num)*0.05)

d <- read_delim("/data/user/yangzz/mapping/fieldergenomecompare/statistic/201_10_nearest/combine_2_D_hete_density", "\t", escape_double = FALSE, col_names = F,trim_ws = TRUE, comment = "#")
ab <- read_delim("/data/user/yangzz/mapping/fieldergenomecompare/statistic/201_10_nearest/combine_2_AB_hete_density", "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE, comment = "#")
dt_all <- read_delim("/data/user/yangzz/mapping/fieldergenomecompare/1.diff_dev/201_final/combine.1M.combinediff", "\t", escape_double = FALSE, col_names = F,trim_ws = TRUE, comment = "#")
#a <- read_delim("/data/user/yangzz/mapping/fieldergenomecompare/2.sample_compare/combine/total_combine_1M", "\t", escape_double = FALSE, col_names = TRUE,trim_ws = TRUE, comment = "#")
nrows(filter(a,X4>=1))
logdt <- log(a[,4]+1,10)
a[is.na(a)]=0

a<-data.frame(d)

###
#mixmdla <- normalmixEM(log(a[,4]+1,10),k=4)

norm_dt<- a[which(a$X4 > 0),]
mixmdla <- normalmixEM(log10(norm_dt$X4),mean.constr=c(0.4,2.2,3.5),k=3) #mean.constr=c(0.55,NA,3.6)
mixmdlab <- mixmdla
#mixmdla <- normalmixEM(DF$V2,k=2)
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,10))
#plot.window(xlim=c(0,4), ylim=c(0,4))
plot(mixmdla, which =2, 
     breaks=100,
     #xlab2 = "Unmatch SNP per 1Mb",ylab2 = "Density",
     xlab2 = "Different SNPs Subtraction Qualified SNPs per 1Mb ",ylab2 = "Density",
     main2 = ""
     ,add = TRUE
)
lab=expression(log[10](frac(diff_snp,1-missing_rate)+1))
#log10((diff_snp/(1-missing_rate))+1)
axis(1,seq(0,4,0.5),cex.axis=1.5); axis(2,cex.axis=1.5)
title(xlab=lab, ylab="Density",cex.lab=1.3,line=5)
title(ylab="Density",cex.lab=1.5,line=3)
title(main="D")
R <- qnorm(0.97,mixmdla$mu[2], mixmdla$sigma[2])
par(mar=c(8, 5, 2, 2))
mu = c(mixmdlab$mu[1], mixmdlab$mu[2], 3.5360584)   #mu 0.6522212 2.1206051 sig 0.2687575 0.3733764 #0.6434359 2.1084451 3.5648158
sigma = c(mixmdlab$sigma[1], mixmdlab$sigma[2], 0.2342772)
#mu = c(0.5715048, mixmdla$mu[2], 3.5360584)
#sigma = c(mixmdla$sigma[1], mixmdla$sigma[2], 0.2342772)
plot(density(rnorm(100000, mean = mu[1], sd =sigma[1])),xlim=c(0,4),ylim=c(0,2.5),main="",xlab="",xaxt = "none",yaxt = "none")
lines(density(rnorm(100000, mean = mu[2], sd =sigma[2])))
axis(1,seq(0,4,0.5),cex.axis=1.5); axis(2,cex.axis=1.5); title(xlab="log10((diff_snp/(1-missing_rate))+1)", ylab="Density",cex.lab=1.5)
abline(v=log(18,10), col='red')



abline(v=qnorm(0.999,mean = 0.5715048, sd = 0.4662130))
abline(v=qnorm(0.001,mean = 2.1242645, sd = 0.2992025))
num = log(90,10)
area1 <- pnorm(num,mean = mu[2], sd = sigma[2]) 
area2 <- 1-pnorm(num,mean = mu[1], sd = sigma[1])
area2-area1  #1.38 (log(23))  <22 same  
#1.20 (log(16))  >48 same
abline(v=log(32,10), col='red')
abline(v=log(23,10), col=color[2])
abline(v=log(49,10), col=color[1])

for (i in seq(22,48)) {
  num=log(i+1,10)
  area1 <- pnorm(num, mean = 2.1242645, sd = 0.2992025) 
  area2 <- 1-pnorm(num,mean = 0.5715048, sd = 0.4662130)
  #cat(i,":[",area2,",",area1,"],")
  print(area1-9*area2)
}

##draw bayes method
dt_1 <- as.data.frame(density(rnorm(100000, mean = 2.1242645, sd = 0.2992025))$x)
dt_1$V2 <- density(rnorm(100000, mean = 2.1242645, sd = 0.2992025))$y
colnames(dt_1) <- c("V1","V2")
dt_tmp=dt_1[which(dt_1$V1<=log(33,10)), ]

dt <- as.data.frame(density(rnorm(100000, mean = 0.5715048, sd = 0.4662130))$x)
dt$V2 <- density(rnorm(100000, mean = 0.5715048, sd = 0.4662130))$y
colnames(dt) <- c("V1","V2")
dt_tmp_1=dt[which(dt$V1>=log(33,10)), ]
#plot(density(rnorm(100000, mean = 0.5715048, sd = 0.4662130)),xlim=c(0,4),ylim=c(0,1.5),main="")
par(mar=c(3,1,1,1))
plot(x=dt$V1,y=dt$V2,xlim=c(0,3),ylim=c(0,2),type = 'l',xaxt = "n", yaxt = "n")
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)
lines(x=dt_1$V1,y=dt_1$V2)
polygon(x=c(dt_tmp_1[1,1],dt_tmp_1$V1,dt_tmp_1$V1[nrow(dt_tmp_1)]),
        y=c(0,dt_tmp_1$V2,0), col = color[2], border = color[2])

polygon(x=c(dt_tmp[1,1],dt_tmp$V1,dt_tmp$V1[nrow(dt_tmp)]),
        y=c(0,dt_tmp$V2,0), col = color[1], border = color[1])
abline(v=log(23,10), col=color[2])
abline(v=log(49,10), col=color[1])
abline(v=log(33,10), col="black")

legend("topright",c("2.5%","97.5%",'Average'),lty = 1, cex = 0.6,col=c("blue","red",'black'),bty='n')

####

plot(mixmdla$x,mixmdla$posterior[,"comp.1"],type='l')


###

dt_tmp <- data.frame(V1=numeric(), stringsAsFactors=FALSE)

len <- length(all_hist$breaks[-1])

sp=spline(dt_tmp[,1],dt_tmp[,2],n=1000)
plot(x=dt_tmp[,1],y=log(dt_tmp[,2]+1,10),type='l',ylim=c(0,0.002))
sp=spline(dt_tmp[,1],dt_tmp[,3],n=1000)
lines(x=dt_tmp[,1],y=log(dt_tmp[,3]+1,10))
lines(x=dt_tmp[,1],y=log(dt_tmp[,4]+1,10))

par(mar=c(6, 6, 2, 2))
pdf("/data/user/yangzz/pdf/diff_snp_frequency.pdf", height = 4, width = 5)

sp=spline(dt_tmp[,1],dt_tmp[,2],n=500)
plot(x=sp,type='l',ylim=c(0,0.005),xlab = 'log10(diff_homoSNP)',ylab = 'Frequency')
sp=spline(dt_tmp[,1],dt_tmp[,3],n=150)
lines(x=sp,col='red')
sp=spline(dt_tmp[,1],dt_tmp[,4],n=500)
lines(x=sp,col='blue')
legend("topright",c("2.5%","97.5%",'Average'),lty = 1, cex = 0.6,col=c("blue","red",'black'),bty='n')
dev.off()
###
library('readr')
dt <- read_delim("/data/user/yangzz/mapping/fieldergenomecompare/statistic/201_10_nearest/combine_AB_col_density", "\t", escape_double = FALSE,col_names = F, trim_ws = TRUE)
dt <- data.frame(dt)

dt_use <- data.frame(a=c(seq(1,nrow(dt))))

lis_max=c()

for (i in seq(1,ncol(dt))) {
  lis_max=c(max(dt[,i]),lis_max)
  lis_max=c(min(dt[,i]),lis_max)
}
min(lis_max)

#draw diffsnp distribution in multiple hist
pdf("/data/user/yangzz/pdf/diff_snp_multi_lines.pdf", height = 4, width = 5)

### use hist to plot homo snp log10 
se <- seq(from = 0, to = 4, by=0.04)
start <- hist(log(na.omit(dt[,4])+1,10),breaks = se,plot = F)$mids
dat <- data.frame(start1=start)


for (i in seq(1,ncol(dt))) {
  dt_lim <- data.frame(dt[,i])
  dt_lim<- dt_lim[which(log(na.omit(dt_lim[,1])+1,10) < 4),]
  dat[,i+1] <- hist(log(na.omit(dt_lim)+1,10),breaks = se,plot = F)$density
  #hist(log(dt[,i]+1,10),breaks = 100,freq=F,ylim=c(0,5),add=T)
}
#lines(density(log(na.omit(dt[,4])+1,10)),col='red')
dat[,ncol(dt)+2] <- seq(1,nrow(dat))
dat[,ncol(dt)+3] <- seq(1,nrow(dat))
dat[,ncol(dt)+4] <- seq(1,nrow(dat))
for (i in seq(1,nrow(dat))) {
  dat[i,ncol(dt)+2] <- sum(sort(dat[i,2:(ncol(dt)+1)])[1:9])/9
  dat[i,ncol(dt)+3] <- sum(sort(dat[i,2:(ncol(dt)+1)])[(ncol(dt)-8):ncol(dt)])/9
  dat[i,ncol(dt)+4] <-  mean(as.numeric(sort(dat[i,2:(ncol(dt)+1)])))
}

hist(log(na.omit(dt[,3])+1,10),breaks = 100,freq=F,ylim=c(0,15),border = '#595959',xlab='',ylab = "",main='',yaxt='n',xaxt='n')
title(xlab="log10((diff_snp/(1-missing_rate))+1)", ylab="Density",cex.lab=1.5)
axis(1,seq(0,4,0.5),cex.axis=1.5);axis(2,cex.axis=1.5)
for (i in seq(1,ncol(dt))) {
  dt_lim <- data.frame(dt[,i])
  dt_lim<- dt_lim[which(log(na.omit(dt_lim[,1])+1,10) < 4),]
  hist(log(na.omit(dt_lim)+1,10),breaks = se,freq=F,ylim=c(0,5),border = '#595959',add=T)
}
lines(spline(x=dat[,1],y=dat[,ncol(dt)+4],n=100),col='blue')

sp_down = spline(x=dat[,1],y=dat[,ncol(dt)+2],n=100)

sp_up = spline(x=dat[,1],y=dat[,ncol(dt)+3],n=100)


polygon(x=c(sp_up$x,rev(sp_up$x)),y=c(sp_up$y,rev(sp_down$y)),col = rgb(255, 0, 0, 20, maxColorValue=255),border='red')

### rainbow plot
hist(as.numeric(dat[1,2:(ncol(dt)+1)]),plot=F)


###

### use density to plot
dt_tmp <- data.frame(density(log(na.omit(dt[,4])+1,10),bw=0.1)$x)
for (i in seq(1,length(colnames(dt)))) {
  dt_tmp[,i+1] <- density(log(na.omit(dt[,i])+1,10),bw=0.1)$y
  #hist(log(dt[,i]+1,10),breaks = 20,freq=F,ylim=c(0,4),add=T)
}
for (i in seq(1,nrow(dt_tmp))) {
  dt_tmp[i,356] <- sum(sort(dt_tmp[i,2:355])[1:9])/9
  dt_tmp[i,357] <- sum(sort(dt_tmp[i,2:355])[346:354])/9
  dt_tmp[i,358] <-  mean(as.numeric(sort(dt_tmp[i,2:355])))
}
plot(x=dt_tmp[,1],y=dt_tmp[,358],type='l',ylim=c(0,1.5),xlim=c(0,4),pch=20,xlab='log10(diff_homoSNP)',ylab = "Density",lwd=2)

polygon(x=c(dt_tmp[,1],rev(dt_tmp[,1])),y=c(dt_tmp[,357],rev(dt_tmp[,356])),col = rgb(255, 0, 0, 20, maxColorValue=255),border='red',lwd=2)

lines(x=dt_tmp[,1],y=dt_tmp[,357])
lines(x=dt_tmp[,1],y=dt_tmp[,356])
###


hist(log(dt[,1]+1,10),breaks = 70,freq=F,ylim=c(0,30))
plot(density(dt[,1]),ylim=c(0,30))
for (i in seq(1,length(colnames(dt)))) {
  lines(density(dt[,i]))
  #hist(dt[,i],breaks = 70,freq=F,ylim=c(0,4),add=T)
}
lines(x=mixmdla$x,y=mixmdla$posterior[,2],cols='green')

legend("topright",c("CS","Cultivar"),lty = c(1,1), cex = 0.2,col=c("red","black"),bty="n")
axis(1,seq(0, 2, 0.1), las = 1)


dt <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/20200416_201_CNV/chr4D_rht1_DP_18_20",as.is = T, header = F, comment.char = "")
dt1 <- read_delim("/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/cultivar_combine_norm", "\t", escape_double = FALSE, trim_ws = TRUE)

dt <- as.data.frame(dt)
#dt <- dt[-1]
### use hist to plot CNV norm
se <- seq(from = 0, to = 4, by=0.04)
start <- hist(dt[,1],breaks = se,plot = F)$breaks
dat <- data.frame(start=start[-length(start)], end=start[-1])
for (i in seq(1,ncol(dt))) {
  dat[,i+2] <- hist(dt[,i],breaks = se,plot = F)$density
  #hist(log(dt[,i]+1,10),breaks = 20,freq=F,ylim=c(0,4),add=T)
}
for (i in seq(1,nrow(dt))) {
  dt[i,204] <- mean(as.numeric(dt[i,4:203]))
  dt[i,205] <- sd(as.numeric(dt[i,4:203]))+mean(as.numeric(dt[i,4:203]))
  dt[i,206] <- mean(as.numeric(dt[i,4:203]))-sd(as.numeric(dt[i,4:203]))
  
  #hist(log(dt[,i]+1,10),breaks = 20,freq=F,ylim=c(0,4),add=T)
}
hete_lis_0.5 <- unique(sort(as.integer(as.numeric(lis)/1000)))*1000+1
colnames(dt)[2] <- 'start'
hete_df <- data.frame(x=hete_lis,y=rep(2,length(hete_lis)))

for (i in hete_lis_0.5){
  lis_0.5 <- c(lis_0.5,dt[which(dt$start == i),"V204"])
}

ggplot(dt)+ geom_ribbon(aes(start,ymin=V206,ymax=V205), fill="grey")+
  geom_line(aes(start,V206), col="grey",size=0.75)+
  geom_line(aes(start,V205), col="grey",size=0.75)+
  geom_line(aes(start,V204), col="black",size=0.75)+
  geom_vline(xintercept = 18781062,col="deepskyblue")+
  geom_point(aes(x=x,y=y),data = hete_df)
  #geom_hline(yintercept = 1.3,col="red")+
  cowplot::theme_cowplot()+  
  scale_y_continuous(limits=c(-7,10),breaks = seq(-5,5,1))+
  #scale_x_continuous(limits=c(18500000,19000000))+
  xlab("Position")+ylab("Read depth")


dat[,ncol(dt)+3] <- seq(1,nrow(dat))
dat[,ncol(dt)+4] <- seq(1,nrow(dat))
dat[,ncol(dt)+5] <- seq(1,nrow(dat))
for (i in seq(1,nrow(dat))) {
  dat[i,ncol(dt)+3] <- sum(sort(dat[i,3:(ncol(dt)+2)])[1:9])/9
  dat[i,ncol(dt)+4] <- sum(sort(dat[i,3:(ncol(dt)+2)])[(ncol(dt)-8):ncol(dt)])/9
  dat[i,ncol(dt)+5] <-  mean(as.numeric(sort(dat[i,3:(ncol(dt)+2)])))
}
sp_mean = spline(x=(dat[,1]+dat[,2])/2,y=dat[,ncol(dt)+5])
plot(sp_mean,type='l',ylim=c(0,20),xlim=c(0,2), pch=20,xlab='log10(diff_homoSNP)',ylab = "Density")
sp_down = spline(x=(dat[,1]+dat[,2])/2,y=dat[,ncol(dt)+3],n=100)
lines(sp_up)
sp_up = spline(x=(dat[,1]+dat[,2])/2,y=dat[,ncol(dt)+4],n=100)
lines(sp_down)

polygon(x=c(sp_up$x,rev(sp_up$x)),y=c(sp_up$y,rev(sp_down$y)),col = rgb(255, 0, 0, 20, maxColorValue=255),border='red')

lines(density(dt[,'CS'])$x,density(dt[,'CS'])$y,col='blue')

###

### use CNV density to plot
dt_tmp <- data.frame(density(dt[,4],bw=0.1)$x)
for (i in seq(1,length(colnames(dt)))) {
  dt_tmp[,i+1] <- density(dt[,i],bw=0.1)$y
  #hist(log(dt[,i]+1,10),breaks = 20,freq=F,ylim=c(0,4),add=T)
}
for (i in seq(1,nrow(dt_tmp))) {
  dt_tmp[i,ncol(dt)+2] <- sum(sort(dt_tmp[i,2:(ncol(dt)+1)])[1:9])/9
  dt_tmp[i,ncol(dt)+3] <- sum(sort(dt_tmp[i,2:(ncol(dt)+1)])[(ncol(dt)-8):ncol(dt)])/9
  dt_tmp[i,ncol(dt)+4] <-  mean(as.numeric(sort(dt_tmp[i,2:(ncol(dt)+1)])))
}
plot(x=dt_tmp[,1],y=dt_tmp[,ncol(dt)+4],type='l',ylim=c(0,20),xlim=c(0,4),pch=20,xlab='log10(diff_homoSNP)',ylab = "Density",lwd=2)

polygon(x=c(dt_tmp[,1],rev(dt_tmp[,1])),y=c(dt_tmp[,ncol(dt)+2],rev(dt_tmp[,ncol(dt)+3])),col = rgb(255, 0, 0, 20, maxColorValue=255),border='red',lwd=2)

lines(x=dt_tmp[,1],y=dt_tmp[,357])
lines(x=dt_tmp[,1],y=dt_tmp[,356])
###
dat[i,ncol(dt)+3] <- sum(sort(dat[i,3:(ncol(dt)+2)])[1:9])/9
dat[i,ncol(dt)+4] <- sum(sort(dat[i,3:(ncol(dt)+2)])[(ncol(dt)-8):ncol(dt)])/9
dat[i,ncol(dt)+5] <-  mean(as.numeric(sort(dat[i,3:ncol(dt)+2])))

