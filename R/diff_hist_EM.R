library(readr)
meta_in <- "/data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_final_headed.txt"
#infile <- "/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/combine.1M.join_norm"
infile <- "/data/user/yangzz/mapping/fieldergenomecompare/1.diff_dev/201_final/combine.1M.combinediff"
meta_data <- read_delim(meta_in, "\t", escape_double = FALSE, trim_ws = TRUE, comment = "#")
dt_all <- read_delim(infile, "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
data_hete_miss <- read_delim(infile, "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
data_only_miss <- read_delim(infile, "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
sample_list <- read_delim("/data/user/yangzz/mapping/fieldergenomecompare/1.diff_dev/200_2_nearest_line.txt", "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
## filter sample
lis <- sample(1:19900,1000,replace=F)
#tmp_index <- colnames(data) %in% meta_data$Sample
#data <- data[,tmp_index]
da_only_miss <- data_only_miss[,as.numeric(sample_list$X1)]
da <- dt_all[,as.numeric(sample_list$X1)]
da_hete_miss <- data_hete_miss[,as.numeric(sample_list$X1)]
da <- da[,-115] #S138_S137
da <- da[,-150] #S205_C23
da <- da[,-84] #S122_C29
da <- da[,-12] #C30_C29
da <- da[,-242] #S251_S105
da <- da[,-83] #S122_C30
######
da <- data[,lis]
se <- seq(from = 1, to = 5, by=0.04)


hist_dt <- hist(log(da[[1]]+10,10),freq = F,breaks =se ,ylim=c(0,10),xlim=c(1,4),border = rgb(0.19,0.19,0.19,alpha = 0.8),xlab='',ylab = "",main='')#,yaxt='n',xaxt='n')
dt_freq <- as.data.frame(hist_dt$counts)
x_list <- hist_dt$mids
title(xlab="log10((diff_snp/(1-missing_rate))+1)", ylab="Density",cex.lab=1.5)
axis(1,seq(0,2,0.5),cex.axis=1.5);axis(2,cex.axis=1.5)
num=1
for (i in colnames(da)) {
  #dt_lim <- data.frame(dt[,i])
  #dt_lim<- dt_lim[which(log(na.omit(dt_lim[,1])+1,10) < 4),]
  hist_dt <- hist(log(da[[i]]+10,10),breaks = se,freq=F,ylim=c(0,10),xlim=c(1,4),border = rgb(0.19,0.19,0.19,alpha = 0.1),add=T)
  #if (max(hist_dt$counts) <= 6000) {
  dt_freq[,num]= hist_dt$counts
  num = num+1
  #}
}
col_num <- ncol(dt_freq)
col_sum <- sum(dt_freq[[1]])
for (i in seq(1:nrow(dt_freq))) {
  row_num <- as.numeric(dt_freq[i,1:col_num])/col_sum
  #dt_freq[i,col_num+1] <- quantile(dt_freq[i,1:col_num])$`25%`
  #dt_freq[i,col_num+2] <- quantile(dt_freq[i,1:col_num])$`75%`
  dt_freq[i,col_num+1] <- mean(row_num)
  dt_freq[i,col_num+2] <- sd(row_num)
  dt_freq[i,col_num+3] <- mean(row_num) + sd(row_num)
  dt_freq[i,col_num+4] <- mean(row_num) - sd(row_num)
}
dt_freq[,col_num+5] <- x_list
names(dt_freq)[col_num+1] <- "mean"
names(dt_freq)[col_num+2] <- "sd"
names(dt_freq)[col_num+3] <- "upper"
names(dt_freq)[col_num+4] <- "lower"
names(dt_freq)[col_num+5] <- "pos"

#dt_freq_filter <- dt_freq[which(dt_freq[["mean"]]>=1), ]
dt_sp <- cbind(as.data.frame(spline(dt_freq_filter$pos,dt_freq_filter$upper,n=100))
                             ,as.data.frame(spline(dt_freq_filter$lower,dt_freq_filter$lower,n=100))[,2]
                             ,as.data.frame(spline(dt_freq_filter$pos,dt_freq_filter$mean,n=100))[,2]
                             )
names(dt_sp)<- c('pos','upper','lower','mean')
dt_freq[which(dt_freq[["lower"]]<0),"lower"]=0
dt_freq_filter[,"mean"]=dt_freq_filter[,"mean"]/sum(dt_freq_filter[,"mean"])
dt_freq_filter[,"upper"]=dt_freq_filter[,"upper"]/sum(dt_freq_filter[,"upper"])
dt_freq_filter[,"lower"]=dt_freq_filter[,"lower"]/sum(dt_freq_filter[,"lower"])
ggplot(dt_freq)+ geom_ribbon(aes(pos,ymin=lower,ymax=upper), fill="grey")+
  geom_line(aes(pos,lower), col="grey",size=0.75)+
  geom_line(aes(pos,upper), col="grey",size=0.75)+
  geom_line(aes(pos,mean), col="black",size=0.75)+
  cowplot::theme_cowplot()+
  # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
  #geom_line(aes(density, mean), size=1.2, col="red") +
  #geom_vline(xintercept = 30) +
  scale_x_continuous(limits=c(1,4),breaks = seq(1,4,0.5))+ #breaks=seq(0,800000000,400000000),labels = c("0(Mb)",400,800)
  #scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1)) +
  ylab("") + xlab("")
#####
dt_col_only_miss <- da_only_miss[,1]
names(dt_col_only_miss) <- 'V1'
for (i in seq(2,ncol(da_only_miss))){
  da_col_only_miss <- da_only_miss[,i]
  names(da_col_only_miss) <- 'V1'
  dt_col_only_miss <-rbind(dt_col_only_miss,da_col_only_miss)
}
dt_col <- da[,1]
names(dt_col) <- 'V1'
for (i in seq(2,ncol(da))){
  da_col <- da[,i]
  names(da_col) <- 'V1'
  dt_col <-rbind(dt_col,da_col)
}
dt_col_hete_miss <- da_hete_miss[,1]
names(dt_col_hete_miss) <- 'V1'
for (i in seq(2,ncol(da_hete_miss))){
  da_col_hete_miss <- da_hete_miss[,i]
  names(da_col_hete_miss) <- 'V1'
  dt_col_hete_miss <-rbind(dt_col_hete_miss,da_col_hete_miss)
}
hist_dt <- hist(log(dt_col[[1]]+10,10),freq = F,breaks =se ,ylim=c(0,3),xlim=c(1,4),border = rgb(0.19,0.19,0.19,alpha = 0.8),xlab='',ylab = "",main='')
hist_dt <- hist(log(dt_col_only_miss[[1]]+10,10),freq = F,breaks =se ,ylim=c(0,3),xlim=c(1,4),border = rgb(0.19,0.19,0.19,alpha = 0.8),xlab='',ylab = "",main='')
hist_dt <- hist(log(dt_col_hete_miss[[1]]+10,10),freq = F,breaks =se ,ylim=c(0,3),xlim=c(1,4),border = rgb(0.19,0.19,0.19,alpha = 0.8),xlab='',ylab = "",main='')


#mixmdla_2 <- normalmixEM(log(dt_col$V1+10,10),mean.constr=c(1.1,2.1,3.5),k=3) #mean.constr=c(0.55,NA,3.6)
mixmdla_2 <- normalmixEM(log(dt_col$V1+10,10),k=3) #mean.constr=c(0.55,NA,3.6)

#####
plot.new()
plot.window(xlim=c(1,5), ylim=c(0,1.5))
plot(mixmdla_2, which =2, 
     breaks=50,
     #xlab2 = "Unmatch SNP per 1Mb",ylab2 = "Density",
     xlab2 = "Different SNPs Subtraction Qualified SNPs per 1Mb ",ylab2 = "Density",
     main2 = "",nclass=20
     ,add = TRUE
)
axis(1,seq(1,4,0.5),cex.axis=1.5); axis(2,seq(0,1.5,0.5),cex.axis=1.5)
mu = c(mixmdla$mu[1], mixmdla$mu[2],mixmdla$mu[3])   #mu 0.6522212 2.1206051 sig 0.2687575 0.3733764 #0.6434359 2.1084451 3.5648158
sigma = c(mixmdla$sigma[1], mixmdla$sigma[2], mixmdla$sigma[3])
lambd <- c(mixmdla$lambda[1], mixmdla$lambda[2], mixmdla$lambda[3])
turn_label_log <- function(x) {paste0("10e-", 6-x)}
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}
p1 <- data.frame(x = mixmdla_2$x) %>% 
  ggplot() +
  geom_histogram(aes(x, ..density..), bins=50, colour = "black", fill = "grey90", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdla_2$mu[1], mixmdla_2$sigma[1], lam = mixmdla_2$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdla$mu[2], mixmdla_2$sigma[2], lam = mixmdla_2$lambda[2]),
                colour = "blue", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdla$mu[3], mixmdla$sigma[3], lam = mixmdla$lambda[3]),
                colour = "green", lwd = 1.5)+
  ylab("Density") + xlab("Variant density (per bp)") +
  cowplot::theme_cowplot()
p1


ta <- merge(as.data.frame(table(round(data[[1]]))),as.data.frame(table(round(data[[2]]))),all=T,by = 'Var1',sort = T)
for (i in seq(3:ncol(data))) {
  ta <- merge(ta,as.data.frame(table(round(data[[i]]))),all=T,by = 'Var1',sort = T)
}
hist(log(da[[1]]+1,10),breaks = 100,xlim = c(0,4),ylim = c(0,10000),border = rgb(0,0,0,alpha = 0.5))
for (i in seq(2,293)) {
  hist(log(da[[i]]+1,10),breaks = 100,xlim = c(0,4),ylim = c(0,10000),border = rgb(0,0,0,alpha = 0.5),add=T)
}
for (i in seq(10,30)) {
  num=log(i+10,10)
  area1 <- mixmdla_2$lambda[1]*pnorm(num, mean = mixmdla_2$mu[1], sd = mixmdla_2$sigma[1],lower.tail = F) 
  area2 <- mixmdla_2$lambda[1]*pnorm(num,mean = mixmdla_2$mu[2], sd = mixmdla_2$sigma[2])
  #cat(i,":[",area2,",",area1,"],")
  print(area1-area2)
}

