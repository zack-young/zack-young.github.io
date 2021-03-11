library(ggplot2)
dt <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/CNV_count/combine_mask_CNV_deletion",as.is = T, header = F, comment.char = "")
v <- paste(dt[,1],dt[,2],sep="_") 
#
dt_freq <- as.data.frame(table(v))
nrow(dt_freq[which(dt_freq[,2]>=1), ])
dt_tmp <- data.frame(V1=numeric(),V2=numeric(), stringsAsFactors=FALSE)
for (i in seq(1:max(dt_freq[,2]))) {
  dt_tmp[i,1]=i
  dt_tmp[i,2]=nrow(dt_freq[which(dt_freq[,2]>=i), ])/14066   #14Gbp
  dt_tmp[i,3]="del"
}
#/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/combine_mask_CNV_duplication_split
dt <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/CNV_count/combine_mask_CNV_duplication",as.is = T, header = F, comment.char = "")
v <- paste(dt[,1],dt[,2],sep="_") 
#
dt_freq <- as.data.frame(table(v))
nrow(dt_freq[which(dt_freq[,2]>=1), ])
dt_tmp_1 <- data.frame(V1=numeric(),V2=numeric(), stringsAsFactors=FALSE)
for (i in seq(1:max(dt_freq[,2]))) {
  dt_tmp_1[i,1]=i
  dt_tmp_1[i,2]=nrow(dt_freq[which(dt_freq[,2]>=i), ])/14066
  dt_tmp_1[i,3]= "dup"
}
dt_total <- rbind(dt_tmp,dt_tmp_1)
#dt_total$V3 <- factor(dt_total$V3)
pdf("/data/user/yangzz/pdf/whole_geno_CNV.pdf", height = 5, width = 7)

ggplot(dt_total,aes(V1,V2,fill=V3))+
  geom_area(position="identity")+
  scale_fill_manual(values = c( "#4889E6","#FF7F7F"))+
  cowplot::theme_cowplot()

plot(x=dt_tmp[,1],y=dt_tmp[,2], xlim =c(1,3),ylim = c(0,0.3) ,col="#4889E6",pch=19,cex=2,
     xlab = 'Number of samples sharing the CNV',ylab='CNV percentage of the genome')
points(x=dt_tmp_1[,1],y=dt_tmp_1[,2], xlim =c(1,3),col="#FF7F7F",pch=19,cex=2)
legend("topright",c("Deletion","Duplication"),lty = c(1,1), cex = 1,col=c("blue","red"),bty="n")
dev.off()
#

# uniqv <- unique(v)
# use_list <- tabulate(match(v, uniqv))
# dt_tmp <- data.frame(V1=numeric(),V2=numeric(), stringsAsFactors=FALSE)
# for (i in seq(1:max(use_list))) {
#   dt_tmp[i,1]=i
#   dt_tmp[i,2]=length(use_list[use_list==i])
# }
# 
# 
# o <- order(dt[,1],dt[,2]) #arrange data frame by the order of two column
# dt <- dt[o,]


dt <- read.table("/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/combine_mask_CNV_deletion_merged",as.is = T, header = F, comment.char = "")
dt <- unique(dt)
v <- dt[,3]-dt[,2]
uniqv <- sort(unique(v))
use_list <- tabulate(match(v, uniqv))
dt_tmp <- data.frame(V1=numeric(),V2=numeric(), stringsAsFactors=FALSE)
for (i in seq(1:length(use_list))) {
  dt_tmp[i,1]=uniqv[i]
  dt_tmp[i,2]=use_list[i]
}
par(mar=c(5.1,5.1,4.1,3.1), oma=c(0,0,0,0))
plot(x=dt_tmp[,1]/1000000,y=log(dt_tmp[,2],10), xlim =c(0,250),#ylim = c(0,3) ,xlim =c(0,200),
     #xaxt="n",
     pch=16,cex=0.5,xlab = 'CNV length',ylab='Number of CNV in genome')
axis(1,seq(0,200,10))
