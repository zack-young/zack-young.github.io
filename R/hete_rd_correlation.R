dt <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/20200416_201_CNV/chr4D_rht1_DP_18_20",as.is = T, header = F, comment.char = "")
dt <- as.data.frame(dt)
for (i in seq(1,nrow(dt))) {
  dt[i,204] <- mean(as.numeric(dt[i,4:203]))
}
colnames(dt)[2] <- 'start'
df_rd <- dt[,c("start","V204")]
rownames(df_rd) <- df_rd[,1]
colnames(df_rd) <- c("start","RD")

infile <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmp.Nf6HH4JtHA" 
tmp <- read.table(infile, header = T, check.names = F, sep = "\t", stringsAsFactors = F)
dt_hete <- data.frame()
num=1
for (i in seq(1,nrow(tmp))){
  dt_hete[i,1] <- as.integer(as.numeric(strsplit(tmp[i,1],split = ";")[[1]][1])/1000)*1000 + 1
  dt_hete[i,2] <- length(which(tmp[i,]==0.5))/ncol(tmp)
  
}
dt_hete_2 <- aggregate(dt_hete[2], dt_hete[-2], FUN = function(X) mean(X))
rownames(dt_hete_2) <- dt_hete_2[,1]
colnames(dt_hete_2) <- c("start","hete")

df_com <- merge(df_rd,dt_hete_2,all = T)
df_com[is.na(df_com)] <- 0
df_com_filter <- df_com[which(df_com[,'hete']!=0),]
plot(jitter(df_com_filter$RD),jitter(df_com_filter$hete),pch=20)
plot(df_com_filter$RD,df_com_filter$hete,pch=20,cex=3,col=rgb(0, 191, 255, alpha=80, maxColorValue = 255))
