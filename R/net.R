library(netview)
library(igraph)
ad_file<-read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/distance_matrix/All_torestdis",stringsAsFactors = F)
# meta_data <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_196_nongke_10+.txt", header = F, sep = "\t", stringsAsFactors = F)
# meta_CN <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/distance_matrix/CN_zone.txt", header = T, sep = "\t", stringsAsFactors = F)
# colnames(meta_CN) <- c('V3','V5')
# meta_all <- merge(meta_CN,meta_data,by = 'V3',all = T)
#meta_all[is.na(meta_all)]='Fo'

meta_1 <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/distance_matrix/all_C_L.txt", header = F, sep = "\t", stringsAsFactors = F)
meta_2 <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/distance_matrix/all_region.txt", header = F, sep = "\t", stringsAsFactors = F)
meta_3 <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/distance_matrix/CN_zone.txt", header = F, sep = "\t", stringsAsFactors = F)
meta_4 <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/distance_matrix/CNL_zone2.txt", header = F, sep = "\t", stringsAsFactors = F)
meta_all <- merge(meta_1,meta_2,by = 'V1',all = T)
meta_all <- merge(meta_all,meta_3,by = 'V1',all = T)
colnames(meta_all) <- c('V1','V2','V3','V4')
meta_all <- merge(meta_all,meta_4,by = 'V1',all = T)
colnames(meta_all) <- c('V1','V2','V3','V4','V5')
#E_c<-get.adjacency(ad_file,sparse=FALSE)
ma_name <- unique(c(ad_file[,2],ad_file[,1]))
mat <- matrix(0, length(ma_name), length(ma_name),dimnames = list(ma_name,ma_name))
for (i in 1:nrow(ad_file)){
  mat[ad_file[i,1],ad_file[i,2]] <- ad_file[i,3]
  mat[ad_file[i,2],ad_file[i,1]] <- ad_file[i,3]
}
diag(mat) <- 1
mat_CN <- mat[meta_CN$V3[1:10],meta_CN$V3[1:10]]
col_anno <- data.frame(Track2=meta_all$V2,Track1=meta_all$V3,row.names  = meta_all$V1)
ann_colors <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/distance_matrix/coloryzz_6cols.txt",comment.char = ":", header = T, sep = "\t", stringsAsFactors = F)
col_lis_all=c(ann_colors$color)
names(col_lis_all) <-  unique(ann_colors$term)  #
col_lis_CN=c("#CE4444","#F88905","#F5F805","#F80680")
names(col_lis_CN) <- c("HH","BB","XN","CJ")
col_lis_CLW=c('#A8A8A8','#494848','#F83333')
names(col_lis_CLW) <- c('C',"L","W")
col_lis_CNL=c('#a6cee3',
  '#1f78b4',
  '#b2df8a',
  '#33a02c',
  '#fb9a99',
  '#e31a1c',
  '#fdbf6f',
  '#ff7f00',
  '#cab2d6',
  '#6a3d9a',
  '#ffff99',
  '#b15928')
names(col_lis_CNL)<-c('BJ','GS','HN',
                    'NX','QH','SC',
                    'SD','SSX','SX',
                      'XJ',
                      'XZ',
                      'YN')
col_lis <- list(Track1=col_lis_all,Track2=col_lis_CLW)

DIS=as.dist(1- mat) 
DIS=as.dist(log(1/mat,2))#dist(log(1/(COR+1),2))
DIS=as.dist((1- mat)^2) 
pheatmap::pheatmap(mat,display_numbers = F ,clustering_method ="ward.D2",   #"ward.D2",
                   annotation_col = col_anno,
                   clustering_distance_rows=DIS,clustering_distance_cols=DIS,annotation_legend = T,
                   annotation_row = col_anno, treeheight_row=100,treeheight_col=100,
                   fontsize_row = 5,fontsize_col = 5,annotation_colors = col_lis,main="1-mat")
#breaks=c(seq(0,0.25,0.05),0.5,0.9,1),legend_breaks=c(seq(0,0.25,0.05),0.5,1),color = c(colorRampPalette(c("#088CDD","#FDF002"))(6),"red","black"),
hc <- hclust(DIS,method = "ward.D2")
dendrogram <- as.dendrogram(hc,hang = -1)
plot(dendrogram, horiz=F)

####++++++++which(mat>0.9&mat<1,arr.ind=T)

voles.mds=cmdscale(1-mat,k=13,eig=T)
sum(abs(voles.mds$eig[1:3]))/sum(abs(voles.mds$eig))
sum((voles.mds$eig[1:2])^2)/sum((voles.mds$eig)^2)
x = voles.mds$points[,1]
y = voles.mds$points[,2]
p=ggplot(data.frame(x,y),aes(x,y,label = colnames(mat)))
p+geom_point(shape=16,size=3,colour='red')+
  geom_text(hjust=-0.1,vjust=0.5,alpha=0.5)
#####++++++++++++++++++
CHR="whole_indel"
sample_group <- read_delim("sample_group.txt", 
                           "\t", escape_double = FALSE, col_names = T, 
                           trim_ws = TRUE)

SM =  read.table( paste0(CHR, ".cluster2") )
TAB = read.table( paste0(CHR, ".mdist") )
colnames(TAB) = SM[,1]
rownames(TAB) = SM[,1]

tmp_index <- rownames(TAB) %in% sample_group$Name
TAB <- TAB[tmp_index, tmp_index]

sample_present <- sample_group[sample_group$Name %in% rownames(TAB),]
LIST <- sample_present$label
names(LIST) <- sample_present$Name

tmp <- dlply(sample_present, "Group")
cls <- lapply(tmp, `[[`, 3)

rownames(TAB) = LIST[rownames(TAB)]

ans <- read_table2("ans.txt", col_names = FALSE)
sample_group <- read_delim("sample_group.txt", 
                           "\t", escape_double = FALSE, col_names = T, 
                           trim_ws = TRUE)
tmp_index <- SM$V1 %in% sample_group$Name
ans <- ans[tmp_index,]

sample_present <- sample_group[sample_group$Name %in% SM$V1,]

metaD <- sample_present[,c(1,4)]
colnames(metaD)[1] <- "ID"
metaD$Colour = metaD$Group
metaD$Colour <- factor(metaD$Colour)

metaD <- read_delim("metaD.csv", "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

oysterOptions <- netviewOptions()
graphs <- netview(as.matrix(TAB), metaD, k=1:50, options=oysterOptions)

for (i in seq(from=5, to = 50, by=5)){
  pdf(file=paste("netview","k",i,".pdf", sep = "_"))
  plot(graphs[[paste0("k",i)]], vertex.size=3, vertex.label.color="grey", vertex.label.cex=0.3)
  dev.off()
}