
dt <-dist(t(data),method = "euclidean")
# reorder sample index
data <- data[, order.dendrogram(dendrogram)]
par(mar=c(10, 8, 8, 12))
## empty figure panel
meta_data = read.table("~/mapping/fieldergenomecompare/metadata_cultivar_final.txt",sep = "\t",head = F) 
for (num in seq(1,7)) {
  for (char in c('A','B','D')) {
  name=paste(char,num,"ma",sep="")
  assign(name,read.csv(paste("~/mapping/fieldergenomecompare/statistic/chr",num,char,"_200_distance",sep=""),sep = "\t",head = T))
  }
}
Ama = (A1ma + A2ma+ A3ma+ A4ma+ A5ma+ A6ma+ A7ma)/7
Bma = (B1ma + B2ma+ B3ma+ B4ma+ B5ma+ B6ma+ B7ma)/7
Dma = (D1ma + D2ma+ D3ma+ D4ma+ D5ma+ D6ma+ D7ma)/7
ABma = (Ama + Bma)/2
ABDma = (Ama + Bma + Dma)/3

M2=ABDma
M2b=B4ma
M2d=D4ma
M2<-as.matrix(M2)


#
M2[M2 >= 0.5]=4
M2[(M2 >= 0.25)&(M2<0.5)]=3
M2[(M2 >= 0.125)&(M2<0.25)]=2
#M2[(M2 >= 0.0625)&(M2<0.125)]=1
#M2[M2<0.125]=0
diag(M2) <- 0  #0
M2[upper.tri(M2)] <- t(M2)[upper.tri(M2)]
# M2[(M2 > 0)&(M2 != 1)]=-log2(M2[(M2 > 0)&(M2 != 1)])
# M2[M2 ==0]=10
# M2[M2 ==1]=0
row.names(M2) = colnames(M2)
#tmp_index <- rownames(M2) %in% meta_data$V1
sample_present <- meta_data[meta_data$V1 %in% rownames(M2),]
colnames(sample_present) <- c("Sample","Name","Label","Group")
LIST <- sample_present$Label
names(LIST) <- sample_present$Sample #assign name to LIST
rownames(M2) <- LIST[rownames(M2)]
colnames(M2) <- LIST[colnames(M2)]
sample_list_b <- read.csv("~/mapping/fieldergenomecompare/statistic/Norin_4B_4D_linked.txt",sep="\t",head = F)
sample_list_b <- as.character(sample_list_b[[1]])



M3 <- M2[sample_list_b,sample_list_b] #round
lis <- c()
for (num in seq(1,ncol(M2))) {
  #M3[,num] <- as.integer(rank(M3[,num]))
  lis <- c(lis,names(M2[,num][M2[,num]==4]))
}
write.csv(lis,"~/mapping/fieldergenomecompare/statistic/sample_ratio4.csv",quote = F)
M4 <- M3 + t(M3)
write.csv(M3,"~/mapping/fieldergenomecompare/statistic/ABD_200_rht_sim_ratio.csv",quote = F)

library(stringr)
M3 <- M2[stringr::str_detect(string =rownames(M2),pattern = 'CNC' ),stringr::str_detect(string =rownames(M2),pattern = 'CNC' )]
M3_lw <- M3[upper.tri(M3)]

#
pdf("/data/user/yangzz/pdf/prolamin_6B_heatmap.pdf", height = 11, width = 15)
par(mar=c(10, 8, 8, 20), oma=c(0,4,2,4))
pheatmap(M2,color=colorRampPalette(c("#4E7DB8","#FEFEBD","#D9352A"))(1000),fontsize = 10,fontsize_row=5,fontsize_col=5,cellheight=5,cellwidth=5)
#plot(pheatmap(M2,color=colorRampPalette(c("#4E7DB8","#FEFEBD","#D9352A"))(1000),fontsize = 10,fontsize_row=3,fontsize_col=3,cellheight=5,cellwidth=5)$tree_row,hang=-1,cex=0.5)
dev.off()
#
#MDS
library(ggrepel)
color_list <- read.csv("~/mapping/fieldergenomecompare/coloryzz_6cols.txt",sep="\t",head = T)
color_lis <- color_list$color
names(color_lis) <- color_list$term
M2[M2 == 10]=100
M2[M2 == 1]=0
dt<-as.dist(M2)
mds_x = cmdscale(dt)#,eig=TRUE, x.ret=TRUE)
mds_x = data.frame(mds_x)
mds_list <- meta_data$V4
names(mds_list) <- meta_data$V3
xy = cbind(mds_x, X3=mds_list[as.character(rownames(M2))],X4=rownames(M2))
xy = cbind(xy, X5=color_lis[as.character(xy$X3)])

ggplot(xy)+geom_point(aes(x=X1,y=X2,color = factor(X3)))+
  geom_text_repel(aes(X1, X2, label=rownames(xy),color =factor(X3)))+
  scale_color_manual(breaks = color_list$term,values = as.character(color_list$color),
                     name = "MPG")#geom_text_repel(aes(X1, X2, label=rownames(xy)),)
pdf("/data/user/yangzz/pdf/D_sub_sim_tred_MDS.pdf", height = 11, width = 15)
plot(x=0,xlim=c(min(xy$X1)-0.5,max(xy$X1)+0.5), ylim=c(min(xy$X2)-0.5,max(xy$X2)+0.5),type="n")
points(x=xy$X1,y=xy$X2, cex = 2.5, col = c(as.character(xy$X4)),pch=19)
text(x=xy$X1,y=xy$X2,  labels =rownames(M2) ,cex = 0.5)
legend("topright", fill = colors, legend = names(colors), border = F, bty='n', cex = 1, x.intersp = 0,xpd=F)
dev.off()
#plot(dendrogram, horiz=T)#,leaflab="none",  dLeaf = 0 ,edgePar=list(lwd=8), yaxt = "none", yaxs="i")
#
#netview
library("netview")
library(DT)
library(networkD3)
library("igraph")
M2_data <- data.frame(ID =rownames(xy),Group= xy$X3, Colour=xy$X4)
dt<-as.dist(M2)
hc <- hclust(dt)
dendrogram <- as.dendrogram(hc)
phylo <- as.phylo(hc)
oysterOptions <- netviewOptions(selectionTitle="ABD_net", nodeID="ID", nodeGroup="Group", nodeColour="Colour", communityAlgorithms=c("Walktrap", "Infomap", "Fast-Greedy"))
graphs <- netview(M2, M2_data, k=1:60, cluster = TRUE, options=oysterOptions)
k10=graphs[[paste0("k",15)]]
plot(k10, vertex.size=3, vertex.label.color="grey", vertex.label.cex=0.3,mark.groups=communities(k10$walktrap))
for (i in seq(from=5, to = 50, by=5)){
  pdf(file=paste("/data/user/yangzz/pdf/ABD_netview","k",i,"nolabel.pdf", sep = "_"), height = 11, width = 15)
  k10=graphs[[paste0("k",i)]]
  plot(k10, vertex.size=3, vertex.label=NA,mark.groups=communities(k10$walktrap))
  
  #plot(k10, vertex.size=3, vertex.label.color="grey", vertex.label.cex=0.3,mark.groups=communities(k10$walktrap))
  legend("topright", fill = colors, legend = names(colors), border = F, bty='n', cex = 1, x.intersp = 0,xpd=F)
  dev.off()
}  

#igraph nwt
library(igraph)
M2=ABDma
M2<-as.matrix(M2)
M2[M2 >= 0.5]=4
M2[(M2 >= 0.125)&(M2<=0.25)]=2
M2[(M2 >= 0.25)&(M2<0.5)]=3
M2[(M2 > 0.0625)&(M2<=0.125)]=1
M2[M2<=0.0625]=0
diag(M2) <- 0  #0
M2[upper.tri(M2)] <- t(M2)[upper.tri(M2)]
row.names(M2) = colnames(M2)
#tmp_index <- rownames(M2) %in% meta_data$V1
sample_present <- meta_data[meta_data$V1 %in% rownames(M2),]
colnames(sample_present) <- c("Sample","Label","Name","Variety")
LIST <- sample_present$Name
names(LIST) <- sample_present$Sample #assign name to LIST
rownames(M2) <- LIST[rownames(M2)]
colnames(M2) <- LIST[colnames(M2)]
sample_present = cbind(sample_present, Colors=colors[as.character(sample_present$Variety)])
node <- data.frame(Id=sample_present$Name,Label=sample_present$Name,Country=sample_present$Variety,Color=sample_present$Colors)
write.csv(M2,"~/mapping/fieldergenomecompare/statistic/ABD_200_adjacent.csv")
write.csv(node,"~/mapping/fieldergenomecompare/statistic/ABD_200_node.csv")

library(MCL)
library(igraph)

mcl_dt <- mcl(x = M2, addLoops = TRUE, allow1 = F,ESM = T)$Equilibrium.state.matrix
gu <- graph.adjacency(mcl_dt, mode="undirected" )
plot(gu,layout=layout.fruchterman.reingold,vertex.size=4)

#adjacency <-mcl(x = adjacency, addLoops=TRUE, inflation = 1.01, max.iter = 100)
rownames(mcl_dt) <- rownames(M2)
colnames(mcl_dt) <- colnames(M2)
write.csv(mcl_dt,"~/mapping/fieldergenomecompare/statistic/4D_3M_MCL_e2r2.csv",quote = F)
write.table(gu,"~/mapping/fieldergenomecompare/statistic/4D_3M_MCL_e2r2.txt",quote = F)


a_test <-M2
a_test[upper.tri(a_test)] =0
a_test <- as.data.frame(a_test)
lis <- list()
for (i in colnames(a_test)){
  lis[[i]] <- rownames(a_test[which(a_test[[i]] > 2), ])
}
lis <- as.matrix(lis)
write.table(lis,"~/mapping/fieldergenomecompare/statistic/ABD_similar_3.csv")


g <- graph.adjacency(
  M2,
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)
gD <- simplify(g)
#
colors <- colorFile[,2]
names(colors) <- colorFile[,1]

NJ_use <- njs(dt)
DF_CNV <- full_join(as_tibble(NJ_use), sample_present, by = 'label')
hc <- hclust(dt)
dendrogram <- as.dendrogram(hc)
phylo <- as.phylo(hc)
col.vector <- vector(mode="character",length=nrow(phylo$edge))
n.tips <- length(phylo$tip.label)
col.vector[phylo$edge[,2]>n.tips] <- "black"
edge.data <- as.data.frame(phylo$edge)
for(i in seq_along(phylo$tip.label)){
  edge.row <- as.numeric(rownames(edge.data[edge.data$V2==i,]))
  col.vector[edge.row] <- colors[as.character(DF_CNV$country[match(phylo$tip.label[i], DF_CNV$label)])]
}
tip.color <- colors[as.character(DF_CNV$country)]
names(tip.color) <- as.character(DF_CNV$label)
plot(phylo, type = 'u',show.tip.label = T,tip.color=tip.color,
     edge.color = col.vector,cex = 0.5)
legend("bottom", fill = colors, legend = names(colors), border = F, box.col = "white", horiz=T, cex = 1, x.intersp = 0,xpd=T)
nodelabels() #add node numbers
tiplabels() #add tip number


col.vector <- vector(mode="character",length=nrow(NJ_use$edge))
n.tips <- length(NJ_use$tip.label)
col.vector[NJ_use$edge[,2]>n.tips] <- "black"
edge.data <- as.data.frame(NJ_use$edge)
for(i in seq_along(NJ_use$tip.label)){
  edge.row <- as.numeric(rownames(edge.data[edge.data$V2==i,]))
  col.vector[edge.row] <- colors[DF_CNV$Group[match(NJ_use$tip.label[i], DF_CNV$label)]]
}
tip.color <- colors[DF_CNV$Group]
names(tip.color) <- as.character(DF_CNV$label)
plot(NJ_use, type = "u", show.tip.label = T, edge.color = col.vector, cex=0.2)
legend("bottom", fill = colors, legend = names(colors), border = F, box.col = "white", horiz=T, cex = 1, x.intersp = 0,xpd=T)
  
#Norin10 Rht1

Rht8 <- read.csv("~/mapping/fieldergenomecompare/statistic/Rht8/chr2B_combine_homo_hmm_snp_level_rht8_dist",sep="\t",head = T)
rhtd <- read.csv("~/mapping/fieldergenomecompare/statistic/Norin10_RHTB1/chr4D_combine_homo_hete_miss_undefined_snp_level_0M_dist",sep="\t",head = T)
meta_data = read.table("~/mapping/fieldergenomecompare/metadata_cultivar_final_headed.txt",sep = "\t",head = T) 

#pro <- read.csv("~/mapping/fieldergenomecompare/statistic/prolamin/chr6B_prolamin_distance",sep="\t",head = T)
#sample_list <- read.csv("~/mapping/fieldergenomecompare/statistic/Norin_4B_4D_linked.txt",sep="\t",head = F)
sample_list <- as.character(sample_list[[1]])
Ma <- Rht8
Ma[Ma==11] <- 30
Ma[Ma==10] <- 0
Ma[Ma==11] <- 20
M2 <- (Ma)
M2<-as.matrix(M2)
row.names(M2) = colnames(M2)
M2[M2==8] <- 2

M2[M2==0] <- 10
M2[M2==11] <- 1
diag(M2) <- 0  #0
M2[upper.tri(M2)] <- t(M2)[upper.tri(M2)]

sample_present <- meta_data[meta_data$Sample %in% rownames(M2),]

#
colnames(sample_present) <- c("Sample","Label","Name","Variety")
LIST <- sample_present$Name
names(LIST) <- sample_present$Sample #assign name to LIST
rownames(M2) <- LIST[rownames(M2)]
colnames(M2) <- LIST[colnames(M2)]
M2 <- M2[sample_list,sample_list]
gu <- graph.adjacency(M2, mode="undirected" )
plot( gu ,layout=layout.fruchterman.reingold,vertex.size=3)
mcl(x = M2, addLoops = TRUE, allow1 = F,ESM = T)
mcl_dt <- mcl(x = M2, addLoops = TRUE ,ESM = TRUE)$Equilibrium.state.matrix
rownames(mcl_dt) <- rownames(M2)
colnames(mcl_dt) <- colnames(M2)
for (i in seq(nrow(mcl_dt))) {
  for (j in seq(nrow(mcl_dt))) {
    if (mcl_dt[i,j]>0) {
      mcl_dt[j,i] <- mcl_dt[j,i]
    }
  }
  
}

write.csv(M2,"~/mapping/fieldergenomecompare/statistic/chr2B_Rht8_undefined_dist.csv",quote = F)
write.csv(mcl_dt,"~/mapping/fieldergenomecompare/statistic/4D_rht1_0M_hmm_Norin_mcl.csv",quote = F)
#
colnames(sample_present) <- c("Sample","Name","country","label","Group")
rhtb_matrix <- M2
rhtd_matrix <- M2
mat <- rhtd_matrix 

sample_list <- rownames(mat[mat["EA_NongLin10",]>0,])
for (item_1 in sample_list) {
  for (item_2 in sample_list) {
    if (item_1 != item_2) {
    #print(paste(item_1,"_",item_2,sep=""))
      mat[item_1,item_2]=0
    }
  
  }
}





write.csv(M2,"~/mapping/fieldergenomecompare/statistic/Prolamin_6B.csv")

