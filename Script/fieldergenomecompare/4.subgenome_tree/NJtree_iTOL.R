#!/usr/bin/env Rscript

# generate file for iTOL

# libraries
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-d", "--distF"), dest = "distF", default = "",
              help = "dist matrix file, n rows by n cols."),
  make_option(c("-s", "--sampleF"), dest = "sampleF", default = "",
              help = "sample name that matches dist file."),
  make_option(c("-a", "--annoF"), dest = "annoF", default = "",
              help = "annotation file, which contains sample name to show, sample group."),
  make_option(c("-p", "--sampleP"), dest = "sampleP", default = "",
              help = "only sample name in this file will be plot."),
  make_option(c("-g", "--tmpG"), dest = "tmpG", default = NULL,
              help = "tmp group, if it exist, groups in this file will replace that in metadata."),
  make_option(c("-c", "--colorF"), dest = "colorF", default = "colorFile_clade.txt",
              help = "group name and color."),
  make_option(c("-t", "--title"), dest = "titleP", default = NULL,
              help = "group name and color."),
  make_option(c("-P", "--plot"), dest = "plotT", default = F,
              help = "plot tree on pdf."),
  make_option(c("-o", "--out"), dest = "out", default = "output",
              help = "output prefix.")
)

parser <- OptionParser(usage = "Rscript treeWholeNJ.R [options] Tip group files",
                       description = 
                       "
if no dist and sample name files provided, use whole genome data. 
Tip group files is used to plot points on branch ends, if none, no points plotted.
                       ",
                       option_list = option_list)

arguments <- parse_args(parser, positional_arguments=c(0,Inf))

distF <- arguments$options$distF
sampleF <- arguments$options$sampleF
annoF <- arguments$options$annoF
sampleP <- arguments$options$sampleP
colorF <- arguments$options$colorF
tmpG <- arguments$options$tmpG
plotT <- arguments$options$plotT
titleP <- arguments$options$titleP
out <- arguments$options$out


distF <- '/data/user/yangzz/mapping/fieldergenomecompare/4.subgenome_tree/AABBDD_setname_maf01_miss2.mdist'
sampleF <- '/data/user/yangzz/mapping/fieldergenomecompare/4.subgenome_tree/meta_201_samplelist'
annoF <- '/data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_196_nongke_10+.txt'
sampleP <- '/data/user/yangzz/mapping/fieldergenomecompare/4.subgenome_tree/meta_201_samplelist'
colorF <- '/data/user/yangzz/mapping/fieldergenomecompare/coloryzz_6cols.txt'
tmpG <- NULL
plotT <- 'T'
titleP <-'AABB_setname_maf003_miss05'
out <- '/data/user/yangzz/mapping/fieldergenomecompare/4.subgenome_tree/AABB_setname_maf003_miss05'
# 
# 
# colorF <- '/data2/rawdata2/tree/combineXN_190803/colorF_6cols.txt'
# sampleP <- '/data2/rawdata2/tree/combineXN_190803/D_sub_190809/03.plot/sample_order_hex.txt'
# annoF <- '/data2/rawdata2/sample_metadata/withXN/metadata_all_v3_5groups.txt'
# sampleF <- '/data2/rawdata2/tree/combineXN_190803/D_sub_190809/03.plot/sample_order_hex.txt'
# distF <- '/data2/rawdata2/tree/combineXN_190803/tryruns/AABB_setname_maf003_miss05.mdist'


# read in data
dist_mat <- read.table(distF)

SM <- read.table(sampleF, header = F, stringsAsFactors = F)
colnames(dist_mat) = SM[[1]]
rownames(dist_mat) = SM[[1]]

meta_data <- read.table(annoF, header = F, stringsAsFactors = F, sep="\t")
colnames(meta_data) <- c('Sample','Name','label','Region')
meta_data <- meta_data[complete.cases(meta_data),]  # remove lines contain NA

present_sample <- read.table(sampleP, header = F, stringsAsFactors = F)

# pre-process
if(!is.null(tmpG)) {
  tmp_group <- read.table(tmpG, header = F, stringsAsFactors = F)
  meta_data[match(tmp_group[,1], meta_data$Sample),4] <- tmp_group[,2]
}

tmp_index <- rownames(dist_mat) %in% meta_data$Name[na.omit(match(present_sample[[1]], meta_data$Name))] # match return the position of x in table
dist_mat <- dist_mat[tmp_index, tmp_index]

sample_present <- meta_data[meta_data$Name %in% rownames(dist_mat),]

LIST <- sample_present$label
names(LIST) <- sample_present$Name #assign name to LIST
rownames(dist_mat) = LIST[rownames(dist_mat)]
colnames(dist_mat)=rownames(dist_mat)
dist_mat1 <- 1- dist_mat
dis_matrix <- as.dist(dist_mat1) #dist try to caculate diatance; as.dist try to coerce a object to a distance matrix


pheatmap::pheatmap(dist_mat1,display_numbers = F ,clustering_method = "ward.D2",
                   annotation_col = col_anno,
                   annotation_row = col_anno,
                   fontsize_row = 5,fontsize_col = 5,annotation_colors = col_lis)

hclust(dis_matrix,)

# define colors
colorFile <- read.table(colorF, header = T, stringsAsFactors = F, sep = "\t", comment.char = "")
colors <- colorFile[,2]
names(colors) <- colorFile[,1]

# NJ-tree

NJ_tree <- njs(dis_matrix)
# NJ_tree <- root(NJ_tree, "Spelt_Baulander_S115")

sample_present_tree <- sample_present[,c(3,1,2,4)]
DF <- full_join(as_tibble(NJ_tree), sample_present_tree, by = 'label')

if(!plotT){
  tree <- as.treedata(DF)
  write.beast(tree, paste0(out, ".nexus"))

  branchcol <- DF %>% filter(!is.na(label)) %>% mutate(recode(Group, !!!colors))  %>% mutate(type = "branch", line = "normal") %>% select(c(4,9,8,10))

  write("TREE_COLORS\nSEPARATOR TAB\nDATA", paste0(out, ".branchcol.txt"))
  write.table(branchcol, paste0(out, ".branchcol.txt"), append = T, col.names = F, row.names = F, quote = F, sep = "\t")
} else {
  col.vector <- vector(mode="character",length=nrow(NJ_tree$edge))
  n.tips <- length(NJ_tree$tip.label)
  col.vector[NJ_tree$edge[,2]>n.tips] <- "black"
  edge.data <- as.data.frame(NJ_tree$edge)
  for(i in seq_along(NJ_tree$tip.label)){
    edge.row <- as.numeric(rownames(edge.data[edge.data$V2==i,]))
    col.vector[edge.row] <- colors[DF$Region[match(NJ_tree$tip.label[i], DF$label)]]
  }
  tip.color <- colors[DF$Region]
  names(tip.color) <- as.character(DF$label)
  pdf( paste0(out, ".pdf"), width = 10, height = 10)
  plot(NJ_tree, type = "phylogram", show.tip.label = T, edge.color = col.vector, cex=0.2)  #plot.phylo (ape)
  legend("bottom", fill = colors, legend = names(colors), border = F, box.col = "white", horiz=T, cex = 1, x.intersp = 0)
  if(!is.null(titleP)) title(titleP, cex.main = 2)
  dev.off()
}
if(!count_nearst){
  mat1 <- as.matrix(dis_matrix)
  mat1[upper.tri(mat1)] <- t(mat1)[upper.tri(mat1)]
  mat1 <- as.data.frame(mat1)
  dt_sample <- data.frame(V1=seq(1,11), stringsAsFactors=FALSE)
  LIST <- sample_present$Sample
  names(LIST) <- sample_present$label #assign name to LIST
  colnames(mat1) = LIST[colnames(mat1)]
  
  for (i in colnames(mat1)){
    dt_sample[[i]]<-colnames(mat1)[order(mat1[[i]],decreasing=F)[1:11]]
  }
  dt_sample<- dt_sample[-1]
  t1 <- t(data.frame(dt_sample))
  write.table(t1, "/data/user/yangzz/mapping/fieldergenomecompare/4.subgenome_tree/nearest10.txt", append = F, col.names = F, row.names = F, quote = F, sep = "\t")
}

ann_colors = list(
  Time = c("white", "firebrick"),
  CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
  GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)

annotation_col<-data.frame(variety=factor(meta_data$Variety))
rownames(annotation_col) <- meta_data$label
colorFile <- t(colorFile)
pheatmap(dist_mat1,annotation_col=annotation_col,annotation_colors = col_lis)

