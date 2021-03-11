#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(dplyr))
# Arguments
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="input file"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file name. [default: gene name]"),
  make_option(c("-a", "--anno"), dest = "anno", default = "metadata_all.append.txt",
              help = "[opt] annotation file. [default: metadata_all.txt]"),
  make_option(c("-g", "--groupf"), dest = "groupf", default = "",
              help = "[opt] group file, sample-name TAB group, one per line"),
  make_option(c("-n", "--gname"), dest = "gname", default = "",
              help="the name of gene to plot"),
  make_option(c("-r", "--region"), dest = "region", default = "",
              help = "example: chr1A:1-100"),
  make_option(c("-c", "--clus"), dest = "clus", default = T,
              help = "reorder samples in each groups according to their genotypes or not"),
  make_option(c("-R", "--revxy"), dest = "revxy", default = F,
              help = "filp x and y axis or not")
)

parser <- OptionParser(usage = "Rscript haplotypePlot.R [options]",
                       option_list=option_list)

arguments <- parse_args(parser, positional_arguments=c(0,Inf))

infile <- arguments$options$infile
gname <- "gname"
#gname <- arguments$options$gname
outfile <- arguments$options$outfile
region <- arguments$options$region
posi <- arguments$args
# check arguments
if(infile == ""){ # default, STDIN
  infile <- file("stdin")
} else { # user specified
  if( file.access(infile) == -1){ # file not exists
    print_help(parser)
  }
}

if(outfile == ""){ # default, "$gene-name.pdf"
  outfile <- paste0(gname, ".pdf")
} else { # user specified
  outfile <- gsub(".pdf$|.png$", "", outfile, perl=T)
}

# input file
#infile <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmp.GORxntdmCK" 
tmp <- read.table(infile, header = T, check.names = F, sep = "\t", stringsAsFactors = F)

cn <- colnames(tmp)[-1]
rn <- as.character(tmp$ANN)
tmp <- as.data.frame(t(tmp[,-1]))
colnames(tmp) <- rn
#tmp$Sample <- cn

# group file
# group_name <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmp.0auvrslX2J" 
# groupf <- read.table(group_name, col.names = c("Sample", "Group"), sep = "\t", stringsAsFactors = F)

groupf <- read.table(arguments$options$groupf, col.names = c("Sample", "Group"), sep = "\t", stringsAsFactors = F)

# metadata

# meta_name <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/metadata_cultivar_noheader.txt" 
# meta_data <- read.table(meta_name, header = T, sep = "\t", stringsAsFactors = F)

meta_data <- read.table(arguments$options$anno, header = T, sep = "\t", stringsAsFactors = F)

# add info
#sample_present <- meta_data[meta_data$Name %in% tmp$Sample,]
sample_present <- meta_data[meta_data$Name %in% rownames(tmp),]
LIST <- sample_present$label
names(LIST) <- sample_present$Name
#tmp$Group <- groupf$Group[match(tmp$Sample, groupf$Sample)]
#tmp$Sample = LIST[tmp$Sample]
rownames(tmp) <- LIST[rownames(tmp)]
# if(arguments$options$revxy){
#   tmp$Sample <- factor(tmp$Sample, levels = rev(tmp$Sample))
# } else{
#   tmp$Sample <- factor(tmp$Sample, levels = tmp$Sample)
# }

#if(arguments$options$clus){
  # re-order  #tryCatch catch error or warning and use function to solve
  # tryCatch(
  #   error = function(cnd) {
  #     tmp$Sample <- factor(tmp$Sample)
  #   },
  #   {
  #     suppressWarnings(out.clust <- hclust(dist(tmp[,-1])))
  #     #suppressWarnings(row.order <- hclust(dist(tmp[,-1]))$order)
  #     tmp$Sample <- factor(tmp$Sample, levels = tmp$Sample[out.clust$order])
  #   }
  # )
#}
#
out.clust <- hclust(dist(tmp))
NJ_tree <- njs(dist(tmp))
#out.clust$labels <- tmp$Sample
#dendrogram <- as.dendrogram(out.clust,hang = -1)
#pdf(paste(outfile,".pdf",sep=""), height = 6, width = nrow(tmp) * 0.1 + 7)
pdf(paste(outfile,".pdf",sep=""), width = 6, height = nrow(tmp) * 0.1 + 7)
#plot(dendrogram)#, hang = -1,labels=NULL)
#plot(out.clust, hang = -1,labels=NULL)



plot(NJ_tree,type='p',cex=0.5,use.edge.length = F)







dev.off()
