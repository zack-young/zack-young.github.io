#!/usr/bin/env Rscript

# generate file for iTOL

# libraries
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-d", "--distF"), dest = "distF", default = "AABB_all.mdist",
              help = "dist matrix file, n rows by n cols. [default: whole_all.mdist]"),
  make_option(c("-s", "--sampleF"), dest = "sampleF", default = "sample_info.txt",
              help = "sample name that matches dist file. [default: id_group.txt]"),
  make_option(c("-a", "--annoF"), dest = "annoF", default = "../metadata2.txt",
              help = "annotation file, which contains sample name to show, sample group. [default: metadata_all.txt]"),
  make_option(c("-p", "--sampleP"), dest = "sampleP", default = "AABBDD_all_sample.txt",
              help = "only sample name in this file will be plot [default: whole_all.txt]"),
  make_option(c("-c", "--colorF"), dest = "colorF", default = "colorFile_itol.txt",
              help = "group name and color [colorFile_itol.txt]"),
  make_option(c("-o", "--out"), dest = "out", default = "output",
              help = "output prefix [default: output]")
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
out <- arguments$options$out

# distance matrix and its info
TAB <- read.table(distF)
SM <- read.table(sampleF, stringsAsFactors = F)
colnames(TAB) = SM[[1]]
rownames(TAB) = SM[[1]]

# metadata about samples
meta_data <- read.table(annoF,
                        sep = "\t",
                        comment.char = "#",
                        header = T, stringsAsFactors = F)
SP <- read.table(sampleP)

tmp_index <- rownames(TAB) %in% meta_data$Name[match(SP[[1]], meta_data$Sample)]
TAB <- TAB[tmp_index, tmp_index]

sample_present <- meta_data[meta_data$Name %in% rownames(TAB),]
LIST <- sample_present$label
names(LIST) <- sample_present$Name

rownames(TAB) = LIST[rownames(TAB)]
dis_matrix <- as.dist(TAB)

# define colors

colorFile <- read.table(colorF, header = F, stringsAsFactors = F, sep = "\t", comment.char = "")
colors <- colorFile[,2]
names(colors) <- colorFile[,1]
# colors = c("#5fb951", "#E74D4A", "#D9B460", "#AB6193","#71ACDF", "#F194BE", "#424E9F")
# names(colors) = c(unique(sample_present$Group))
# colors <- colors[!is.na(names(colors))]
# colors <- c("#b7a47c", "#5fb951", "#b00d28", "#0076b9")
# names(colors) <- c("NTC", "TS", "TL", "NC")

# NJ-tree

NJ_tree <- nj(dis_matrix)
# NJ_tree <- root(NJ_tree, "Spelt_Baulander_S115")

sample_present_tree <- sample_present[,c(3,1,2,4)]
DF <- full_join(as_data_frame(NJ_tree), sample_present_tree, by = 'label')

tree <- full_join(as_tibble(NJ_tree), sample_present_tree, by = 'label') %>% as.treedata
write.beast(tree, paste0(out, ".nexus"))

df1 <- DF %>% filter(!is.na(label)) %>% mutate(recode(Group, !!!colors))  %>% mutate(type = "branch", line = "normal")

write("TREE_COLORS\nSEPARATOR TAB\nDATA", paste0(out, ".branchcol.txt"))
write.table(df1[,c(4,9,8,10)], paste0(out, ".branchcol.txt"), append = T, col.names = F, row.names = F, quote = F, sep = "\t")

