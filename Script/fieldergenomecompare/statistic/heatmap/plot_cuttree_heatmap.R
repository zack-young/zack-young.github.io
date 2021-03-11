#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggdendro))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(cowplot))
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
gname <- arguments$options$gname
outfile <- arguments$options$outfile
region <- arguments$options$region
posi <- arguments$args


posi<-'not'

# check arguments
if(infile == ""){ # default, STDIN
  infile <- file("stdin")
} else { # user specified
  if( file.access(infile) == -1){ # file not exists
    print_help(parser)
  }
}

if(outfile == ""){ # default, "$gene-name.pdf"
  outfile <- gname
}


# input file
#infile <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmpfile" 
tmp <- read.table(infile, header = T, check.names = F, sep = "\t", stringsAsFactors = F)
pos_data <- tmp$ANN
cn <- colnames(tmp)[-1]  #-1 delete first column
rn <- as.character(tmp$ANN)
tmp <- as.data.frame(t(tmp[,-1]))  #t() is to transposition a matrix
colnames(tmp) <- rn
tmp$Sample <- cn

# group file
#groupf <- read.table("/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmpfileG", col.names = c("Sample", "Group"), sep = "\t", stringsAsFactors = F)

groupf <- read.table(arguments$options$groupf, col.names = c("Sample", "Group"), sep = "\t", stringsAsFactors = F)

# metadata
#meta_data <- read.table("/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/metadata_all.append.txt", header = T, sep = "\t", stringsAsFactors = F)

meta_data <- read.table(arguments$options$anno, header = T, sep = "\t", stringsAsFactors = F)

# add info
sample_present <- meta_data[meta_data$Name %in% tmp$Sample,]
LIST <- sample_present$label
names(LIST) <- sample_present$Name
tmp$Group <- groupf$Group[match(tmp$Sample, groupf$Sample)]
tmp$Sample = LIST[tmp$Sample]

if(arguments$options$revxy){
  tmp$Sample <- factor(tmp$Sample, levels = rev(tmp$Sample))
} else{
  tmp$Sample <- factor(tmp$Sample, levels = tmp$Sample)
}


# re-order  #tryCatch catch error or warning and use function to solve
# tryCatch(
#   error = function(cnd) {
#     tmp$Sample <- factor(tmp$Sample)
#   },
#   {
data_set <- subset(tmp, select = -c(Sample,Group))
row.names(data_set) <- tmp$Sample
suppressWarnings(tmp_clust <- hclust(dist(data_set)))   # why tmp[,-1]
row.order <- tmp_clust$order
tmp$Sample <- factor(tmp$Sample, levels = tmp$Sample[row.order])
dend <- as.dendrogram(tmp_clust)
# Cut the dendrogram
depth_cutoff <- 15
h_c_cut <- cut(dend, h = depth_cutoff) # this cut() is in dendextend to cut dendrogram
dend_cut <- as.dendrogram(h_c_cut$upper)
dend_cut <- hang.dendrogram(dend_cut)
# Format to extend the branches (optional)
dend_cut <- hang.dendrogram(dend_cut, hang = -1) 
dend_data_cut <- dendro_data(dend_cut)

# Extract the names assigned to the clusters (e.g., "Branch 1", "Branch 2", ...)
cluster_names <- as.character(dend_data_cut$labels$label)
# Extract the names of the haplotypes that belong to each group (using
# the 'labels' function)
lst_genes_in_clusters <- h_c_cut$lower %>% 
  lapply(labels) %>% 
  setNames(cluster_names)

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
  segment(dend_data_cut), 
  data.frame(x = y, y = x, xend = yend, yend = xend))
# Extract the positions of the clusters (by getting the positions of the 
# leafs); data is already in the same order as in the cluster name
cluster_positions <- segment_data[segment_data$xend == 0, "y"]
cluster_pos_table <- data.frame(y_position = cluster_positions, 
                                cluster = cluster_names)

# Specify the positions for the genes, accounting for the clusters
sample_pos_table <- lst_genes_in_clusters %>%
  ldply(function(ss) data.frame(Sample = ss), .id = "cluster") %>%
  mutate(y_center = 1:nrow(.), 
         height = 1)
write.csv(subset(sample_pos_table, select = c(cluster,Sample)), file=paste0("/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmp/",outfile,'.csv'))
# > head(gene_pos_table, 3)
#    cluster gene y_center height
# 1 Branch 1  g11        1      1
# 2 Branch 1  g20        2      1
# 3 Branch 1  g12        3      1
# Table to position the samples
base_pos_table <- data.frame(Mutation = pos_data) %>%
  mutate(x_center = 1:nrow(.), width = 1)
# Coordinates for plotting rectangles delimiting the clusters: aggregate
# over the positions of the genes in each cluster
cluster_delim_table <- sample_pos_table %>%
  group_by(cluster) %>%
  summarize(y_min = min(y_center - 0.5 * height), 
            y_max = max(y_center + 0.5 * height)) %>%
  as.data.frame() %>%
  mutate(x_min = with(base_pos_table, min(x_center - 0.5 * width)), 
         x_max = with(base_pos_table, max(x_center + 0.5 * width)))

# Limits for the vertical axes (genes / clusters)
gene_axis_limits <- with(
  sample_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1
# melt
heatmap_data <-tmp %>% 
  reshape2::melt( id.vars= c("Sample", "Group")) 
colnames(heatmap_data) <- c("Sample", "Group", "Mutation", "Type")

df <- left_join(heatmap_data,sample_pos_table) %>%
  left_join(base_pos_table)

df$Group <- factor(df$Group, levels = unique(df$Group))

# turn discrete
df$Type <- as.factor(df$Type)
levels(df$Type)[levels(df$Type)=="-1"] <- "Missing"
levels(df$Type)[levels(df$Type)=="0"] <- "None"
levels(df$Type)[levels(df$Type)=="0.5"] <- "Heter"
levels(df$Type)[levels(df$Type)=="1"] <- "Homo"

#
# Heatmap plot
plt_hmap <- ggplot(df, 
                   aes(x = x_center, y = y_center, fill = Type, 
                       height = height, width = width)) + 
  geom_tile() +
  geom_rect(data = cluster_delim_table, 
            aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), 
            fill = NA, colour = "black", inherit.aes = FALSE) + 
  scale_fill_manual(values = c("gray95", "grey", "deepskyblue", "blue")) +
  scale_x_continuous(breaks = base_pos_table$x_center, 
                     labels = base_pos_table$Mutation, 
                     expand = c(0.01, 0.01)) + 
  scale_y_continuous(breaks = sample_pos_table$y_center, 
                     labels = sample_pos_table$Sample, 
                     limits = gene_axis_limits, 
                     expand = c(0, 0), 
                     position = "right") + 
  labs(x = "Sample", y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(1), hjust = 1, angle = 45), 
        # margin: top, right, bottom, and left
        plot.margin = unit(c(1, 0.2, 0.2, -0.1), "cm"), 
        panel.grid.minor = element_blank())

# Dendrogram plot
plt_dendr <- ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse(expand = c(0, 0.5)) + 
  scale_y_continuous(breaks = cluster_pos_table$y_position, 
                     labels = cluster_pos_table$cluster, 
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Distance", y = "", colour = "", size = "") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank())
#
p <-plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(1, 1))
#outfile <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/test22.pdf"
ggsave(p, filename = paste0("/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmp/",outfile,'.pdf'), height =nrow(tmp) * 0.2 + 14, width =nrow(tmp) * 0.2 + 7 , limitsize = FALSE)
#  }
#)
#
