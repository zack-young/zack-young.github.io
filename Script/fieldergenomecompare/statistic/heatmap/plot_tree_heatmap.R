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

#gname <- 'gname'
#region <- 'chr1A:1-10'
#posi<-'not'

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
infile <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmp.Jso2GusaGX" 
#infile <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmpfile" 

tmp <- read.table(infile, header = T, check.names = F, sep = "\t", stringsAsFactors = F)
pos_data <- tmp$ANN
cn <- colnames(tmp)[-1]  #-1 delete first column
rn <- as.character(tmp$ANN)
tmp <- as.data.frame(t(tmp[,-1]))  #t() is to transposition a matrix
colnames(tmp) <- rn
tmp$Sample <- cn

# group file
groupf <- read.table("/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmp.urHQBYDgAM", col.names = c("Sample", "Group"), sep = "\t", stringsAsFactors = F)
#groupf <- read.table("/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmpfileG", col.names = c("Sample", "Group"), sep = "\t", stringsAsFactors = F)

#groupf <- read.table(arguments$options$groupf, col.names = c("Sample", "Group"), sep = "\t", stringsAsFactors = F)

# metadata
meta_data <- read.table("/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/metadata_cultivar_noheader.txt", header = T, sep = "\t", stringsAsFactors = F)

#meta_data <- read.table(arguments$options$anno, header = T, sep = "\t", stringsAsFactors = F)

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

if(arguments$options$clus){
# re-order  #tryCatch catch error or warning and use function to solve
tryCatch(
  error = function(cnd) {
    tmp$Sample <- factor(tmp$Sample)
  },
  {
    data_set <- subset(tmp, select = -c(Sample,Group))
    row.names(data_set) <- tmp$Sample
    suppressWarnings(tmp_clust <- hclust(dist(data_set)))   # why tmp[,-1]
    row.order <- tmp_clust$order
    dend <- as.dendrogram(tmp_clust)
    tmp$Sample <- factor(tmp$Sample, levels = tmp$Sample[row.order])
    
    # Obtain the dendrogram
    dend_data <- dendro_data(dend)  #Extract Cluster Data From A Model Into A List Of Data Frames.
    # Setup the data, so that the layout is inverted (this is more 
    # "clear" than simply using coord_flip()) !
    segment_data <- with(             #with() makes segment to be a environment
      segment(dend_data), 
      data.frame(x = y, y = x, xend = yend, yend = xend))
    # Use the dendrogram label data to position the gene labels
    sample_pos_table <- with(
      dend_data$labels, 
      data.frame(y_center = x, Sample = as.character(label), height = 1))
    # Table to position the samples
    base_pos_table <- data.frame(Mutation = pos_data) %>%
      mutate(x_center = (1:n()), 
             width = 1)
    # Limits for the vertical axes
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
    plt_dendr <- ggplot(segment_data) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      scale_x_reverse(expand = c(0, 0.5)) + 
      scale_y_continuous(breaks = sample_pos_table$y_center, 
                         labels = sample_pos_table$Sample, 
                         limits = gene_axis_limits, 
                         expand = c(0, 0)) + 
      labs(x = "Distance", y = "", colour = "", size = "") +
      theme_bw() + 
      theme(panel.grid.minor = element_blank())
    
    plt_hmap <- ggplot(df, 
                       aes(x = x_center, y = y_center, fill = Type, 
                           height = height, width = width)) + 
      geom_tile() +
      scale_fill_manual(values = c("gray95", "grey", "deepskyblue", "blue")) +
      scale_x_continuous(breaks = base_pos_table$x_center, 
                         labels = base_pos_table$Mutation, 
                         expand = c(0, 0)) + 
      # For the y axis, alternatively set the labels as: gene_position_table$gene
      scale_y_continuous(breaks = sample_pos_table[, "y_center"], 
                         labels = rep("", nrow(sample_pos_table)),
                         limits = gene_axis_limits, 
                         expand = c(0, 0)) + 
      labs(x = "Sample", y = "") +
      theme_bw() +
      theme(axis.text.x = element_text(size = rel(1), hjust = 1, angle = 45), 
            # margin: top, right, bottom, and left
            plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
            panel.grid.minor = element_blank())
    #
    p <-plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(0.1, 1))
    ggsave(p, filename = paste0("/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmp/",outfile,'.pdf'), height =nrow(tmp) * 0.2 + 14, width =nrow(tmp) * 0.2 + 7 , limitsize = FALSE)
    
  }
)
}

