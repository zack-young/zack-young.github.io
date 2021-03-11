#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))

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
#gname <- "gname"
gname <- arguments$options$gname
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
  outfile <-  paste0(outfile , ".pdf", sep="")
  #gsub(".pdf$|.png$", "", outfile, perl=T)
}

# input file
#infile <- "/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmp.GFfNDBqmLe" 
tmp <- read.table(infile, header = T, check.names = F, sep = "\t", stringsAsFactors = F)

cn <- colnames(tmp)[-1]
rn <- as.character(tmp$ANN)
tmp <- as.data.frame(t(tmp[,-1]))
colnames(tmp) <- rn
tmp$Sample <- cn

# group file
#groupf <- read.table("/data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/tmp.vk3FelgMAc", col.names = c("Sample", "Group"), sep = "\t", stringsAsFactors = F)

groupf <- read.table(arguments$options$groupf, col.names = c("Sample", "Group"), sep = "\t", stringsAsFactors = F)

# metadata
#meta_data <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_final_headed.txt", header = T, sep = "\t", stringsAsFactors = F)

meta_data <- read.table(arguments$options$anno, header = F, sep = "\t", stringsAsFactors = F)
colnames(meta_data) <- c("Sample"	,"Name"	,"label"	,"Group")

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
    suppressWarnings(row.order <- hclust(dist(tmp[,]))$order)
    tmp$Sample <- factor(tmp$Sample, levels = tmp$Sample[row.order])
  }
)
}

# melt
df <- reshape2::melt(tmp, id.vars = c("Sample", "Group"))
colnames(df) <- c("Sample", "Group", "Mutation", "Type")
df$Group <- factor(df$Group, levels = unique(df$Group))

# turn discrete
df$Type <- as.factor(df$Type)
levels(df$Type)[levels(df$Type)=="-1"] <- "Missing"
levels(df$Type)[levels(df$Type)=="0"] <- "None"
levels(df$Type)[levels(df$Type)=="0.5"] <- "Heter"
levels(df$Type)[levels(df$Type)=="1"] <- "Homo"

# plot
if(arguments$options$revxy){
p <- ggplot(df,aes(Mutation,Sample,fill=Type))+
  geom_tile() +
  scale_fill_manual(values = c("gray95", "grey", "deepskyblue", "blue")) +
  #scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
  scale_x_discrete(limits = rev(levels(df$Mutation))) +
  theme_minimal()+
  labs(title= paste(gname, region, gsub("\\|", "\n",posi), sep="\n")) +
  facet_grid(rows = vars(Group), scales = "free", space = "free") +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.justification = "right",
        plot.title = element_text(size = 17,hjust=0.95),
        axis.text.x = element_text(angle = 70, hjust = 1),
        axis.title = element_text(size=25),
        strip.text = element_text(size=25),
        aspect.ratio = 1) +
  NULL
ggsave(p, filename = outfile,height = nrow(tmp) * 0.2 + 7, width = ncol(tmp) * 0.2 + 5, limitsize = FALSE)
} else {
  p <- ggplot(df,aes(Sample,Mutation,fill=Type))+
    geom_tile() +
    scale_fill_manual(values = c("gray95", "grey", "deepskyblue", "blue")) +
    #scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
    scale_y_discrete(limits = levels(df$Mutation)) +
    theme_minimal()+
    labs(title= paste(gname, region, gsub("\\|", "\n",posi), sep="\n")) +
    facet_grid(cols = vars(Group), scales = "free", space = "free") +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.justification = "right",
          plot.title = element_text(size = 17,hjust=0.95),
          axis.text.x = element_text(angle = 70, hjust = 1),
          axis.title = element_text(size=25),
          strip.text = element_text(size=25),
          aspect.ratio = 1) +
    NULL
  ggsave(p, filename = outfile,height = ncol(tmp) * 0.2 + 5, width = nrow(tmp) * 0.2 + 7, limitsize = FALSE)
}
