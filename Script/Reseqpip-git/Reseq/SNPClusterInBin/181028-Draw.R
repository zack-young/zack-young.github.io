#!/usr/bin/env Rscript

is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

if(!is.installed("optparse")){
    warning("Detect package \"optparse\" is not installed in your R enviroment.")
    warning("Trying to install the \"optparse\" package.")
    warning("If failed, please try to install it mannually.")
    install.packages("optparse")
}


## libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dendextend))

# Arguments
option_list <- list(
    make_option(c("-i", "--infile"), dest = "infile", default = "",
                            help="input file"),
    #make_option(c("-r", "--methylation_right"), dest = "methylation.right", default = TRUE,
    #                        help="[opt] show average methylation level on the right panel."),
    make_option(c("-o", "--outfile"), dest = "outfile", default = "",
                            help = "[opt] output file name. [default: Clustmap.SysDate.pdf]"),

    make_option(c("-c", "--commentfile"), dest = "commentfile", default = "",
                            help="sample name commit file. Two columns in the file, name in the front, with comment latter.\nFile should be read-able by read.table() in R."),

    make_option(c("-m", "--mapfile"), dest = "mapfile", default = "",
                            help="sample name transfer file. Origon name in the second column, with target name in the front.\nFile should be read-able by read.table() in R."),

    make_option(c("-t", "--title"), dest = "title", default = "",
                            help = "[opt] title of whole file. [default: ]"),

    make_option(c("-W","--width"), dest = "figure.width", default = 20,
                            help = "[opt] width of figure (inch). [default: 20]"),

    make_option(c("-H","--height"), dest = "figure.height", default = 40,
                            help = "[opt] height of figure (inch). [default: 40]"),

    make_option(c("-w","--treewidth"), dest = "treewidth", default = 2,
                            help = "[opt] width of clust tree, while the width of heatmap is 5. [default: 2] \nDefault value means clust tree have 2/7 of the whole width. while heatmap have 5/7."),

    make_option(c("-f","--format"), dest = "figure.format", default = "pdf",
                            help = "[opt] format of output figure. Alternative: png. [default: pdf]"),

    make_option(c("-R","--resolution"), dest = "figure.resolution", default = 300,
                            help = "[opt] Resolution in ppi. Only available for png format. [default: 300]"),

    make_option(c("-r","--rowhight"), dest = "row_hight", default = 0.9,
                            help = "[opt] Hight of every row of heatmap. [default: 0.9]"),

    make_option(c("-e","--extra"), dest = "extra", default = 0,
                            help = "[opt] extra message for debug. 0 for shut. [default: 0]"),

    make_option(c("-s","--mainsample"), dest = "main_sample", default = "",
                            help = "[opt] main sample list to mark, like -S1,S2,S3-, without \"-\". [default: ]")

)

parser <- OptionParser(usage = "mapdrawer [options]",
                                             option_list = option_list, description = "\
Description: Clust data drawer\
Author : Wang Wenxi\
Last update: 2018-11\
")

## check arguments
arguments <- parse_args(parser)

infile <- arguments$infile

if(infile == ""){ # default, STDIN
    infile <- file("stdin")
} else { # user specified
    if( file.access(infile) == -1){ # file not exists
        print_help(parser)
    }
}

outfile <- arguments$outfile

if(outfile == ""){ # default, "FragRegView.Date"
    outfile <- paste("Clustmap", Sys.Date(), sep=".")
} else { # user specified
    outfile <- gsub(".pdf$|.png$", "", outfile, perl=T)
}

figure.width <- arguments$figure.width

figure.height <- arguments$figure.height

figure.format <- arguments$figure.format

if(! figure.format %in% c("pdf", "png")){ # format not support
    print_help(parser)
} else {
    outfile <- paste(outfile, figure.format, sep = ".")
}

figure.resolution <- arguments$figure.resolution

data <- read.table(file = infile, header = T, stringsAsFactors = F, check.names = FALSE)

data[,1] <- NULL
data[,1] <- NULL
data[,1] <- NULL

## figure device
if (figure.format == "png"){
    png(outfile, height = figure.height, width = figure.width, res = figure.resolution, units = "in")
} else if (figure.format == "pdf"){
    pdf(outfile, height = figure.height, width = figure.width)
}

#initlize graphic 
if(arguments$title == ""){
    if(arguments$extra != 0){
        par(mar=c(4,0,0,0))
    }else{
        par(mar=c(4, 0, 0, 0), oma=c(0,4,2,3))
    }
}else{
    if(arguments$extra != 0){
        par(mar=c(4,0,1,1))
    }else{
        par(mar=c(4, 0, 1, 0), oma=c(0,4,2,3))
    }
}
layout(matrix(1:2, nrow=1), widths=c(arguments$treewidth,5))

row.height <- arguments$row_hight

#clust
clust_data <- data.frame()
for(i in 1:ncol(data)){
    v <- c()
    for(j in 1:ncol(data)){
        sam1 <- data[,i]
        sam2 <- data[,j]
        count <- 0
        for(k in 1:length(sam1)){
            if(sam1[k] != sam2[k]){
                count <- count+1
            }
        }
        v <- c(v,count/length(sam1))
    }
    clust_data <-rbind(clust_data, v) 
}
names(clust_data) <- names(data)
hc <- hclust(d = as.dist(clust_data), method = "ward.D2")
den <- as.dendrogram(hc)

#Change ylim to change the hight of clust tree
if(arguments$extra != 0){
    plot(den, horiz=T, yaxt = "none", yaxs="i",ylim=c(0, ncol(data)*(row.height + 0.1)+(ncol(data) - 3)*0.1 - 0.5))
}else{
    plot(den, horiz=T, leaflab="none", yaxt = "none", dLeaf = 0, yaxs="i",ylim=c(0, ncol(data)*(row.height + 0.1)+(ncol(data) - 3)*0.1 - 0.5))
}

data <- data[, order.dendrogram(den)]


#sample name transfer
mapfile <- arguments$mapfile
sample_name_output <- names(data)
if(mapfile != ""){
    map_file <- read.table(mapfile, header = F, as.is = T)
    for(i in 1:length(sample_name_output)){
        for(j in 1:nrow(map_file)){
            if(sample_name_output[i] == map_file[j,2]){
                sample_name_output[i] <- map_file[j,1]
                break
            }
        }
    }
}

#commit sample name
commentfile <- arguments$commentfile
#sample_name_output <- names(data)
if(commentfile != ""){
    comment_file <- read.table(commentfile, header = F, as.is = T)
    for(i in 1:length(sample_name_output)){
        for(j in 1:nrow(comment_file)){
            if(sample_name_output[i] == comment_file[j,1]){
                sample_name_output[i] <- paste(sample_name_output[i], "(", comment_file[j,2], ")", sep="")
                break
            }
        }
    }
}

#heat_map
color_pad <- c("#000000", "#e5e5e5", "#00FFFF", "#1E90FF", "#66CD00", "#FFFF00", "#FFC1C1", "#CD5C5C", "#8b658b", "#EE00EE", "#FF0000", "#DAA520", "#FF8C00", "#FFA07A", "#FF4500", "#BC8F8F", "#D2691E", "#FFB6C1", "#DC143C", "#DB7093", "#C71585", "#8B008B", "#6A5ACD", "#191970", "#2F4F4F", "#006400", "#556B2F", "#800000", "#4B0082", "#4F4F4F", "#68228B", "#CDB38B", "#0000CD", "#708090")

if(arguments$title == ""){
    plot(x=0, type="n", bty="n", yaxt="n",# xaxt="n",
        xlab="", ylab="", 
        xlim=c(1, nrow(data) * 1.3), ylim=c(1, ncol(data)*(row.height + 0.1)+ncol(data)*0.1), 
        xaxs="i", yaxs="i")
}else{
    plot(x=0, type="n", bty="n", yaxt="n",# xaxt="n",
        xlab="", ylab="", 
        xlim=c(1, nrow(data) * 1.3), ylim=c(1, ncol(data)*(row.height + 0.1)+ncol(data)*0.1), 
        xaxs="i", yaxs="i", main=arguments$title)
}



rightPanel.extension <- 1
main_sample <- arguments$main_sample
main_sample <- strsplit(main_sample, split = ",")
main_sample <- main_sample[[1]]
for(i in 1:ncol(data)){
	y.i <- i
	color <- c()
	swi <- F
	for (j in 1:length(data[,i])){
		if((-as.numeric(data[j,i])+2) > length(color_pad)){
			color <- c(color, "#000000")
		}else{
			color <- c(color, color_pad[(-as.numeric(data[j,i])+2)])
		}
	}

	rect(xleft = 1:nrow(data) - 1, ybottom = rep(y.i*(row.height + 0.1), nrow(data)) + 0.5, xright = 1:nrow(data), ytop = rep(y.i*(row.height + 0.1), nrow(data)) + row.height-0.1 + 0.5, col = color, border = color)
    if(length(main_sample) != 0){
        for(j in 1:length(main_sample)){
            if(main_sample[j] == names(data)[i]){
                rect(xleft = mean(as.numeric(data[,i]), na.rm = T) * nrow(data)*(rightPanel.extension-1) + nrow(data) * 1.01 + 1,
                    ybottom = y.i*(row.height + 0.1) + 0.5,
                    xright = mean(as.numeric(data[,i]), na.rm = T) * nrow(data)*(rightPanel.extension-1) + nrow(data) * 1.01 + 1 + 15,
                    ytop = y.i*(row.height + 0.1) + row.height - 0.1 + 0.5,
                    col = color_pad[j+2],
                    border = color_pad[j+2]
                    )
                break
            }
        }
    }
    text(x=mean(as.numeric(data[,i]), na.rm = T) * nrow(data)*(rightPanel.extension-1) + nrow(data) * 1.01 + 1,
		y=y.i*(row.height + 0.1)+0.1+row.height/2 + 0.5,
		sample_name_output[i],#bar后面的文字
		cex = 0.7, adj = c(0, 0.5)
		)
}
#text(x=5,y=ncol(data)*(row.height + 0.1),string_otp,cex = 0.7, adj = c(0, 0.5))

invisible(dev.off())