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

    make_option(c("-e","--extra"), dest = "extra", default = 0,
                            help = "[opt] extra message for debug. 0 for shut. [default: 0]"),

    make_option(c("-s","--mainsample"), dest = "main_sample", default = "",
                            help = "[opt] main sample list to mark, like -S1,S2,S3-, without \"-\". [default: ]")

)

parser <- OptionParser(usage = "GroupFregCount [options]",
                                             option_list = option_list, description = "\
Description: Clusted data same group fregments count\
Author : Wang Wenxi\
Last update: 2019-01\
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

data <- read.table(file = infile, header = T, stringsAsFactors = F, check.names = FALSE)

data[,1] <- NULL
data[,1] <- NULL
data[,1] <- NULL

main_sample <- arguments$main_sample
main_sample <- strsplit(main_sample, split = ",")
main_sample <- main_sample[[1]]

#sample name transfer
mapfile <- arguments$mapfile
sample_name_output <- main_sample
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


main_sample <- arguments$main_sample
main_sample <- strsplit(main_sample, split = ",")
main_sample <- main_sample[[1]]

res_fra <- data.frame(sample=sample_name_output, count=rep(0, length(sample_name_output)))
for(i in 1:ncol(data)){
    for(j in 1:nrow(data)){
        if(-as.numeric(data[j,i]) > 0){
            res_fra[-as.numeric(data[j,i]), 2] <- res_fra[-as.numeric(data[j,i]), 2] + 1
        }
    }
}

write.table(res_fra, outfile, quote = F, row.names = F)