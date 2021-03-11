suppressPackageStartupMessages(library(sf, lib.loc = "/data/user/shinyug/R/x86_64-redhat-linux-gnu-library/3.5/"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-g", "--geoDB"), dest = "geoDB", default = "/data2/rawdata2/hapmap/script/sample_location_intro.txt",
              help = "geographic data file, accession location lon lat [default: sample_location.txt]"),
  make_option(c("-i", "--sampleF"), dest = "sampleF", default = "/data2/rawdata2/tetraintro/200621/analyze/intro_diffusion/group_3.txt",
              help = "sample group file, accession group [default: id_group.txt]"),
  make_option(c("-r", "--GeoRange"), dest = "GeoRange", default = "-20,150,0,60",
              help = "longitude and latitude range, sep by comma [-70,140,-20,71]"),
  make_option(c("-t", "--titleP"), dest = "titleP", default = "geographic distribution",
              help = "output prefix [default: output]"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "group5.pdf",
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

geoDB <- arguments$options$geoDB
sampleF <- arguments$options$sampleF
GeoRange <- arguments$options$GeoRange
titleP <- arguments$options$titleP
outfile <- arguments$options$outfile

# longitude and latitude range, four number
GeoRange <- as.numeric(strsplit(GeoRange, ",")[[1]])

# sample geographical location
geoDB <- read.table(geoDB, header = F, stringsAsFactors = F)

# genotype data of one site eg. chr1A.1:1000 (a snp site)
sample <- read.table(sampleF, header = F, stringsAsFactors = F)

b <- merge(geoDB,sample,by = 'V1')
names(b) <- c("sample", "location", "Longitude", "Latitude", "Ratio")

# world <- ne_countries(scale = "medium", returnclass = "sf")
# save(file = "world", world)
load("world")

ggplot(data = world) +
  geom_sf(col="white") +
  geom_jitter(data = b, aes(x = Longitude, y = Latitude, size = Ratio, col = Ratio)) +
  scale_color_viridis_c(limits=c(0.01,1), oob = scales::squish, trans = "log2",breaks=c(0.001,0.01,0.1,1)) +
  scale_size(range = c(0, 5)) +
  coord_sf(xlim=c(GeoRange[1],GeoRange[2]),ylim=c(GeoRange[3],GeoRange[4])) + 
  labs(title = titleP) +
  guides(size = FALSE) +
  cowplot::theme_cowplot() +
  
ggsave(outfile, width = 7, height = 3.5)
