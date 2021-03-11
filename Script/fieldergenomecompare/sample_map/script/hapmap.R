suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggmap))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-g", "--geoDB"), dest = "geoDB", default = "/data2/rawdata2/hapmap/script/sample_location_WGSwithXN_WEC1026.txt",
              help = "geographic data file, accession location lon lat [default: sample_location.txt]"),
  make_option(c("-i", "--sampleF"), dest = "sampleF", default = "/data2/rawdata2/altitude/30.5_altitude_adaptive_gene/DD1B-3B_group/clade_map_input.txt",
              help = "sample group file, accession group [default: id_group.txt]"),
  make_option(c("-c", "--colorF"), dest = "colorF", default = "/data2/rawdata2/altitude/pub/colorFile_clade.txt",
              help = "color file, group color [colorFile_itol.txt]"),
  make_option(c("-r", "--GeoRange"), dest = "GeoRange", default = "-70,140,-20,71",
              help = "longitude and latitude range, sep by comma [-70,140,-20,71]"),
  make_option(c("-T", "--capP"), dest = "capP", default = NA,
              help = "plot caption [default: output]"),
  make_option(c("-m", "--mergeD"), dest = "mergeD", default = 5,
              help = "merge samples within distance [5]"),
  make_option(c("-s", "--pointS"), dest = "pointS", default = 0.3,
              help = "point size factor [0.3]"),
  make_option(c("-t", "--titleP"), dest = "titleP", default = "geographic distribution",
              help = "output prefix [default: output]"),
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

geoDB <- arguments$options$geoDB
sampleF <- arguments$options$sampleF
colorF <- arguments$options$colorF
GeoRange <- arguments$options$GeoRange
mergeD <- as.numeric(arguments$options$mergeD)
pointS <- as.numeric(arguments$options$pointS)
titleP <- arguments$options$titleP
capP <- arguments$options$capP
out <- arguments$options$out

# longitude and latitude range, four number
GeoRange <- as.numeric(strsplit(GeoRange, ",")[[1]])

# merge locs between which dist less than disthr
mergeloc <- function(locdf, disthr = 0.1){

  # merge same location, ggmap use location as group
  if(disthr < 0.1){disthr = 0.1}
  
  newlocdf <- locdf[FALSE,]
  for (i in seq_len(nrow(locdf))){
    if (i == 1){
      newlocdf <- rbind(newlocdf, locdf[i,])
    }else{
      locnum <- nrow(newlocdf)
      n = 0
      for (j in seq_len(nrow(newlocdf))){
        if (as.numeric(dist(rbind(locdf[i,c(3,4)], newlocdf[j,c(3,4)]))) < disthr){
          newlocdf[j,-c(1:4)] <- newlocdf[j,-c(1:4)] + locdf[i,-c(1:4)]
          break
        }else{
          n=n+1
        }
      }
      if (n == locnum){
        newlocdf <- rbind(newlocdf, locdf[i,])
      }
    }
  }
  return(newlocdf)
}

# sample geographical location
geoDB <- read.table(geoDB, header = TRUE, stringsAsFactors = F)

# genotype data of one site eg. chr1A.1:1000 (a snp site)
sample <- read.table(sampleF, header = TRUE, stringsAsFactors = F)
total_attri <- colnames(sample)[-1]

dff <- inner_join(geoDB, sample, by = 'Accession')
dff <- dff[complete.cases(dff),]

# disthr should change with size of map to avoid points overlaping
dff <- mergeloc(dff, disthr = mergeD)
dff$radius <- log(apply(dff[,-c(1:4)], 1, function(x) sum(x)) +1)*pointS
pie <- dff[,-c(1,2)]

# define colors
if(!is.na(colorF)){
  colorFile <- read.table(colorF, header = T, stringsAsFactors = F, sep = "\t", comment.char = "")
  colors <- colorFile[,2]
  names(colors) <- colorFile[,1]
}

# ggplot map data
world = map_data("world", resolution=0)

filenametag <- NULL
parameter <- NULL

p <- ggplot(data = world, aes(x=long, y=lat, group=group)) + 
  geom_polygon(fill = "#DFE0E1", color = "black",size=0.25) + 
  coord_quickmap(xlim = c(GeoRange[1], GeoRange[2]), ylim = c(GeoRange[3], GeoRange[4])) +
  labs(title = titleP) +
  ylab("Latitude") + 
  xlab("Longitude") + 
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(colour = "grey90", size = 0.5), 
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.spacing.x = unit(0.01, "npc"),
    legend.justification = "right",
    plot.title = element_text(size = 25,hjust=0.5),
    axis.title = element_text(size=20),
    strip.text = element_text(size=15),
    plot.margin = grid::unit(c(0.05,0.05,0.05,0.05), "npc"))

# it seems grob by ggplotGrob cant clip to plot
# grob by *Grob func works fine.
if(!is.na(colorF)){
  piegrob <- function(d){
    d2 <- as.data.frame(d)
    if (max(d2$value)==sum(d2$value)){
      d <- d2[which(d2$value==max(d2$value)),]
    }
    ggplotGrob(ggplot(d, 
                      aes(x = 1, y = value, fill = type)) +
                 geom_col( size=0.0000001,#color = "black",
                          show.legend = FALSE) +
                 coord_polar(theta = "y") +
                 scale_fill_manual(values = colors) +
                 theme_void())
  }
} else {
  piegrob <- function(d){
    d2 <- as.data.frame(d)
    if (max(d2$value)==sum(d2$value)){
      d <- d2[which(d2$value==max(d2$value)),]
    }
    ggplotGrob(ggplot(d, 
                      aes(x = 1, y = value, fill = type)) +
                 geom_col( size=0.000001,#color = "black",
                          show.legend = FALSE) +
                 coord_polar(theta = "y") +
                 theme_void())
  }
}

pie.list <- pie %>% 
  filter(lon > GeoRange[1] & lon < GeoRange[2] & lat > GeoRange[3] & lat < GeoRange[4]) %>%
  tidyr::gather(type, value, -lon, -lat, -radius) %>%
  tidyr::nest(type, value) %>%
  
  # make a pie chart from each row, & convert to grob
  mutate(pie.grob = purrr::map(data, piegrob)) %>%
  
  # convert each grob to an annotation_custom layer. I've also adjusted the radius
  # value to a reasonable size (based on my screen resolutions).
  rowwise() %>%
  mutate(radius = radius * 4) %>%
  mutate(subgrob = list(annotation_custom(grob = pie.grob,
                                          xmin = lon - radius, xmax = lon + radius,
                                          ymin = lat - radius, ymax = lat + radius)))

p <- p + 
  
  # Optional. this hides some tiles of the corresponding color scale BEHIND the
  # pie charts, in order to create a legend for them
  geom_tile(data = pie %>% tidyr::gather(type, value, -lon, -lat, -radius),
            aes(x = lon, y = lat, fill = type), 
            color = "black", width = 0.01, height = 0.01, 
            inherit.aes = FALSE) +
  pie.list$subgrob

if(!is.na(colorF)) p <- p + scale_fill_manual(values = colors)
if(!is.na(capP)) p <- p + labs(caption = capP) + theme(plot.caption = element_text(size = 15))

ggsave(paste0(out,".pdf"), p, width = 12, height = 7)
