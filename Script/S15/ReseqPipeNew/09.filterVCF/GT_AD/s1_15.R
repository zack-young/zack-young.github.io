#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2018-08-08

library(ggplot2)
require(cowplot)
library(Rmisc)
library(lattice)
library(plyr)
files = read.table("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plotfile.pt", as.is = T, header = F, comment.char = "")
pdf("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/chr3A.pdf", height = 20, width = 11.69)
for(j in 1:15){  #nrow(files)
  #data1 = read.table(files[j,1], as.is = T, header = F, comment.char = "")
  sample <- paste('s',j,sep = "")
  datas = read.table(files[j,1], as.is = T, header = F, comment.char = "")
  ratio <- datas[,6]
  chr <- datas[,2]
  sample_name = paste("plot",j,sep="")
  print(sample_name)
  a <- ggplot(datas, aes(x=chr, y=ratio)) + geom_line(linetype="solid") + ylim(0, 0.5) + xlab(sample)
  assign(sample_name,a)#   + xlab(sample))
}
# plot.1 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.2 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.3 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.4 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.5 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.6 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.7 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.8 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.9 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.10 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.11 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.12 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.13 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.14 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")
# plot.15 <- ggplot(data = files, aes(y=ratio, x=region)) +
#   geom_line(linetype="solid")   # ctrl+shift+c plus#
#multiplot(plot4,plot3, layout = matrix(c(1,2), nrow=2, byrow=TRUE))
ggdraw() +
  draw_plot(plot1, 0, 0.1, 1, 0.06) + 
  draw_plot(plot2, 0, 0.15, 1, 0.08) +
  draw_plot(plot3, 0, 0.25, 1, 0.08) 
#  draw_plot(plot4, 0, 0.25, 1, 0.06) +
#  draw_plot(plot5, 0, 0.30, 1, 0.06) +
#  draw_plot(plot6, 0, 0.35, 1, 0.06) +
#  draw_plot(plot7, 0, 0.4, 1, 0.06) +
#  draw_plot(plot8, 0, 0.45, 1, 0.06) +
#  draw_plot(plot9, 0, 0.5, 1, 0.06) +
#  draw_plot(plot10, 0, 0.55, 1, 0.06) +
#  draw_plot(plot11, 0, 0.6, 1, 0.06) +
#  draw_plot(plot12, 0, 0.65, 1, 0.06) +
#  draw_plot(plot13, 0, 0.7, 1, 0.06) +
#  draw_plot(plot14, 0, 0.75, 1, 0.06) +
#  draw_plot(plot15, 0, 0.8, 1, 0.06)
  #draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)
dev.off()
