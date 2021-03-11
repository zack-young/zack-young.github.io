par(mar=c(2, 2, 2, 2))
toplot <- c("A", "B", "D")
tmplist <- list()
list_5 <- data.frame()
karyotype_use <-"/data2/rawdata2/readDepth/DP_bySampel/CS_karyotype.txt"
karyotype <- read.table(karyotype_use, header = T, stringsAsFactors = F)
chr_max_length <- max(karyotype$Length)/1000000

for (i in 1:nrow(karyotype)){
  chrID <- karyotype[i,1]
  DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/crossover_distri/",chrID,".combine_hmm_crossover_count_compensent",sep=""),as.is = T, header = F, comment.char = "")
  colnames(DF) <- c("start1","end1","ratio")
  DF[,'chr'] <- chrID
  #DF$del <- DF$del*(-1)
  list_5 <- rbind(list_5,DF[which(DF$ratio>0.05),])
  p <- ggplot(DF)+ 
    geom_ribbon(aes(start1,ymin=0,ymax=ratio), fill="#4889E6", col="#4889E6") +
    #geom_ribbon(aes(start1,ymin=0,ymax=dup), fill="red", col="red")+
    geom_line(aes(start1,ratio), col="#4889E6",size=1)+
    #geom_line(aes(start1,dup), col="red",size=1)+
    #geom_line(aes(start1,zero), col="black",size=1)+
    cowplot::theme_cowplot()+
    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
    #geom_line(aes(density, mean), size=1.2, col="red") +
    #geom_vline(xintercept = 30) +
    scale_x_continuous(limits=c(0,835000000),breaks=seq(0,800000000,400000000),labels = c("0(Mb)",400,800))+ #breaks=seq(0,800000000,400000000),labels = c("0(Mb)",400,800)
    scale_y_continuous(limits=c(0,0.2), breaks=c(0,0.1,0.2),labels = c(0,0.1,0.2)) +
    ylab("") + xlab("")
  
  # ggsave(paste0(CHR, ".logscale.pdf"), p, height = 2.5, width = 4)
  tmplist[[i]] <- p
}

pall <- cowplot::plot_grid(plotlist = tmplist,ncol = 3)
