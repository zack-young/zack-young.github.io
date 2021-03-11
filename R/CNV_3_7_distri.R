karyotype_use <-"/data2/rawdata2/readDepth/DP_bySampel/CS_karyotype.txt"
karyotype <- read.table(karyotype_use, header = T, stringsAsFactors = F)
toplot <- c("A", "B", "D")
tmplist <- list()
for (i in 1:nrow(karyotype)){
  chrID <- karyotype[i,1]
  #DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/crossover_distri/",chrID,".combine_hmm_crossover_count_compensent",sep=""),as.is = T, header = F, comment.char = "")

  DF1 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/CNV_count/",chrID,".combine_mask_CNV_deletion_compensent",sep=""),as.is = T, header = F, comment.char = "")
  DF2 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/CNV_count/",chrID,".combine_mask_CNV_duplication_compensent",sep=""),as.is = T, header = F, comment.char = "")
  DF1$V1 <- DF1$V1/1000000
  DF1$V2 <- DF1$V2/1000000
  DF2$V1 <- DF2$V1/1000000
  DF2$V2 <- DF2$V2/1000000
  DF <- cbind(DF1,DF2,rep(0,nrow(DF1)))
  colnames(DF) <- c("start1","end1","del","start2","end2","dup","zero")
  DF$del <- DF$del*(-1)
  
  p <- ggplot(DF)+ 
    geom_ribbon(aes(start1,ymin=del,ymax=0), fill="#4889E6", col="#4889E6") +
    geom_ribbon(aes(start1,ymin=0,ymax=dup), fill="red", col="red")+
    geom_line(aes(start1,del), col="#4889E6",size=1)+
    geom_line(aes(start1,dup), col="red",size=1)+
    geom_line(aes(start1,zero), col="black",size=1)+
    theme_bw()+
    theme(panel.grid =element_blank(),panel.border = element_blank(),
          axis.line.y = element_line(size=0.5, colour = "black"),
          axis.ticks.y=element_line(size=0.5, colour = "black"),
          axis.text.y=element_text(size=15,color="black"),
          axis.title.y=element_text(size=16,color="black"),
          axis.ticks.x=element_line(size=0.5, colour = "black"),
          axis.title.x=element_text(size=16,color="black"), #element_text(size=16,color="black"),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.text.x=element_text(size=15,color="black"))+ #top,right,bottom,left    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
    
    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
    #geom_line(aes(density, mean), size=1.2, col="red") +
    #geom_vline(xintercept = 30) +
    #scale_x_continuous(limits=c(0,820000000),breaks = NULL )+ #breaks=seq(0,800000000,400000000),labels = c("0(Mb)",400,800)
    scale_y_continuous(limits=c(-1,1), breaks=c(-1,0,1),labels = c("100%",0,"100%")) +
    #scale_y_continuous(limits=c(0,0.25), breaks=c(0,0.5),labels = c(0,0.5)) +
    ylab("") + xlab(chrID)

  # ggsave(paste0(CHR, ".logscale.pdf"), p, height = 2.5, width = 4)
  tmplist[[i]] <- p
}

pall <- cowplot::plot_grid(plotlist = tmplist,ncol = 3)
pall
