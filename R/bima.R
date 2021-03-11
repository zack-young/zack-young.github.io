toplot <- c("A", "B", "D")
options(scipen=200)
karyo_lis <- c("chr1A", "chr2A" ,"chr3A" ,"chr4A" ,"chr5A" ,"chr6A", "chr7A" ,"chr1B", "chr2B" ,"chr3B" ,"chr4B" ,"chr5B" ,"chr6B" ,"chr7B" ,"chr1D" ,"chr2D" ,"chr3D" ,"chr4D" ,"chr5D" ,"chr6D" ,"chr7D")
qrc_lis <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/QRC_list_filter.txt",as.is = T, header = F, comment.char = "",sep = '\t')
tmplist <- list()
#annotate draw line
for (i in 1:length(karyo_lis)){
  chrID <- karyo_lis[i]
  #DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/crossover_distri/",chrID,".combine_hmm_crossover_count_compensent",sep=""),as.is = T, header = F, comment.char = "")
  qrc_df <- qrc_lis[qrc_lis$V1==chrID,]
  DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare/S213_PH46/",chrID,".1M.delcnv_density",sep=""),as.is = T, header = T, comment.char = "")
  colnames(DF) <- c("chr","start1","end1","ratio")
  DF$ratio <- log(as.numeric(DF[["ratio"]])+10,10)
  END_pos <- max(DF$end1)
  #DF$del <- DF$del*(-1)
  if (nrow(qrc_df)!= 0){
    p <- ggplot(DF)+ 
      #geom_ribbon(aes(start1,ymin=del,ymax=0), fill="#4889E6", col="#4889E6") +
      #geom_ribbon(aes(start1,ymin=0,ymax=ratio), fill="red", col="red")+
      geom_line(aes(start1,ratio), col="black",size=0.5)+
      geom_segment(aes(y= 1.3010,yend=1.3010,x=-Inf,xend=835000000),col="red",size=0.5)+
      cowplot::theme_cowplot()+
      # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
      #geom_line(aes(density, mean), size=1.2, col="red") +
      geom_text(data = qrc_df,aes(V2*1000000,4,label=V4),angle=45,hjust = 0)+
      annotate("rect", xmin = qrc_df$V2*1000000, xmax = qrc_df$V3*1000000, ymin = -Inf, ymax = 4, alpha = 0.5, 
               fill = "#22BBEE")+
      #geom_vline(xintercept = qrc_df$V3*1000000) +
      scale_x_continuous(limits=c(0,835000000),breaks=c(0,400000000,800000000),labels = c("0(Mb)",400,800))+#,breaks = NULL )+ #breaks=seq(0,800000000,400000000),labels = c("0(Mb)",400,800)
      scale_y_continuous(limits=c(1,4.5), breaks=c(1,1.30103,2.70757,4.000434),labels = c(0,10,500,10000)) +
      ylab("") + xlab(chrID)
  }else{
    p <- ggplot(DF)+ 
      geom_line(aes(start1,ratio), col="black",size=0.5)+
      cowplot::theme_cowplot()+
      geom_segment(aes(y= 1.3010,yend=1.3010,x=-Inf,xend=835000000),col="red",size=0.5)+
      # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
      #geom_line(aes(density, mean), size=1.2, col="red") +
      #geom_vline(xintercept = qrc_df$V3*1000000) +
      scale_x_continuous(limits=c(0,835000000),breaks=c(0,400000000,800000000),labels = c("0(Mb)",400,800))+#,breaks = NULL )+ #breaks=seq(0,800000000,400000000),labels = c("0(Mb)",400,800)
      scale_y_continuous(limits=c(1,4.5), breaks=c(1,1.30103,2.70757,4.000434),labels = c(0,10,500,10000))+
      ylab("") + xlab(chrID)
  }
  # ggsave(paste0(CHR, ".logscale.pdf"), p, height = 2.5, width = 4)
    tmplist[[i]] <- p

} 

pall <-cowplot::plot_grid(plotlist = tmplist,ncol = 7)
ggsave(pall, filename = "/data/user/yangzz/mapping/fieldergenomecompare/pdf/bima_14_qtl_ann.pdf", height = 8, width = 35, limitsize = FALSE)
#line
cowplot::plot_grid(plotlist = tmplist,ncol = 7)
pall <- cowplot::plot_grid(plotlist = tmplist,ncol = 3)
#######annotate draw point
tmplist <- list()
karyotype_use <-"/data2/rawdata2/readDepth/DP_bySampel/CS_karyotype.txt"
karyotype <- read.table(karyotype_use, header = T, stringsAsFactors = F)
#for (i in c("A","B","D")){
for (i in 1:nrow(karyotype)){
  chrID <- karyotype[i,1]
  #qrc_df <- qrc_lis[grep(i, qrc_lis$V1, value = F),]
  #qrc_df$V2 <- qrc_df$V2*1000000
  #qrc_df$V3 <- qrc_df$V3*1000000
  #DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare/S213_PH46/",i,".1M.density_level",sep=""),as.is = T, header = F, comment.char = "")
  #S213_PH46
  DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare/S213_CS/",chrID,".1M.density_level",sep=""),as.is = T, header = T, comment.char = "")
  colnames(DF)<-c("V1","V2","V3","V4","V5")
  DF$V2 <- DF$V2/1000000
  DF$V3 <- DF$V3/1000000
  #for (num in seq(1,nrow(qrc_df))){
  #  qrc_df[num,5] <- as.numeric(rownames(DF[which(DF$V1==qrc_df[num,1]&DF$V2==qrc_df[num,2]+1),]))*1000000
  #  qrc_df[num,6] <- as.numeric(rownames(DF[which(DF$V1==qrc_df[num,1]&DF$V3==qrc_df[num,3]+1),]))*1000000
  #} 
  #DF[1,2] <-1
  #DF[1,3] <-1000000
  #DF[2:nrow(DF),2] <-seq(1,nrow(DF)-1,by = 1)*1000000
  #DF[2:nrow(DF),3] <-seq(1,nrow(DF)-1,by = 1)*1000000+1000000
  #DF$V1 <- factor(DF$V1, 
  #                levels = c(paste("chr",seq(1,7),i,sep = "")),
  #                labels = c(rep(c(1,2),3),1)
  #)
  DF$V5[stringr::str_detect(string =DF$V5,pattern = 'CNV' )] <- "CNV"
  #diff_df <- DF[DF$V5=='diff',]
  CNV_df <- DF[DF$V5=='CNV',]
  DF$V5 <- factor(DF$V5, 
                  levels = c("low","mid","high","CNV"),
                  labels = c(1,2,3,4))
  DF$V4 <-  log(DF$V4+10,10)
  p <- ggplot(DF)+
    annotate("rect", xmin = CNV_df$V2, xmax = CNV_df$V3, ymin =1, ymax = 4, alpha = 1, 
             fill = "gray")+
    geom_point(aes(V2,V4,colour=V5),size=0.5)+scale_colour_manual(values=c("#22BAEE","#334899","#58AA61","gray"),breaks = c(1,2,3,4))+
    cowplot::theme_cowplot()+
    #geom_segment(aes(y= 1.3010,yend=1.3010,x=-Inf,xend=835000000),col="red",size=0.5)+
    
    #annotate("rect", xmin = qrc_df$V5, xmax = qrc_df$V6, ymin = 4, ymax = 4.2, alpha = 1, 
    #         fill = "#ED1C24")+
    #geom_text(data = qrc_df,aes(V5,4.2,label=V4),angle=20,hjust = 0,)+
    
    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
    #geom_line(aes(density, mean), size=1.2, col="red") +
    #geom_vline(xintercept = qrc_df$V3*1000000) +
    #scale_x_continuous(limits=c(0,835000000),breaks=c(0,400000000,800000000),labels = c("0(Mb)",400,800))+#,breaks = NULL )+ #breaks=seq(0,800000000,400000000),labels = c("0(Mb)",400,800)
    scale_y_continuous(limits=c(1,4.5), breaks=c(1,1.30103,3,4.000434),labels = c(0,10,1000,10000))+
    ylab("DSr") + xlab(chrID)+theme(legend.position="none")#+theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  p
  #ggsave(p, filename = paste("/data/user/yangzz/mapping/fieldergenomecompare/pdf/",i,"_bima_14_qtl_ann.pdf",sep=""), height = 1.3, width = 15, limitsize = FALSE)
  
  tmplist[[i]] <- p
}
pall <-cowplot::plot_grid(plotlist = tmplist,nrow= 7)
#ggsave(pall, filename = "/data/user/yangzz/mapping/fieldergenomecompare/pdf/bima_14_qtl_ann.pdf", height = 8, width = 15, limitsize = FALSE)


