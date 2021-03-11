library(ggplot2)
karyotype_use <-"/data2/rawdata2/readDepth/DP_bySampel/CS_karyotype.txt"
karyotype_use <-"/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/CS_karyotype.txt"
karyotype <- read.table(karyotype_use, header = T, stringsAsFactors = F)
toplot <- c("A", "B", "D")

qrc_lis <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/QRC_list.txt",as.is = T, header = F, comment.char = "",sep = '\t')
qrc_lis$V2 <- as.numeric(qrc_lis$V2)
qrc_lis$V3 <- as.numeric(qrc_lis$V3)
#gene_lis <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/gene.txt",as.is = T, header = T, comment.char = "",sep = '\t')
#gene_lis$Position <- as.numeric(gene_lis$Position)
#gene_lis <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/important_genes.txt",as.is = T, header = T, comment.char = "",sep = '\t')
#colnames(gene_lis) <- c("Chr","Position","Position2","Gene2","Gene")
gene_lis <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/gene_all_list_merged.txt",as.is = T, header = T, comment.char = "",sep = '\t')
colnames(gene_lis) <- c("Chr","Position","Position2","Gene")

gene_lis$Position <- as.numeric(gene_lis$Position)/1000000

chr_num=1
total_num=1
tmplist <- list()
totallist <- list()
for (i in c('chr1A','chr2A','chr3A','chr4A','chr5A','chr6A','chr7A','chr1B','chr2B','chr3B','chr4B','chr5B','chr6B','chr7B','chr1D','chr2D','chr3D','chr4D','chr5D','chr6D','chr7D')){#seq(1,(2*nrow(karyotype)),2)){

  chrID <-i #chrID <- karyotype[i,1]#karyotype[ceiling(i/2),1]
  #DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/crossover_distri/",chrID,".combine_hmm_crossover_count_compensent",sep=""),as.is = T, header = F, comment.char = "")
  qrc_df <- qrc_lis[qrc_lis$V1==chrID,]
  if (nrow(qrc_df)==0) {
    qrc_df <- data.frame(V1=chrID, V2=0,V3=0,V4='')
  }
  gene_df <- gene_lis[gene_lis$Chr==chrID,]
  if (nrow(gene_df)==0) {
    gene_df <- data.frame(V1=chrID, Position=-100,Gene='')
  }
  # 
  # # ggsave(paste0(CHR, ".logscale.pdf"), p, height = 2.5, width = 4)

    typ = 'CNC'
    files = read_delim(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/",chrID,"_",typ,"_hap_count.txt",sep=""), "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
    files = files[which(files$X4 != 'duplication'&files$X4 != 'deletion'),]
    wide_file <- dcast(files,X1+X2~X4,value.var = 'X3')
    col_len <- ncol(wide_file)
    for (i in seq(1,nrow(wide_file))) {
      if (length(which(as.numeric(wide_file[i,3:(col_len-1)])>5)) > 0){
        wide_file[i,'5%'] = sum(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])>5)+2])
      }else{
        wide_file[i,'5%'] = 0
      }
      if (length(which(as.numeric(wide_file[i,3:(col_len-1)])>2)) > 0){
        wide_file[i,'1%'] = sum(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])>2)+2])          
      }else{
        wide_file[i,'1%'] = 0
      }
      
      wide_file[i,'all']=sum(wide_file[i,3:col_len],na.rm = T)
      
    }
    #aql_all <- melt(files,measure.vars = c("all_1","all_1_5%","all_5%"),id.vars = c( "V2"))
    p1 <- ggplot(wide_file) +
      geom_area(aes(x = X2/1000000, y = all,alpha=0.5),fill='#747070')+
      geom_area(aes(x =X2/1000000, y = `1%`),fill='#84D9DC')+
      geom_area(aes(x = X2/1000000, y = `5%`),fill='#F2675E')+
      geom_text(data = gene_df,aes(Position,60,label=Gene),angle=45,hjust = 0,size=3)+
      geom_segment(data = gene_df, aes(x = Position, xend = Position, 
                   y = rep(0,length(Position)), yend = rep(60,length(Position))
                 ), linetype=2)+
      
      #scale_fill_manual(values = col_lis)+     
      theme_bw()+
      theme(panel.grid = element_blank(), panel.border  = element_blank(),          
            axis.line.y = element_line(size=0.5, colour = "black"),
            axis.ticks.y=element_line(size=0.5, colour = "black"),
            axis.text.y=element_text(size=15,color="black"),
            axis.text.x=element_blank(),
            axis.title.y=element_text(size=16,color="black"),
            axis.ticks.x=element_line(size=0.5, colour = "black"),
            axis.title.x=element_blank(),
            axis.line.x = element_line(size=0.5, colour = "black"),
            legend.position  = "none",plot.margin = unit(c(0,1,1,0), "cm"))+
      #annotate(x=-5, xend=-5, y=0, yend=60, colour="black", lwd=1, geom="segment") +
      # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
      #geom_line(aes(density, mean), size=1.2, col="red") +
      #geom_vline(xintercept = 30) +
      scale_x_continuous(expand = c(0.02,0))+
      scale_y_continuous(expand = c(0,0),limits = c(0,62),breaks = c(0,20,40,60))+ylab("Counts") + xlab(chrID)
    
##    tmplist[[chr_num]] <- p1
##    chr_num <- chr_num +1
    typ = 'CNL'
    files = read_delim(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/",chrID,"_",typ,"_hap_count.txt",sep=""), "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
    files = files[which(files$X4 != 'duplication'&files$X4 != 'deletion'),]
    wide_file <- dcast(files,X1+X2~X4,value.var = 'X3')
    col_len <- ncol(wide_file)
    for (i in seq(1,nrow(wide_file))) {
      if (length(which(as.numeric(wide_file[i,3:(col_len-1)])>5)) > 0){
        wide_file[i,'5%'] = sum(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])>5)+2])
      }else{
        wide_file[i,'5%'] = 0
      }
      if (length(which(as.numeric(wide_file[i,3:(col_len-1)])>2)) > 0){
        wide_file[i,'1%'] = sum(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])>2)+2])          
      }else{
        wide_file[i,'1%'] = 0
      }
      
      wide_file[i,'all']=sum(wide_file[i,3:col_len],na.rm = T)
      
    }
    #aql_all <- melt(files,measure.vars = c("all_1","all_1_5%","all_5%"),id.vars = c( "V2"))
    p2 <- ggplot(wide_file) +
      geom_area(aes(x = X2/1000000, y = all,alpha=0.5),fill='#747070')+
      geom_area(aes(x =X2/1000000, y = `1%`),fill='#84D9DC')+
      geom_area(aes(x = X2/1000000, y = `5%`),fill='#F2675E')+
      geom_vline(xintercept = gene_df$Position,linetype=2)+
      
      #scale_fill_manual(values = col_lis)+     
      theme_bw()+
      theme(panel.grid = element_blank(), panel.border  = element_blank(),          
            axis.line.y = element_line(size=0.5, colour = "black"),
            axis.ticks.y=element_line(size=0.5, colour = "black"),
            axis.text.y=element_text(size=15,color="black"),
            axis.text.x=element_blank(),
            axis.title.y=element_text(size=16,color="black"),
            axis.ticks.x=element_line(size=0.5, colour = "black"),
            axis.title.x=element_blank(),
            axis.line.x = element_line(size=0.5, colour = "black"),
            legend.position  = "none",plot.margin = unit(c(0,1,1,0), "cm"))+
      #
      # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
      #geom_line(aes(density, mean), size=1.2, col="red") +
      #geom_vline(xintercept = 30) +
      scale_x_continuous(expand = c(0.02,0))+
      scale_y_continuous(expand = c(0,0),limits = c(0,60),breaks = c(0,20,40,60))+ylab("Counts") + xlab(chrID)#+
    
    
##    tmplist[[chr_num]] <- p2
##    chr_num <- chr_num +1
    
  DF1 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/entropy/",chrID,"_CNC_entropy",sep=""),as.is = T, header = F, comment.char = "")
  colnames(DF1) <- c("chr1","start1","num1")
  DF2 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/entropy/",chrID,"_CNL_entropy",sep=""),as.is = T, header = F, comment.char = "")
  colnames(DF2) <- c("chr2","start1","num2")
  DF3 <-merge(DF1,DF2,by = 'start1',all = T)
  #DF[,4][which(DF[,4] < 4)]= 4
  #DF[,4][which(DF[,4] > 7)] = 7
  #print(paste("max",max(DF$num),"min",min(DF$num)))
  #da <- smooth.spline(DF$start1,DF$num)
  p3 <- ggplot(DF3)+

    geom_line(aes(start1,num1), col="#619CFF",size=0.7)+
    geom_line(aes(start1,num2), col="#D85D4A",size=0.7)+
    geom_line(aes(start1,num1-num2), col="#6A3D9A",size=0.7)+
    geom_text(data = gene_df,aes(Position*1000000,5,label=Gene),angle=45,hjust = 0,size=3)+
    #geom_text(data = gene_df,aes(Position*1000000,6.8,label=Gene),angle=45,hjust = 0,size=3)+
    geom_vline(xintercept = gene_df$Position*1000000,linetype=2)+
    #annotate("rect", xmin = qrc_df$V2*1000000, xmax = qrc_df$V3*1000000, ymin = -Inf, ymax = 4, alpha = 0.5,
    #         fill = "#22BBEE")+
    theme_bw()+
    theme(panel.grid =element_blank(),panel.border = element_blank(),
          axis.line.y = element_line(size=0.5, colour = "black"),
          axis.ticks.y=element_line(size=0.5, colour = "black"),
          axis.text.y=element_text(size=15,color="black"),
          axis.title.y=element_text(size=16,color="black"),
          axis.ticks.x=element_line(size=0.5, colour = "black"),
          axis.title.x=element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.text.x=element_blank(),
          plot.margin = unit(c(1,1,1,0), "cm"))+ #top,right,bottom,left
    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
    #geom_line(aes(density, mean), size=1.2, col="red") +
    #geom_vline(xintercept = 30) +
    scale_x_continuous(expand = c(0.02,0))+ #limits=c(0,860000000), breaks = seq(0,860000000,200000000),
    scale_y_continuous(expand = c(0,0))+ 
    #scale_y_continuous(limits=c(4,7), breaks=c(4,5,6,7),labels = c(4,5,6,7)) +
    ylab("Entropy") + xlab("")  
  # chr_num=chr_num+1
    tmplist[[chr_num]] <- p3
    chr_num=chr_num+1
  
  DF1 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/pi_stat/",chrID,".CNC.1M_nohet.windowed.pi",sep=""),as.is = T, header = T, comment.char = "")
  colnames(DF1) <- c("chr1","start1","end1","num1","pi1")
  DF2 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/pi_stat/",chrID,".CNL.1M_nohet.windowed.pi",sep=""),as.is = T, header = T, comment.char = "")
  colnames(DF2) <- c("chr2","start1","end2","num2","pi2")
  DF3 <-merge(DF1,DF2,by = 'start1',all = T)
  #DF[,5] <- -log(DF[,5],10)
  #DF[,5][which(DF[,5] > 0.005)]= 0.005
  #DF[,5][which(DF[,5] < 2.5)] = 2.5
  p4 <- ggplot(DF3)+
    #geom_ribbon(aes(start1,ymin=0,ymax=num), fill="red", col="red")+
    geom_line(aes(start1,pi1), col="#619CFF",size=0.7)+
    geom_line(aes(start1,pi2), col="#D85D4A",size=0.7)+
    geom_line(aes(start1,pi1-pi2), col="#6A3D9A",size=0.7)+
    #geom_line(aes(start1,zero), col="black",size=1)+
    #cowplot::theme_cowplot()+
    #geom_text(data = qrc_df,aes(V2*1000000,2,label=V4),angle=45,hjust = 0,size=3)+
    geom_vline(xintercept = gene_df$Position*1000000,linetype=2)+
    #geom_vline(xintercept = qrc_df$V2*1000000)+
    #annotate("rect", xmin = qrc_df$V2*1000000, xmax = qrc_df$V3*1000000, ymin = -Inf, ymax = 4, alpha = 0.5,
    #         fill = "#22BBEE")+
    theme_bw()+
    theme(panel.grid =element_blank(),panel.border = element_blank(),
           axis.line.y = element_line(size=0.5, colour = "black"),
           axis.ticks.y=element_line(size=0.5, colour = "black"),
           axis.text.y=element_text(size=15,color="black"),
           axis.title.y=element_text(size=16,color="black"),
           axis.ticks.x=element_line(size=0.5, colour = "black"),
           axis.title.x=element_text(size=15,color="black"),
           axis.line.x = element_line(size=0.5, colour = "black"),
           axis.text.x=element_text(size=15,color="black"),
           plot.margin = unit(c(0,1,0,0), "cm"))+ #top,right,bottom,left
    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
    #geom_line(aes(density, mean), size=1.2, col="red") +
    #geom_vline(xintercept = 30) +
    scale_x_continuous(expand = c(0.02,0))+ #limits=c(0,860000000), breaks = seq(0,860000000,200000000),
    scale_y_continuous(expand = c(0,0))+#coord_fixed(ratio=2)+
    # scale_y_continuous(limits=c(0,0.005), breaks=c(0,0.002,0.005),labels = c(0,0.002,0.005)) +
    ylab(expression(pi)) + xlab(chrID)
  #
  # # ggsave(paste0(CHR, ".logscale.pdf"), p, height = 2.5, width = 4)
  # chr_num=chr_num+1
  tmplist[[chr_num]] <- p4
  chr_num=chr_num+1
#  DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/pi_stat/",chrID,".CNC.1M_nohet.Tajima.D",sep=""),as.is = T, header = T, comment.char = "")
#  colnames(DF) <- c("chr","start1","end1","num")
  #DF[,4][which(DF[,4] < -3)]= -3
  #DF[,4][which(DF[,4] > 3)] = 3
#  p <- ggplot(DF)+
    #geom_ribbon(aes(start1,ymin=0,ymax=num), fill="red", col="red")+
#    geom_line(aes(start1+1,num), col="#4889E6",size=0.5)+
#    geom_vline(xintercept = gene_df$Position*1000000,linetype=2,size=0.5)+
    #geom_line(aes(start1,dup), col="red",size=1)+
    #geom_line(aes(start1,zero), col="black",size=1)+
    #cowplot::theme_cowplot()+
    #geom_text(data = qrc_df,aes(V2*1000000,2,label=V4),angle=45,hjust = 0,size=3)+
    #geom_vline(xintercept = qrc_df$V2*1000000)+
    #annotate("rect", xmin = qrc_df$V2*1000000, xmax = qrc_df$V3*1000000, ymin = -Inf, ymax = 4, alpha = 0.5,
    #         fill = "#22BBEE")+
#   theme_bw()+
#   theme(panel.grid =element_blank(),panel.border = element_blank(),
#         axis.line.y = element_line(size=0.5, colour = "black"),
#         axis.ticks.y=element_line(size=0.5, colour = "black"),
#         axis.text.y=element_text(size=15,color="black"),
#         axis.title.y=element_text(size=16,color="black"),
#         axis.ticks.x=element_line(size=0.5, colour = "black"),
#         axis.title.x=element_text(size=15,color="black"),
#         axis.line.x = element_line(size=0.5, colour = "black"),
#         axis.text.x=element_blank(),#element_text(size=15,color="black"),
#         plot.margin = unit(c(0,1,1,0), "cm"))+ #top,right,bottom,left    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
    #geom_line(aes(density, mean), size=1.2, col="red") +
    #geom_vline(xintercept = 30) +
#    scale_x_continuous(expand = c(0.02,0))+ #limits=c(0,860000000), breaks = seq(0,860000000,200000000),
#    scale_y_continuous(expand = c(0,0))+
   #scale_y_continuous(limits=c(-3,3), breaks=c(-3,0,3),labels = c(-3,0,3)) +
#    ylab("Tajima`s D") + xlab(chrID)

  # ggsave(paste0(CHR, ".logscale.pdf"), p, height = 2.5, width = 4)
#  tmplist[[chr_num]] <- p
#  chr_num=chr_num+1
  pall <- cowplot::plot_grid(plotlist = tmplist,ncol = 1,byrow=F,align = 'v')
  chr_num=1
  tmplist=list()
  totallist[[total_num]] <- pall
  total_num =total_num+ 1
  # 
}

pall <- cowplot::plot_grid(plotlist = totallist,ncol = 3,byrow=F,align = 'v',scale = 0.9)
pdf("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/pdf/H_pi_delta.pdf", height = 35, width = 25)

pall
dev.off()
da2 <- spline(DF$start1,DF$num,n = 100)
plot(DF$start1,DF$num)
lines(da2$x,da2$y,col='red')
lines(da$x,da$y,col='red')


