library(reshape2)
library(readr)

gene_lis <- read.table("/data2/rawdata2/tetraintro/bin/important_genes.txt",as.is = T, header = T, comment.char = "",sep = '\t')
colnames(gene_lis) <- c("Chr","Position","Position2","Gene2","Gene")
gene_lis$Position <- as.numeric(gene_lis$Position)/1000000
qrc_lis <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/QRC_list.txt",as.is = T, header = F, comment.char = "",sep = '\t')
qrc_lis$V2 <- as.numeric(qrc_lis$V2)
qrc_lis$V3 <- as.numeric(qrc_lis$V3)


gene_lis <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/QRC_list.txt",as.is = T, header = F, comment.char = "",sep = '\t')


col_lis <- c("#E31A1C" ,"#F7D000" ,"#F71F60" ,"#EB832E" ,
             "#14ACCC" ,"#2FD7BA" ,"#33A02C", "#B2DF8A", "#A6CEE3",
             "#F7568B", "#FDBF6F" ,"#FF7F00", "#FB9A99",
             "#FFFF99", "#6A3D9A","#858687")
names(col_lis)=c(seq(1,15),"exclusive")

qrc_lis <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/Ppd1_FT2_5A.txt",as.is = T, header = F, comment.char = "",sep = '\t')
colnames(qrc_lis) <- c("Chr","Position","Position2","Gene")  #var_cwi5D_bas1_gene.txt
qrc_lis$Position <- as.numeric(qrc_lis$Position)/1000000

tmplist <- list()
chr_num=1
for (i in unique(qrc_lis$Chr)){#seq(1,(2*nrow(karyotype)),2)){
  #DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/crossover_distri/",chrID,".combine_hmm_crossover_count_compensent",sep=""),as.is = T, header = F, comment.char = "")
  chrID <-i
  qrc_df <- qrc_lis[qrc_lis$Chr==chrID,]
#  if (nrow(qrc_df)==0) {
#    qrc_df <- data.frame(V1=chrID, V2=0,V3=0,V4='')
#  }
#  gene_df <- gene_lis[gene_lis$Chr==chrID,]
#  if (nrow(gene_df)==0) {
#    gene_df <- data.frame(V1=chrID, Position=-100,Gene='')
#  }
   #chrID <- karyotype[i,1]#karyotype[ceiling(i/2),1]
  for (a in seq(1,nrow(qrc_df)) ){
    
            # files = files[which(files$X4!='exclusive'&files$X4 != 'duplication'&files$X4 != 'deletion'),]
      # files <- files[which(files$X2>=(qrc_df[a,2]-25)*1000000&files$X2-1<=(qrc_df[a,2]+25)*1000000),]
      # wide_file <- dcast(files,X1+X2~X4,value.var = 'X3')
      # col_len <- ncol(wide_file)
      # for (i in seq(1,nrow(wide_file))) {
      #   if (length(which(as.numeric(wide_file[i,3:(col_len-1)])>5)) > 0){
      #     wide_file[i,'5%'] = sum(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])>5)+2])
      #   }else{
      #     wide_file[i,'5%'] = 0
      #   }
      #   if (length(which(as.numeric(wide_file[i,3:(col_len-1)])>2)) > 0){
      #     wide_file[i,'1%'] = sum(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])>2)+2])          
      #   }else{
      #     wide_file[i,'1%'] = 0
      #   }
      #   
      #   wide_file[i,'all']=sum(wide_file[i,3:col_len],na.rm = T)
      #   
      # }
      #aql_all <- melt(files,measure.vars = c("all_1","all_1_5%","all_5%"),id.vars = c( "V2"))
      files = read_delim(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/",chrID,"_CNC_hap_count.txt",sep=""), "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
      files <- files[which(files$X2>=(qrc_df[a,2]-150)*1000000&files$X2-1<=(qrc_df[a,2]+150)*1000000),]
      DF <- files[which(files$X4 != 'duplication'&files$X4 != 'deletion'),]
      uniq_lis <- unique(sort(as.numeric(DF[which(DF$X4!='exclusive'),]$X4),decreasing = T))
      DF$X4<- factor(DF$X4, levels = c('exclusive',uniq_lis)) #,labels = seq(length(unique(sort(DF$V4))),1))
      
      p <- ggplot(DF) +
        geom_bar( aes(x = X2/1000000, y = X3,fill=X4),stat = "identity") +
        scale_fill_manual(values = col_lis)+   
        geom_text(data=qrc_df[a,],aes(Position,55,label=Gene),hjust = 0,size=6)+
        geom_vline(xintercept = qrc_df[a,]$Position,linetype=2,size=0.3)+
        
        theme_bw()+
        theme(panel.grid = element_blank(), panel.border = element_blank(),          
              axis.line.y = element_line(size=0.5, colour = "black"),
              axis.ticks.y=element_line(size=0.5, colour = "black"),
              axis.text.y=element_text(size=15,color="black"),
              axis.text.x=element_text(size=15,color="black"),
              axis.title.y=element_text(size=16,color="black"),
              axis.ticks.x=element_line(size=0.5, colour = "black"),
              axis.title.x=element_text(size=15,color="black"),
              axis.line.x = element_line(size=0.5, colour = "black"),#legend.position  = "none",
              plot.margin = unit(c(1,1,0,0), "cm"))+
        # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
        #geom_line(aes(density, mean), size=1.2, col="red") +
        #geom_vline(xintercept = 30) +
        scale_x_continuous(expand = c(0.02,0))+
        scale_y_continuous(expand = c(0,0),limits = c(0,60))+ylab("") + xlab(chrID)#,limits = c(0,6),breaks = c(0,2,4,6))+
      
      #tmplist[[chr_num]] <- p
      #chr_num <- chr_num +1
      
      DF1 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/entropy/",chrID,"_CNC_entropy",sep=""),as.is = T, header = F, comment.char = "")
      colnames(DF1) <- c("chr1","start1","num1")
      DF2 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/entropy/",chrID,"_CNL_entropy",sep=""),as.is = T, header = F, comment.char = "")
      colnames(DF2) <- c("chr2","start1","num2")
      DF3 <-merge(DF1,DF2,by = 'start1',all = T)

      DF3 <- DF3[which(DF3$start1>=(qrc_df[a,2]-150)*1000000&DF3$start1-1<=(qrc_df[a,2]+150)*1000000),]
      p3 <- ggplot(DF3)+
        
        geom_line(aes(start1,num1), col="#619CFF",size=0.7)+
        geom_line(aes(start1,num2), col="#D85D4A",size=0.7)+
        geom_line(aes(start1,num1-num2), col="#6A3D9A",size=0.7)+
        geom_text(data = qrc_df,aes(Position*1000000,5,label=Gene),angle=45,hjust = 0,size=3)+
        #geom_text(data = gene_df,aes(Position*1000000,6.8,label=Gene),angle=45,hjust = 0,size=3)+
        geom_vline(xintercept = qrc_df[a,]$Position*1000000,linetype=2)+
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
        #scale_x_continuous(limits=c(0,860000000), breaks = seq(0,860000000,200000000) )+
        #scale_y_continuous(limits=c(0,5), breaks=c(0,2.5,5),labels = c(0,2.5,5)) +
        #ylab(expression(H["CNC"]^"'")) 
      # chr_num=chr_num+1
      tmplist[[chr_num]] <- p3
      chr_num=chr_num+1
      DF1 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/pi_stat/",chrID,".CNC.1M_nohet.windowed.pi",sep=""),as.is = T, header = T, comment.char = "")
      colnames(DF1) <- c("chr1","start1","end1","num1","pi1")
      DF2 <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/pi_stat/",chrID,".CNL.1M_nohet.windowed.pi",sep=""),as.is = T, header = T, comment.char = "")
      colnames(DF2) <- c("chr2","start1","end2","num2","pi2")
      DF3 <-merge(DF1,DF2,by = 'start1',all = T)

      DF3 <- DF3[which(DF3$start1>=(qrc_df[a,2]-150)*1000000&DF3$start1-1<=(qrc_df[a,2]+150)*1000000),]
      p4 <- ggplot(DF3)+
        #geom_ribbon(aes(start1,ymin=0,ymax=num), fill="red", col="red")+
        geom_line(aes(start1,pi1), col="#619CFF",size=0.7)+
        geom_line(aes(start1,pi2), col="#D85D4A",size=0.7)+
        geom_line(aes(start1,pi1-pi2), col="#6A3D9A",size=0.7)+
        #geom_line(aes(start1,zero), col="black",size=1)+
        #cowplot::theme_cowplot()+
        #geom_text(data = qrc_df,aes(V2*1000000,2,label=V4),angle=45,hjust = 0,size=3)+
        geom_vline(xintercept = qrc_df$Position*1000000,linetype=2)+
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
      
      DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/pi_stat/",chrID,"_CNC.windowed.pi",sep=""),as.is = T, header = T, comment.char = "")
      colnames(DF) <- c("chr","start1","end1","num","pi1")
      DF <- DF[which(DF$start1>=(qrc_df[a,2]-150)*1000000&DF$start1-1<=(qrc_df[a,2]+150)*1000000),]
      #DF[,5] <- -log(DF[,5],10)
      DF[,5][which(DF[,5] > 0.003)]= 0.003
      p <- ggplot(DF)+
        #geom_ribbon(aes(start1,ymin=0,ymax=num), fill="red", col="red")+
        geom_line(aes(start1/1000000,pi1), col="#8FCA86",size=0.5)+
        #geom_line(aes(start1,dup), col="red",size=1)+
        #geom_line(aes(start1,zero), col="black",size=1)+
        #cowplot::theme_cowplot()+
        geom_text(data = qrc_df,aes(V2*1000000,2,label=V4),angle=45,hjust = 0,size=3)+
        geom_vline(xintercept = qrc_df[a,]$Position,linetype=2,size=0.3)+
        #annotate("rect", xmin = qrc_df$V2*1000000, xmax = qrc_df$V3*1000000, ymin = -Inf, ymax = 0.2, alpha = 0.5,
        #         fill = "#22BBEE")+
        theme_bw()+
        theme(panel.grid =element_blank(),panel.border = element_blank(),
              axis.line.y = element_line(size=0.5, colour = "black"),
              axis.ticks.y=element_line(size=0.5, colour = "black"),
              axis.text.y=element_text(size=15,color="black"),
              axis.title.y=element_text(size=16,color="black"),
              axis.ticks.x=element_line(size=0.5, colour = "black"),
              axis.title.x=element_text(size=16,color="black"), #element_text(size=16,color="black"),
              axis.line.x = element_line(size=0.5, colour = "black"),
              axis.text.x=element_blank(),#element_blank(),
              plot.margin = unit(c(0,1,1,0), "cm"))+ #top,right,bottom,left    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
        # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
        #geom_line(aes(density, mean), size=1.2, col="red") +
        #geom_vline(xintercept = 30) +
        scale_x_continuous(expand = c(0.02,0))+
        scale_y_continuous(expand = c(0,0),limits=c(-4,4))+ 
        #scale_x_continuous(limits=c(0,860000000), breaks = seq(0,860000000,200000000) )+
        #scale_y_continuous(limits=c(0,0.003), breaks=c(0,0.001,0.002,0.003),labels = c(0,0.001,0.002,0.003)) +
        ylab(expression(pi)) + xlab(chrID)
      #ylab(expression(H["all"]^"'"-H["CN"]^"'")) + xlab("")
      # # ggsave(paste0(CHR, ".logscale.pdf"), p, height = 2.5, width = 4)
      # chr_num=chr_num+1
      #tmplist[[chr_num]] <- p
      #chr_num=chr_num+1
      
      
      # DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/CNV_CNC_count/",chrID,".combine_mask_CNV_compensent",sep=""),as.is = T, header = T, comment.char = "")
      # #colnames(DF) <- c(CHR     START   Del     Dup     Com)
      # DF <- DF[which(DF$START>=(qrc_df[a,2]-25)*1000000&DF$START-1<=(qrc_df[a,2]+25)*1000000),]
      # DF[,4][which(DF[,4] > 0.7)] = 0.7
      # p <- ggplot(DF)+
      #   geom_line(aes(START,Del), col="#9373AD",size=0.5)+
      #   geom_line(aes(START,Dup), col="#EA7559",size=0.5)+
      #   geom_vline(xintercept = qrc_df$Position*1000000,linetype=2,size=0.3)+
      #   theme_bw()+
      #   theme(panel.grid =element_blank(),panel.border = element_blank(),
      #         axis.line.y = element_line(size=0.5, colour = "black"),
      #         axis.ticks.y=element_line(size=0.5, colour = "black"),
      #         axis.text.y=element_text(size=15,color="black"),
      #         axis.title.y=element_text(size=16,color="black"),
      #         axis.ticks.x=element_line(size=0.5, colour = "black"),
      #         axis.title.x=element_text(size=16,color="black"),#element_blank(),
      #         axis.line.x = element_line(size=0.5, colour = "black"),
      #         axis.text.x=element_blank(),#element_text(size=15,color="black"),
      #         plot.margin = unit(c(0,1,1,0), "cm"))+ #top,right,bottom,left    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
      #   ylab("CNV ratio") + xlab(chrID)
      # tmplist[[chr_num]] <- p
      # chr_num=chr_num+1
    }
}

pall <- cowplot::plot_grid(plotlist = tmplist,nrow = 6,byrow=F,align = 'v')
pdf("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/pdf/H_hap_CNC_gene.pdf", height = 10, width = 20)

pall
dev.off()

  # p <- ggplot(DF)+
  #   #geom_ribbon(aes(start1,ymin=0,ymax=num), fill="red", col="red")+
  #   geom_line(aes(start1,pi2), col="#8FCA86",size=0.5)+
  #   #geom_line(aes(start1,dup), col="red",size=1)+
  #   #geom_line(aes(start1,zero), col="black",size=1)+
  #   #cowplot::theme_cowplot()+
  #   #geom_text(data = qrc_df,aes(V2*1000000,2,label=V4),angle=45,hjust = 0,size=3)+
  #   #geom_text(data = gene_df,aes(Position*1000000,1,label=Gene),angle=30,hjust = 0,size=3)+
  #   geom_vline(xintercept = gene_df$Position*1000000,linetype=2)+
  #   #annotate("rect", xmin = qrc_df$V2*1000000, xmax = qrc_df$V3*1000000, ymin = -Inf, ymax = 4, alpha = 0.5,
  #   #         fill = "#22BBEE")+
  #   theme_bw()+
  #   theme(panel.grid =element_blank(),panel.border = element_blank(),
  #         axis.line.y = element_line(size=0.5, colour = "black"),
  #         axis.ticks.y=element_line(size=0.5, colour = "black"),
  #         axis.text.y=element_text(size=15,color="black"),
  #         axis.title.y=element_text(size=16,color="black"),
  #         axis.ticks.x=element_line(size=0.5, colour = "black"),
  #         axis.title.x=element_blank(),
  #         axis.line.x = element_line(size=0.5, colour = "black"),
  #         axis.text.x=element_blank(),#element_text(size=15,color="black"),
  #         plot.margin = unit(c(0,1,1,0), "cm"))+ #top,right,bottom,left    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
  #   # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
  #   #geom_line(aes(density, mean), size=1.2, col="red") +
  #   #geom_vline(xintercept = 30) +
  #   scale_x_continuous(limits=c(0,860000000), breaks = seq(0,860000000,200000000) )+
  #   #   scale_y_continuous(limits=c(2.5,4.5), breaks=c(2.5,3.5,4.5),labels = c(2.5,3.5,4.5)) +
  #   ylab(expression(H["all"]^"'"-H["CN"]^"'")) + xlab("")
  # # # ggsave(paste0(CHR, ".logscale.pdf"), p, height = 2.5, width = 4)
  # # chr_num=chr_num+1
  # tmplist[[chr_num]] <- p
  # chr_num=chr_num+1
  DF <- read.table(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/",chrID,"_CNC_entropy",sep=""),as.is = T, header = F, comment.char = "")

  
  #DF[,6] <- DF[,4]-DF[,5]
  #DF[,6][which(DF[,6] > 4)]= 4
  #DF[,6][which(DF[,6] < 0.5)] = 0.5
  #colnames(DF) <- c("chr","start1","end1","pi1","pi2","pi3")
  colnames(DF) <- c("chr","start1","end1","pi1")
  p <- ggplot(DF)+
    #geom_ribbon(aes(start1,ymin=0,ymax=num), fill="red", col="red")+
    geom_line(aes(start1,pi1), col="#8FCA86",size=0.5)+
    #geom_line(aes(start1,dup), col="red",size=1)+
    #geom_line(aes(start1,zero), col="black",size=1)+
    #cowplot::theme_cowplot()+
    #geom_text(data = qrc_df,aes(V2*1000000,2,label=V4),angle=45,hjust = 0,size=3)+
    #geom_text(data = gene_df,aes(Position*1000000,1,label=Gene),angle=30,hjust = 0,size=3)+
    geom_vline(xintercept = gene_df$Position*1000000,linetype=2,size=0.3)+
    annotate("rect", xmin = qrc_df$V2*1000000, xmax = qrc_df$V3*1000000, ymin = -Inf, ymax = 4, alpha = 0.5,
             fill = "#22BBEE")+
    theme_bw()+
    theme(panel.grid =element_blank(),panel.border = element_blank(),
          axis.line.y = element_line(size=0.5, colour = "black"),
          axis.ticks.y=element_line(size=0.5, colour = "black"),
          axis.text.y=element_text(size=15,color="black"),
          axis.title.y=element_text(size=16,color="black"),
          axis.ticks.x=element_line(size=0.5, colour = "black"),
          axis.title.x=element_blank(),
          axis.line.x = element_line(size=0.5, colour = "black"),
          axis.text.x=element_blank(),#element_text(size=15,color="black"),
          plot.margin = unit(c(0,1,1,0), "cm"))+ #top,right,bottom,left    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
    #geom_vline(xintercept = 30) +
    scale_x_continuous(limits=c(0,860000000), breaks = seq(0,860000000,200000000) )+
    #scale_y_continuous(limits=c(0.5,4), breaks=c(1,2.5,4),labels = c(1,2.5,4)) +
    #ylab(expression(H["all"]^"'"-H["CN"]^"'")) + 
    ylab(expression(H["CNC"]^"'")) + 
    xlab("")
  # # ggsave(paste0(CHR, ".logscale.pdf"), p, height = 2.5, width = 4)
  # chr_num=chr_num+1
  tmplist[[chr_num]] <- p
  chr_num=chr_num+1
  # p <- ggplot(DF)+
  #   #geom_ribbon(aes(start1,ymin=0,ymax=num), fill="red", col="red")+
  #   geom_line(aes(start1,pi2), col="#8FCA86",size=0.5)+
  #   #geom_line(aes(start1,dup), col="red",size=1)+
  #   #geom_line(aes(start1,zero), col="black",size=1)+
  #   #cowplot::theme_cowplot()+
  #   #geom_text(data = qrc_df,aes(V2*1000000,2,label=V4),angle=45,hjust = 0,size=3)+
  #   #geom_text(data = gene_df,aes(Position*1000000,1,label=Gene),angle=30,hjust = 0,size=3)+
  #   geom_vline(xintercept = gene_df$Position*1000000,linetype=2)+
  #   #annotate("rect", xmin = qrc_df$V2*1000000, xmax = qrc_df$V3*1000000, ymin = -Inf, ymax = 4, alpha = 0.5,
  #   #         fill = "#22BBEE")+
  #   theme_bw()+
  #   theme(panel.grid =element_blank(),panel.border = element_blank(),
  #         axis.line.y = element_line(size=0.5, colour = "black"),
  #         axis.ticks.y=element_line(size=0.5, colour = "black"),
  #         axis.text.y=element_text(size=15,color="black"),
  #         axis.title.y=element_text(size=16,color="black"),
  #         axis.ticks.x=element_line(size=0.5, colour = "black"),
  #         axis.title.x=element_text(size=16,color="black"),
  #         axis.line.x = element_line(size=0.5, colour = "black"),
  #         axis.text.x=element_blank(),#element_text(size=15,color="black"),
  #         plot.margin = unit(c(0,1,1,0), "cm"))+ #top,right,bottom,left    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
  #   # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
  #   #geom_line(aes(density, mean), size=1.2, col="red") +
  #   #geom_vline(xintercept = 30) +
  #   scale_x_continuous(limits=c(0,860000000), breaks = seq(0,860000000,200000000) )+
  #   #   scale_y_continuous(limits=c(2.5,4.5), breaks=c(2.5,3.5,4.5),labels = c(2.5,3.5,4.5)) +
  #   ylab(expression(H["all"]^"'"-H["CN"]^"'")) + xlab(chrID)
  # # # ggsave(paste0(CHR, ".logscale.pdf"), p, height = 2.5, width = 4)
  # # chr_num=chr_num+1
  # tmplist[[chr_num]] <- p
  # chr_num=chr_num+1
}

pall <- cowplot::plot_grid(plotlist = tmplist,ncol = 3,byrow=F,align = 'v')
pdf("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/pdf/H_cnv_CNC.pdf", height = 35, width = 25)

pall
dev.off()


