##draw hap types density
library(reshape2)

tmplist <- list()
chr_num=1
for (i in c('chr1A','chr2A','chr3A','chr4A','chr5A','chr6A','chr7A','chr1B','chr2B','chr3B','chr4B','chr5B','chr6B','chr7B','chr1D','chr2D','chr3D','chr4D','chr5D','chr6D','chr7D')){#seq(1,(2*nrow(karyotype)),2)){
  chrID <-i #chrID <- karyotype[i,1]#karyotype[ceiling(i/2),1]
  typ = 'CNC'
  files = read_delim(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/",chrID,"_",typ,"_hap_count.txt",sep=""), "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
  files = files[which(files$X4 != 'duplication'&files$X4 != 'deletion'),]
  wide_file <- dcast(files,X1+X2~X4,value.var = 'X3')
  col_len <- ncol(wide_file)
  for (i in seq(1,nrow(wide_file))) {
    
    wide_file[i,'5%'] = length(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])>5)+2])

    wide_file[i,'1%'] = length(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])>2 & as.numeric(wide_file[i,3:(col_len-1)])<=5)+2])          
 
    wide_file[i,'ex'] = length(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])<=2)+2]) + wide_file[i,col_len]
  }
  p1 <- ggplot(wide_file) +
    #geom_line(aes(x = X2/1000000, y = ex),col='#747070')+
    geom_line(aes(x =X2/1000000, y = `1%`),col='#84D9DC')+
    geom_line(aes(x = X2/1000000, y = `5%`),col='#F2675E')+
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
    scale_y_continuous(expand = c(0,0),limits = c(0,10),breaks = c(0,2,4,6,8,10))+ylab("Counts") + xlab(chrID)
  
    tmplist[[chr_num]] <- p1
    chr_num <- chr_num +1
    
    typ = 'CNL'
    files = read_delim(paste("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/",chrID,"_",typ,"_hap_count.txt",sep=""), "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
    files = files[which(files$X4 != 'duplication'&files$X4 != 'deletion'),]
    wide_file <- dcast(files,X1+X2~X4,value.var = 'X3')
    col_len <- ncol(wide_file)
    for (i in seq(1,nrow(wide_file))) {
      
      wide_file[i,'5%'] = length(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])>5)+2])
      
      wide_file[i,'1%'] = length(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])>2 & as.numeric(wide_file[i,3:(col_len-1)])<=5)+2])          
      
      wide_file[i,'ex'] = length(wide_file[i,which(as.numeric(wide_file[i,3:(col_len-1)])<=2)+2]) + wide_file[i,col_len]
    }
    p1 <- ggplot(wide_file) +
      #geom_line(aes(x = X2/1000000, y = ex),col='#747070')+
      geom_line(aes(x =X2/1000000, y = `1%`),col='#84D9DC')+
      geom_line(aes(x = X2/1000000, y = `5%`),col='#F2675E')+
      #scale_fill_manual(values = col_lis)+     
      theme_bw()+
      theme(panel.grid = element_blank(), panel.border  = element_blank(),          
            axis.line.y = element_line(size=0.5, colour = "black"),
            axis.ticks.y=element_line(size=0.5, colour = "black"),
            axis.text.y=element_text(size=15,color="black"),
            axis.text.x=element_text(size=16,color="black"),
            axis.title.y=element_text(size=16,color="black"),
            axis.ticks.x=element_line(size=0.5, colour = "black"),
            axis.title.x=element_text(size=16,color="black"),
            axis.line.x = element_line(size=0.5, colour = "black"),
            legend.position  = "none",plot.margin = unit(c(0,1,0,0), "cm"))+
      #annotate(x=-5, xend=-5, y=0, yend=60, colour="black", lwd=1, geom="segment") +
      # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
      #geom_line(aes(density, mean), size=1.2, col="red") +
      #geom_vline(xintercept = 30) +
      scale_x_continuous(expand = c(0.02,0))+
      scale_y_continuous(expand = c(0,0),limits = c(0,10),breaks = c(0,2,4,6,8,10))+ylab("Counts") + xlab(chrID)
    
    tmplist[[chr_num]] <- p1
    chr_num <- chr_num +1
}  

pall <- cowplot::plot_grid(plotlist = tmplist,ncol = 3,byrow=F,align = 'v',scale = 0.9)
pdf("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/pdf/CNC_CNL_hapnum.pdf", height = 35, width = 25)

pall
dev.off()
