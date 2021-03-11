files = read.table("/data/user/yangzz/mapping/fieldergenomecompare/ppd1/CNL_CNC_FC_FL_hap.txt", as.is = T, header = T, comment.char = "")
files$hap = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17')
files$CNC = files$CNC
files$CNL = files$CNL
files$FC = files$FC
files$FL = files$FL
files_melt <- melt(files,measure.vars = c('CNL','FC','FL','CNC'),id.vars = c("hap"))

col_lis <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/24_distinct_colors.txt", as.is = T, header = F, comment.char = "")
files$col <- col_lis$V1[1:17]
tmplist <- list()
chr_num=1
name_lis <- c("CNC","CNL","FC","FL")
for (line in seq(1,4)) {
  files_tmp <- files[which(files[,line]!=0),][,c(line,5,6)]
  colnames(files_tmp) <- c('num','hap','col')
  #files_tmp = files_tmp[order(files_tmp[,1], decreasing = TRUE),]
  files_tmp=files_tmp %>% mutate(prop=round(num/sum(num),2)) %>%
    mutate(type=factor(hap,levels=hap))
  
  p <- ggplot(files_tmp, aes(x="", y=prop, fill=type)) +
    geom_bar(stat="identity", width=1, color="white")+ 
    coord_polar("y", start=0) + theme_void()+
    scale_fill_manual(values=files_tmp$col,
                      name="hap type",
                      labels=sprintf("%s",files_tmp$hap))+
    geom_text(aes(x=1.6,label = num),
              position = position_stack(vjust = 0.5))+
    #geom_text(aes(label = hap), 
    #          position = position_stack(vjust = 0.5))+
    #scale_fill_discrete()+
    ggtitle(name_lis[line])+
    theme(plot.title = element_text(hjust = 0.5))
  #files[which(files$CNL!=0),]%>%
  #  mutate(hap = factor(hap, levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16')),
  #         cumulative = cumsum(CNL),
  #         midpoint = cumulative-,
  #         label = paste0("hap",hap, " ", round(CNL / sum(CNL) * 100, 1), "%"))
  
  
  tmplist[[chr_num]] <- p
  chr_num=chr_num+1
  
}
pall <- cowplot::plot_grid(plotlist = tmplist,byrow=F,align = 'v')
pall
#geom_line(aes(x=hap,y=CNC),col='#3171B8')+
#geom_line(aes(x=hap,y=CNL),col="#EB746A")+
#geom_line(aes(x=hap,y=FL),col="#BC0C3E")+
#geom_line(aes(x=hap,y=FC),col="#59AF57")+

library(dplyr)
data=data %>% mutate(prop=round(num/sum(num),2))

