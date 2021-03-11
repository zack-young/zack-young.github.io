library(ggplot2)#导入ggplot2包
a<-read.table(file = "C:/Users/lenovo/Desktop/sevenkind.txt", header=F,sep="\t",stringsAsFactors=F,na.strings = "")#读入数据
options(scipen=200)

#g<-sprintf("%s%s","Chr",unique(a[,1]))#选出染色体号，同时每个号码前加上Chr前缀，以备后面用来标记x轴
a[,6]<-a[,3]-a[,2]#提取每段长度
colors<-c("28","purple4","red1","gray","royalblue1","magenta","hotpink","white")#三种分段，分别标记的颜色为红蓝绿，自己可以随意选取。
p<-ggplot(a,aes(x=a[,1],y=a[,6],fill=a[,5],group=a[,4]))+geom_bar(stat = "identity",position = "stack",width = 0.5)+    #x轴选的是染色体号加减0.2的数据，以样本类别分组，以分段类别来填充。
  #geom_text(aes(label = a[,6]), size = 3, hjust = 0, vjust = 0,position = "stack",angle=45)+ #每个柱子的最上面添加注释样品类别信息
  theme_bw() + theme(panel.grid.major = element_blank(),
                     legend.position = 'right',legend.title=element_blank(), axis.text.x = element_text(hjust = 0.5,vjust = 0.5,size = 12),axis.text.y = element_text(hjust = 0.3,vjust = 0.5,size = 12))+#去背景，同时调整legend的位置，以及去除legend的title
  scale_fill_manual(values = colors)#+coord_fixed(ratio=2/1)#设置填充的颜色  position 可设置图例位置
p+labs(xlab(""))
p+labs(ylab(""))
p+scale_x_discrete(breaks = c(1:21),labels = g) 
p+coord_flip()#x y轴互换

#scale_y_continuous(labels = c(0,"200MB","400MB","600MB","800MB")) #设置y轴的刻度信息。


pdf(file = "result.pdf",height = 12,width = 15)
print(p)
dev.off()
png(filename = "result.png",type = "cairo",height = 700,width = 1000)
print(p)
dev.off() 
```