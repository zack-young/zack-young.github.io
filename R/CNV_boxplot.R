D_CNV_du <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/D_duplication_number", as.is = T, header = F, comment.char = "")
A_CNV_du <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/A_duplication_number", as.is = T, header = F, comment.char = "")
B_CNV_du <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/B_duplication_number", as.is = T, header = F, comment.char = "")

D_CNV_de <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/subgenome_CNV_count/D_deletion_number", as.is = T, header = F, comment.char = "")
A_CNV_de <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/subgenome_CNV_count/A_deletion_number", as.is = T, header = F, comment.char = "")
B_CNV_de <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/subgenome_CNV_count/B_deletion_number", as.is = T, header = F, comment.char = "")
B_CNV_no1B1R <- read.table("/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/subgenome_CNV_count/B_no1Bs_deletion_number", as.is = T, header = F, comment.char = "")

A_CNV <- A_CNV_du
B_CNV <- B_CNV_du
D_CNV <- D_CNV_du


D_CNV[,1] <- (D_CNV[,1]*1000000*100)/3951074735
A_CNV[,1] <- (A_CNV[,1]*1000000*100)/4934891648
B_CNV[,1] <- (B_CNV[,1]*1000000*100)/5180314468
#B_CNV_no1B1R[,1] <- (B_CNV_no1B1R[,1]*1000000*100)/4941314468

D_CNV[,3]='D'
A_CNV[,3]='A'
B_CNV[,3]='B'
#B_CNV_no1B1R[,3]='B_de'
combine_CNV <- rbind(D_CNV,A_CNV,B_CNV)
#combine_CNV <- B_CNV
#combine_CNV[,3] ='B_de'
# combine_CNV_no1B1R <- combine_CNV[-c(grep(pattern = "S140",x = combine_CNV[,2]),
#                                      grep(pattern = "NZ05",x = combine_CNV[,2]),
#                                      grep(pattern = "S3331",x = combine_CNV[,2]),
#                                      grep(pattern = "S64",x = combine_CNV[,2]),
#                                      grep(pattern = "S135",x = combine_CNV[,2]),
#                                      grep(pattern = "S209",x = combine_CNV[,2]),
#                                      grep(pattern = "S67",x = combine_CNV[,2]),
#                                      grep(pattern = "S6554",x = combine_CNV[,2]),
#                                      grep(pattern = "S131",x = combine_CNV[,2]),
#                                      grep(pattern = "S212",x = combine_CNV[,2]),
#                                      grep(pattern = "S221",x = combine_CNV[,2])),]


combine_CNV <- rbind(D_CNV,A_CNV,B_CNV,B_CNV_no1B1R)

z.test2sam = function(a, b){
  n.a = length(a)
  n.b = length(b)
  var.a = var(a)
  var.b = var(b)
  zeta = (mean(a) - mean(b)) / (sqrt(var.a/n.a + var.b/n.b))
  return(zeta)
}

z.test2sam(A_CNV[,1],D_CNV[,1])
z.test2sam(A_CNV[,1],B_CNV[,1])
z.test2sam(B_CNV[,1],D_CNV[,1])
mean(D_CNV[,1])

D_test=combine_CNV_no1B1R[which(combine_CNV_no1B1R[,3]=="D"),][,1]
A_test=combine_CNV_no1B1R[which(combine_CNV_no1B1R[,3]=="A"),][,1]
B_test=combine_CNV_no1B1R[which(combine_CNV_no1B1R[,3]=="B"),][,1]



D_test=combine_CNV[which(combine_CNV[,3]=="D"),][,1]
A_test=combine_CNV[which(combine_CNV[,3]=="A"),][,1]
B_test=combine_CNV[which(combine_CNV[,3]=="B"),][,1]
B_no1B1R=combine_CNV[which(combine_CNV[,3]=="B_de"),][,1]

t.test(A_test,B_test,var.equal=T)
t.test(A_test,D_test,var.equal=T)
t.test(B_test,D_test,var.equal=T)
par(mar=c(4, 4, 2, 4))





pdf("/data/user/yangzz/mapping/fieldergenomecompare/pdf/ABD_duplication.pdf", height = 4, width = 5)
#ggplot(combine_deletion, aes(x=V3, y=V1)) + geom_boxplot(outlier.colour = NA,width = 0.5)+  theme(legend.position="right")+
#    labs(x="Subgenome", y = "Percentage") + scale_y_continuous(limits=c(0, 10))+theme_bw()+
#    theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
#          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))
ylims <- combine_CNV %>%
  group_by(V3) %>%
  summarise(Q1 = quantile(V1, 1/4), Q3 = quantile(V1, 3/4))
ylims <- mutate(ylims,V4=Q1-1.5*(Q3-Q1))
ylims <- mutate(ylims,V5=Q1+1.5*(Q3-Q1))
  #ungroup() %>%
  #get lowest Q1 and highest Q3
  #summarise(lowQ1 = min(Q1), highQ3 = max(Q3))
D_CNV=D_CNV[which(D_CNV$V1>ylims[[4]][3]&D_CNV$V1<ylims[[5]][3]),]
B_CNV=B_CNV[which(B_CNV$V1>ylims[[4]][2]&B_CNV$V1<ylims[[5]][2]),]
A_CNV=A_CNV[which(A_CNV$V1>ylims[[4]][1]&A_CNV$V1<ylims[[5]][1]),]
combine_CNV_fi <- rbind(D_CNV,A_CNV,B_CNV)

ggplot(combine_CNV_fi) + 
  geom_violin(aes(x=V3, y=V1),width = 0.75,colour = "black",fill="#C15065",adjust = 2)+  
  geom_boxplot(data=combine_CNV, aes(x=V3, y=V1),width = 0.075,colour = "black")+  theme(legend.position="right")+ #outlier.colour = NA
  labs(x="Subgenome", y = "Percentage") + scale_y_continuous(limits=c(0, 8), breaks=seq(0,8,2),labels = c("0","20%","40%","60%","80%"))+
  #scale_y_continuous(breaks =seq(0, 7,1))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15,colour = "black"),
        axis.text.y = element_text(size = 15,colour = "black"),
        axis.title.x = element_text(size = 15,colour = "black"),
        axis.title.y = element_text(size = 15,colour = "black"))
#C15065
#83A0CA
dev.off()
