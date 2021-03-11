specificity <- c(100, 99.9, 99.9, 99.9,99.9, 100, 100, 99.9)
Group <- c(rep('deletion',4),rep('duplication',4))
percentage <- c(90, 60, 30, 10, 90, 60, 30, 10)
df <- data.frame(percentage = percentage, specificity=specificity, Group = Group)
ggplot(data = df, mapping = aes(x = percentage, y = specificity, colour = Group)) + geom_line()+
  scale_y_continuous(breaks = seq(99.8,100,0.1))+scale_x_continuous(breaks=c(10, 30, 60,90))+
  theme_bw()+scale_color_manual(values = c('steelblue','darkred'))+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))

sensitivity <- c(98.8, 99.8, 100, 98.3,100, 100, 96.6, 100)
Group <- c(rep('deletion',4),rep('duplication',4))
percentage <- c(90, 60, 30, 10, 90, 60, 30, 10)
df <- data.frame(percentage = percentage, sensitivity =sensitivity , Group = Group)


ggplot(data = df, mapping = aes(x = percentage, y = sensitivity , colour = Group)) + geom_line()+
  scale_y_continuous(limits = c(96,100))+
  #scale_y_continuous(breaks = c(95,96,97,98,99,100))+
  scale_x_continuous(breaks=c(10, 30, 60,90))+
  theme_bw()+scale_color_manual(values = c('steelblue','darkred'))+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))
