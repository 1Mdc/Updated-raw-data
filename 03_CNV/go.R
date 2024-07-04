####################################
BP <- read.csv(file = "GO_BP_chr6_B&cells_up_select.csv",header = T)
BP$description<-gsub("^.{1}.*?~", "", BP$Term)
BP  <- na.omit(BP )
BP<- BP[order(BP$Count,decreasing = F),]
BP$description<-levels(factor(BP$description))
p1<-ggplot(data = BP, mapping = aes(x=description,y=Count,fill=-1*log10(PValue))) +
  geom_bar(stat="identity")+ 
  coord_flip()+ #
  scale_fill_gradient(low="grey80",high ='#CCCC33')+
  theme(legend.title = element_text(size = 15, face = 2))+
  theme(legend.key.size=unit(1,'cm'))+
  theme(legend.text = element_text( size = 15,face = 'bold'))+
  labs(x=NULL,y = "Count",title = "chr6 amplification genes ")+
  theme(plot.title = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.x = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.y = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.text.x = element_text(size = 12, family = "myFont", color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 12, family = "myFont", color ="black", face = "bold"))+
  theme_bw()
p1
ggsave("GO_BP_chr6_3 &up gennes of B cells 1.pdf",p1,height = 5,width = 6)


######################################
BP <- read.csv(file = "GO_BP_chr21_B_select.csv",header = T)
BP$description<-gsub("^.{1}.*?~", "", BP$Term)
BP  <- na.omit(BP )
BP<- BP[order(BP$Count,decreasing = F),]
BP$description<-levels(factor(BP$description))
p2<-ggplot(data = BP, mapping = aes(x=description,y=Count,fill=-1*log10(PValue))) +
  geom_bar(stat="identity")+ 
  coord_flip()+ #
  scale_fill_gradient(low="grey80",high ='#CC9933')+
  theme(legend.title = element_text(size = 15, face = 2))+
  theme(legend.key.size=unit(1,'cm'))+
  theme(legend.text = element_text( size = 15,face = 'bold'))+
  labs(x=NULL,y = "Count",title = "chr21 amplification genes ")+
  theme(plot.title = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.x = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.y = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.text.x = element_text(size = 12, family = "myFont", color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 12, family = "myFont", color ="black", face = "bold"))+
  theme_bw()
p2
ggsave("GO_BP_chr21_3.pdf",p2,height = 5,width = 6)
