mydata.B.cells <- readRDS("Z:/users/liy/20240129_Pre-B_leukemia_scRNA/01_landscape/mydata.B.cells.RDS")
dir.create("Z:/users/liy/20240129_Pre-B_leukemia_scRNA/02_B_cells_subcluster")
setwd("Z:/users/liy/20240129_Pre-B_leukemia_scRNA/02_B_cells_subcluster")
library(Seurat)
library(tidydr)
library(dplyr)
library(stringr)
#BiocManager::install("hdf5r")    #
library(hdf5r)
library(Seurat) 
library(ggplot2)
library(Matrix)
library(SingleCellExperiment)
library(cols4all)
library(harmony)
library(dplyr)
library(stringr)
#####################################################################################
scRNA <- mydata.B.cells
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
scRNA <- CellCycleScoring(scRNA, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
scRNA = SCTransform(scRNA, vars.to.regress=c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"), verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
scRNA = RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(scRNA,ndims = 50)
scRNA <- FindNeighbors(scRNA, dims=1:20, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:20, reduction="harmony")
colors = c4a('tableau.superfishel_stone',9)
p = DimPlot(scRNA, reduction="umap", group.by="orig.ident", pt.size=1.8)+
  theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("UMAP_Sample.pdf", p, height=6, width=6)
##################################### 
c4a_gui()
colors <- sample(c4a('brewer.dark2',7),size = 5,replace = F)
mydata <- FindClusters(scRNA, resolution=0.3)
UMAPPlot(mydata, pt.size=1.6, label=T, cols=colors, label.size=5)+NoLegend()
ggsave("cluster_umap_03.pdf",height = 10,width = 10)
markers <- FindAllMarkers(mydata,group_by=cell_type,only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "./All_celltype_DEG_03.txt", col.names=T, row.names=T, quote=F, sep="\t")
table(mydata$seurat_clusters)
############## ##################
VlnPlot(mydata, features=c('TNFAIP3','KDM5B','IGLL5',"CXXC5"), pt.size=0, cols=colors,group.by ='seurat_clusters')+NoLegend()+theme(axis.title.x=element_blank())
#0:B cells 1 :TNFAIP3,KDM5B
#1:B cells 2:
#2:B cells 3:IGLL5,CXXC5
#3:B cells 4
#4:B cells 5
####小提琴图的gene marker
VlnPlot(mydata, features=c('TNFAIP3','KDM5B','IGLL5',"CXXC5"), pt.size=0, cols=colors,group.by ='cell_type',ncol = 2)+NoLegend()+theme(axis.title.x=element_blank())
#mydata = subset(mydata, seurat_clusters %in% c(0,1,2,3,4,6,7,8))
#mydata@meta.data$seurat_clusters = droplevels(mydata@meta.data$seurat_clusters, exclude=c(5, 8, 9, 11, 12, 13))
# 
cell_label = c('B cells 1','B cells 2','B cells 3','B cells 4','B cells 5')
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
#c4a_gui()
#colors = c4a('carto.pastel',11)
p = UMAPPlot(mydata, pt.size=1.3, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)
ggsave("UMAP_subcluster.pdf", p, width=6, height=6)
# genes = c()
# p = DotPlot(mydata, features=genes)+coord_flip()+scale_color_gradientn(colors=c("black", "dodgerblue", "white", "orange", "firebrick1"))+theme_minimal()+
#   theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
# ggsave("dotplot_gene_marker.pdf",p,height = 6,width = 8)
###########################B####################################
setwd
#
BP <- read.csv(file = "GO_B_cells_1_subcluster_select.csv",header = T)
BP$description<-gsub("^.{1}.*?~", "", BP$Term)
BP  <- na.omit(BP )
BP<- BP[order(BP$Count,decreasing = F),]
BP$description<-levels(factor(BP$description))

#cluster_1 
BP1 <- read.csv(file = "GO_B_cells_2_subcluster_select.csv",header = T)
BP1$description<-gsub("^.{1}.*?~", "", BP1$Term)
BP1  <- na.omit(BP1 )
BP1<- BP1[order(BP1$Count,decreasing = F),]
BP1$description<-levels(factor(BP1$description))

#cluster_2 
BP2 <- read.csv(file = "GO_B_cells_3_subcluster_select.csv",header = T)
BP2$description<-gsub("^.{1}.*?~", "", BP2$Term)
BP2  <- na.omit(BP2 )
BP2<- BP2[order(BP2$Count,decreasing = F),]
BP2$description<-levels(factor(BP2$description))
#cluster_3 
BP3 <- read.csv(file = "GO_B_cells_4_subcluster_select.csv",header = T)
BP3$description<-gsub("^.{1}.*?~", "", BP3$Term)
BP3  <- na.omit(BP3 )
BP3<- BP3[order(BP3$Count,decreasing = F),]
BP3$description<-levels(factor(BP3$description))
#cluster_4 
BP4 <- read.csv(file = "GO_B_cells_5_subcluster_select.csv",header = T)
BP4$description<-gsub("^.{1}.*?~", "", BP4$Term)
BP4 <- na.omit(BP4 )
BP4<- BP4[order(BP4$Count,decreasing = F),]
BP4$description<-levels(factor(BP4$description))

p1<-ggplot(data = BP, mapping = aes(x=description,y=Count,fill=-1*log10(PValue))) +
  geom_bar(stat="identity")+ 
  coord_flip()+ #
  scale_fill_gradient(low="grey80",high =colors[1])+
  theme(legend.title = element_text(size = 15, face = 2))+
  theme(legend.key.size=unit(1,'cm'))+
  theme(legend.text = element_text( size = 15,face = 'bold'))+
  labs(x=NULL,y = "Count",title = "B cells 1")+
  theme(plot.title = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.x = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.y = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.text.x = element_text(size = 12, family = "myFont", color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 12, family = "myFont", color ="black", face = "bold"))+
  theme_bw()
p1
ggsave("GO_B_cells_1.pdf",p1,height = 5,width = 6)
p2<-ggplot(data = BP1, mapping = aes(x=description,y=Count,fill=-1*log10(PValue))) +
  geom_bar(stat="identity")+ 
  coord_flip()+ #
  scale_fill_gradient(low="grey80",high =colors[2])+
  theme(legend.title = element_text(size = 15, face = 2))+
  theme(legend.key.size=unit(1,'cm'))+
  theme(legend.text = element_text( size = 15,face = 'bold'))+
  labs(x=NULL,y = "Count",title = "B cells 2")+
  theme(plot.title = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.x = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.y = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.text.x = element_text(size = 12, family = "myFont", color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 12, family = "myFont", color ="black", face = "bold"))+
  theme_bw()
ggsave("GO_B_cells_2.pdf",p2,height = 6,width = 6)
p3<-ggplot(data = BP2, mapping = aes(x=description,y=Count,fill=-1*log10(PValue))) +
  geom_bar(stat="identity")+ 
  coord_flip()+ #
  scale_fill_gradient(low="grey80",high =colors[3])+
  theme(legend.title = element_text(size = 15, face = 2))+
  theme(legend.key.size=unit(1,'cm'))+
  theme(legend.text = element_text( size = 15,face = 'bold'))+
  labs(x=NULL,y = "Count",title = "B cells 3")+
  theme(plot.title = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.x = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.y = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.text.x = element_text(size = 12, family = "myFont", color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 12, family = "myFont", color ="black", face = "bold"))+
  theme_bw()
ggsave("GO_B_cells_3.pdf",p3,height = 6,width = 6)

p4<-ggplot(data = BP3, mapping = aes(x=description,y=Count,fill=-1*log10(PValue))) +
  geom_bar(stat="identity")+ 
  coord_flip()+ #
  scale_fill_gradient(low="grey80",high =colors[4])+
  theme(legend.title = element_text(size = 15, face = 2))+
  theme(legend.key.size=unit(1,'cm'))+
  theme(legend.text = element_text( size = 15,face = 'bold'))+
  labs(x=NULL,y = "Count",title = "B cells 4")+
  theme(plot.title = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.x = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.y = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.text.x = element_text(size = 12, family = "myFont", color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 12, family = "myFont", color ="black", face = "bold"))+
  theme_bw()
ggsave("GO_B_cells_4.pdf",p4,height = 6,width = 6)

p5<-ggplot(data = BP4, mapping = aes(x=description,y=Count,fill=-1*log10(PValue))) +
  geom_bar(stat="identity")+ 
  coord_flip()+ #
  scale_fill_gradient(low="grey80",high =colors[5])+
  theme(legend.title = element_text(size = 15, face = 2))+
  theme(legend.key.size=unit(1,'cm'))+
  theme(legend.text = element_text( size = 15,face = 'bold'))+
  labs(x=NULL,y = "Count",title = "B cells 5")+
  theme(plot.title = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.x = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.title.y = element_text(size = 15, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.text.x = element_text(size = 12, family = "myFont", color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 12, family = "myFont", color ="black", face = "bold"))+
  theme_bw()
ggsave("GO_B_cells_5.pdf",p5,height = 6,width = 6)

###########################
VlnPlot(mydata, features=c('TNFAIP3','KDM5B','IGLL5',"CXXC5"), pt.size=0, cols=colors,group.by ='cell_type',ncol = 2)+NoLegend()+theme(axis.title.x=element_blank())
ggsave("B_cells_subcluster_gene_marker.pdf",height = 6,width = 8)
####################################
bar <-  as.data.frame(with(mydata@meta.data, table(type, cell_type)))
ggplot(data=bar, aes(x=type, y=Freq, fill=cell_type))+ 
  geom_bar(stat="identity", position=position_fill())+
  scale_fill_manual(values=colors)+theme_classic()+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12, angle=0, hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")
ggsave("cell_type_number.pdf",  width=6, height=6)
