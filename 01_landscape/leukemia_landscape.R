library(Seurat)
library(tidydr)
library(dplyr)
library(stringr)
#BiocManager::install("hdf5r")    #
library(Seurat) 
library(ggplot2)
library(Matrix)
library(SingleCellExperiment)
library(cols4all)
library(harmony)
library(dplyr)
library(stringr)


scRNAlist <- list()
file <- list.files()
dir <- paste0("./",file)
samples_name = file
for(i in 1:length(file ) ) {

  counts <- Read10X(data.dir=dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])  
   
    #
    if(T){ 
      scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
    }
    #
    if(T){
      scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
    }
    #
    if(T){
      HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
      HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
      scRNAlist[[i]][["percent.HB"]]=PercentageFeatureSet(scRNAlist[[i]],features=HB.genes)

  }
}
### 
names(scRNAlist) <- samples_name
# 
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
c4a_gui()
colors = c4a('tableau.superfishel_stone',9)
p = VlnPlot(scRNA, features=c("nFeature_RNA","nCount_RNA",'percent.mt'), pt.size=0.0, cols=colors)
scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<5000&percent.mt<10)  #
ggsave("after_nFeature_count_mt.pdf",height = 6,width = 18)
#####################################################################################
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
scRNA <- CellCycleScoring(scRNA, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
scRNA = SCTransform(scRNA, vars.to.regress=c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"), verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
scRNA = RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(scRNA,ndims = 50)
scRNA <- FindNeighbors(scRNA, dims=1:20, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:20, reduction="harmony")
p = DimPlot(scRNA, reduction="umap", group.by="orig.ident", pt.size=1)+
  theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("UMAP_Sample.pdf", p, height=6, width=6)
###################################### 
c4a_gui()
colors <- sample(c4a('poly.sky24',24),size = 20,replace = F)
mydata <- FindClusters(scRNA, resolution=0.1)
UMAPPlot(mydata, pt.size=1, label=T, cols=colors, label.size=5)+NoLegend()
ggsave("cluster_umap.pdf",height = 10,width = 10)
markers <- FindAllMarkers(mydata,group_by=cell_type,only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "./All_celltype_DEG_cell_type.txt", col.names=T, row.names=T, quote=F, sep="\t")
table(mydata$seurat_clusters)
# 0    1    2    3    4    5    6    7    8 
# 4199 3949 2299 1782 1105  707  512  507  114 
############## ##################

VlnPlot(mydata, features=c("SPINK2"), pt.size=0, cols=colors,group.by ='seurat_clusters')+NoLegend()+theme(axis.title.x=element_blank())
#0: Naive T cells: TCF7,  CD3E, CD3D
#1: B cells: PCDH9,QRSL1,TLE1, MS4A1
#2: Hematopoietic Stem cells 1: KLF1,AURKB,GTSE1
#3: Hematopoietic Stem cells 2: KLF1,GATA1,SOX6
#4: Myeloid cells: CDA,TREM1,MPO
#5: Hematopoietic Stem cells 3: SPINK2,CD34
#6: Cytotoxic NK/T cells: GZMH ,PRF1,GNLY,KLRD1,GZMA
#7: Plasma B cells: TNFRSF13B,ARHGAP24,FCRLA,IGHG3
#8: pDCs: CLEC4C,ZFAT,MDFIC


#mydata = subset(mydata, seurat_clusters %in% c(0,1,2,3,4,6,7,8))
#mydata@meta.data$seurat_clusters = droplevels(mydata@meta.data$seurat_clusters, exclude=c(5, 8, 9, 11, 12, 13))
# ################
cell_label = c("Naive T cells", "B cells ", "Hematopoietic Stem cells 1", "Hematopoietic Stem cells 2", "Myeloid cells",
               "Hematopoietic Stem cells 3", "Cytotoxic NK/T cells", "Plasma B cells","pDCs")
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
c4a_gui()
colors = c4a('carto.pastel',11)
p = UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)
ggsave("UMAP_subcluster.pdf", p, width=6, height=6)
genes = c('TCF7', 'CD3E', 'CD3D',
          'PCDH9','QRSL1','TLE1',
          'AURKB','GTSE1',
          'KLF1','GATA1','SOX6',
          'CDA','TREM1','MPO',
          'SPINK2','CD34',
          'GZMH','PRF1','GNLY',
          'TNFRSF13B','FCRLA','IGHG3',
         'CLEC4C','ZFAT','MDFIC')
p = DotPlot(mydata, features=genes)+coord_flip()+scale_color_gradientn(colors=c("black", "dodgerblue", "white", "orange", "firebrick1"))+theme_minimal()+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("dotplot_gene_marker.pdf",p,height = 6,width = 8)
#
mydata@meta.data <- mydata@meta.data %>% 
  mutate(type= str_extract(orig.ident, "^[^_]*"))
mydata$type <- factor(mydata$type,levels = c('HHD','PBMMC'))
# ##############################################################################
# top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# DoHeatmap(mydata, features = genes,slot ='scale.data',group.colors=colors,label = F) + NoLegend()+
#  scale_fill_gradientn(colors=c("white", "lightgray", "red"))
#########marker#####################################################
# set = c("ICOS", "MZB1", 'CD79A',"CD8A", "S100A12", "GNLY", "HBB", "CPA3")
# FeaturePlot(mydata, features=set, cols=c("snow", "purple"), ncol=4)
# ggsave('FeaturePlot.pdf',height = 4,width = 10)
gene =  c('TCF7',
                  'PCDH9',
                  'AURKB',
                  'KLF1',
                  'CDA',
                  'SPINK2',
                  'GZMH',
                  'IGHG3',
                  'CLEC4C')
VlnPlot(mydata, features=gene, pt.size=0, cols=colors, ncol=3)+NoLegend()+theme(axis.title.x=element_blank())
ggsave('violin_Plot.pdf',height = 12,width =10)

########### ######################################################## 
bar <-  as.data.frame(with(mydata@meta.data, table(type, cell_type)))
ggplot(data=bar, aes(x=type, y=Freq, fill=cell_type))+ 
  geom_bar(stat="identity", position=position_fill())+
  scale_fill_manual(values=colors)+theme_classic()+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12, angle=15, hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")
ggsave("cell_type_number.pdf",  width=6, height=6)
####################################################################  
Type_label = c("PBMMC", "HHD")
bar$type = factor(bar$type, levels=Type_label)
bar = bar %>% group_by(type) %>% mutate(percent=100*Freq/sum(Freq))

ggplot(data=bar, aes(x=cell_type, y=percent, fill=type,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c('#009999','#FF6600'))+theme_classic()+
  ggtitle("Percent(%)")+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
ggsave("barplot_pair_number3.pdf",  width=10, height=5)
###############################################################
######HHD
#install.packages('lessR')
library(lessR)
mata.data.HHD <- subset(mydata@meta.data,type=="HHD")
celltype.HHD <-  as.data.frame(mata.data.HHD[,'cell_type'])
colnames(celltype.HHD) <-'cell_type' 
rownames(celltype.HHD) <-rownames(mata.data.HHD)
pdf("HHD_PieChart.pdf",height =6,width = 10)
PieChart(cell_type, data = celltype.HHD,
         hole =0.5,
         main="Cell Proportion of HHD",
         main_cex =1.3,fill=colors)
dev.off()
#######PBMMC
mata.data.PBMMC <- subset(mydata@meta.data,type=="PBMMC")
celltype.PBMMC <-  as.data.frame(mata.data.PBMMC[,'cell_type'])
colnames(celltype.PBMMC) <-'cell_type' 
rownames(celltype.PBMMC) <-rownames(mata.data.PBMMC)
pdf("PBMMC_PieChart.pdf",height =6,width = 10)
PieChart(cell_type, data = celltype.PBMMC,
         hole =0.5,
         main="Cell Proportion of PBMMC",
         main_cex =1.3,fill=colors)
dev.off()
saveRDS(mydata,file ="mydata.RDS" )
subcell_type <- 'B cells '
mydata.B.cells <- subset(mydata,cell_type==subcell_type)
saveRDS(mydata.B.cells,file ='mydata.B.cells.RDS' )

