library(Seurat)
library(doParallel)
library(Seurat)
library(SCopeLoomR)
setwd("/home/pub252/users/liy/20240129_Pre-B_leukemia_scRNA/04_SCENIC")
b_cells_1 <- readRDS()
exprMat <- as.matrix(b_cells_1@assays$RNA@counts)
write.table(exprMat,file = "exprMat.txt",sep = "\t",col.names = T,row.names = T)
cellInfo = b_cells_1@meta.data[, c(16,2,3)]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
write.table(cellInfo,file = "metadata.txt",sep = "\t",col.names = T,row.names = T,quote = F)

exprMat <- as.matrix(read.table("exprMat.txt",sep="\t",header = T,row.names = 1))

### Initialize settings
library(SCENIC)
# 
# 
mydbDIR <- "./database"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather", "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="hgnc", dbDir=mydbDIR, dbs=mydbs, nCores=4, datasetTitle="Epithelial cells") 
saveRDS(scenicOptions, "int/scenicOptions.rds")

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene=3*0.01*ncol(exprMat), minSamples=ncol(exprMat)*0.01)
exprMat_filtered <- exprMat[genesKept, ]
## 
runCorrelation(exprMat_filtered, scenicOptions)
## TF-Targets
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts=20)
## 

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top10perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") 

