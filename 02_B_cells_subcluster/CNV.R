# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("infercnv")

library(infercnv)#
count.matrix=readRDS("tumor_and_ref_count.Rds")
class(count.matrix)#
load("data/gene.loc.RData")#
head(gene.loc)
class(gene.loc)
cell.anno<-readRDS("all_cell_barcode.Rds")#
head(cell.anno)
#
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=count.matrix,gene_order_file=gene.loc,
                                     annotations_file=cell.anno,
                                     ref_group_names=c("Oligodendrocyte","Macrophage","T-cell"),
                                     delim = "\t",max_cells_per_group = NULL,
                                     min_max_counts_per_cell = c(100, +Inf),
                                     chr_exclude = c("chrX", "chrY", "chrM"))

infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir=tempfile(), 
                              cluster_by_groups=TRUE, 
                              denoise=TRUE,
                              HMM=TRUE)