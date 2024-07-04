setwd("/home/pub252/users/liy/20240129_Pre-B_leukemia_scRNA/04_B_cells_TF")
##
AUCmatrix <- readRDS("3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
saveRDS(AUCmatrix,'scRNAauc.rds')

###
selected_columns <- !grepl("extended", colnames(AUCmatrix))
AUCmatrix <- AUCmatrix[,selected_columns]
write.table(AUCmatrix,file = "AUCmatrix.txt",sep = "\t",col.names = T,row.names = T,quote = F)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
# 
df_long <-AUCmatrix %>%
  gather(key = "TF", value = "TF_AUC")
#
df_long <- df_long[df_long$TF_AUC!=0,]
# 
average_auc <- aggregate(TF_AUC ~ TF, data = df_long , FUN = mean)

# 
average_auc_sorted <- average_auc[order(-average_auc$TF_AUC), ]

# 
df_long$TF <- factor(df_long$TF, levels = average_auc_sorted$TF)
df_long <- df_long[df_long$TF_AUC!=0,]
ggplot(df_long, aes(x =TF , y = TF_AUC)) + 
  
  geom_violin(fill = "#99CCFF",
              color = 'grey', alpha = 0.8,
              scale = 'width',
              linewidth = 0.6, #
              trim = TRUE) + 
  geom_boxplot( color = 'white',
                outlier.color = 'black',
                width = 0.4, #
                size = 0.8, #
                fill = NA, outlier.size = 0)+
  theme_bw() +
  rotate_x_text(angle = 30)
  #geom_point(aes())+
  #facet_wrap(~cell_type, scales = "free") + 
  # labs(title = "Gene Expression by Cell Type and Disease Status",
  #      x = "Gene",
  #      y = "Expression Level") +
  # scale_fill_manual(values = c("pre" = "#0099CC", "sepsis" = "#FF9966"))+
  #stat_compare_means(aes(group = TF), label = "p.format",method ='kruskal' )+
  
ggsave('TF_B_cells.pdf',width = 8,height = 5)

