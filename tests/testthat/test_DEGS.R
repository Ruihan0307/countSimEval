
library(Seurat)
library(cluster)
library(scater)
library(SingleCellExperiment)
library(stats)
library(zellkonverter)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(scDesign3)

rm(list = ls())
source("R/calcDEGs.R")

#hierarchicell只能生成control组和case组
#splat无法指定细胞类型生成，无法评价
#scDD只能生成对应的两个condition，对于多celltype的无法评价聚类效果

origin_sce <- readH5AD("E:/Datas/PythonDatas/InHouse/origin/filtered_rna.h5ad")

simu_sce1 <- readH5AD("E:/Datas/PythonDatas/InHouse/simu/scDiffusion_sim_all_num_as_origin.h5ad")
#simu_sce2 <- readH5AD("E:/Datas/PythonDatas/MouseBrain/simu/scGAN_simu_MouseBrain.h5ad")
simu_sce3 <- readH5AD("E:/Datas/PythonDatas/InHouse/simu/scVI_simu_Inhouse.h5ad")

#simu_sce2 <- subset(simu_sce2, , cluster %in% c(0, 1))

colData(simu_sce2)
#修改顺序
table(origin_sce$cell_type)
table(simu_sce1$label)
#table(simu_sce2$cluster)
table(simu_sce3$celltype)

dim(origin_sce)
dim(simu_sce1)
#dim(simu_sce2)
dim(simu_sce3)
#新排序
#simu_sce1 <- simu_sce1[, order(colData(simu_sce1)$cell_ontology_class)]
#simu_sce2 <- simu_sce2[, order(colData(simu_sce2)$cell_ontology_class)]

#检验两个sce的细胞类型是否完全相同
#all(simu_sce1$cell_ontology_class == simu_sce2$cell_ontology_class)

origin_sce$group <- origin_sce$cell_type
simu_sce1$group <- simu_sce1$label
#simu_sce2$group <- simu_sce2$cluster
simu_sce3$group <- simu_sce3$celltype

table(origin_sce$group)
table(simu_sce1$group)
#table(simu_sce2$group)
table(simu_sce3$group)



assayNames(origin_sce)[assayNames(origin_sce) == "X"] <- "counts"
assayNames(simu_sce1)[assayNames(simu_sce1) == "X"] <- "counts"
#assayNames(simu_sce2)[assayNames(simu_sce2) == "X"] <- "counts"
assayNames(simu_sce3)[assayNames(simu_sce3) == "X"] <- "counts"



#名称记得修改
ddsList <- list(origin = origin_sce, scDiffusion = simu_sce1, scVI = simu_sce3)

#names(result$results)[names(result$results) == "scDesign2"] <- "ESCO"

result <- calculateDEGsParameters(ddsList)

#write.csv(result$results,"E:/R_Files/DEGs_parameters/Python_results/InHouse/InHouse_DEGs.csv" , row.names = FALSE)
#saveRDS(result,"E:/R_Files/DEGs_parameters/Python_results/InHouse/InHouse_DEGs.rds")
