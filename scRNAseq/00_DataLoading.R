# R script for scRNA-seq loading and preprocessing
library(Seurat)
library(harmony)
library(ggplot2)

source(file = "./Function.R")

#construct the Seurat object and preprocess
scRNA_metaData <- data.table::fread("scRNA_matrix_metadata/scRNA_metaData.csv")
scRNA_metaData <- as.data.frame(scRNA_metaData)
rownames(scRNA_metaData) <- scRNA_metaData[[1]]
scRNA_metaData <- scRNA_metaData[,-1]

scRNA_RawCount <- data.table::fread("scRNA_matrix_metadata/scRNA_RawCount.csv")
scRNA_RawCount <- as.data.frame(scRNA_RawCount)
colnames(scRNA_RawCount) <- gsub("^X", "", colnames(scRNA_RawCount))
rownames(scRNA_RawCount) <- scRNA_RawCount[[1]]
scRNA_RawCount <- scRNA_RawCount[,-1]

object_scRNA <- CreateSeuratObject(counts = scRNA_RawCount,
                                   meta.data = scRNA_metaData,
                                   project = "Trem1")

object_scRNA <- ScRNA_data_processing(object_scRNA,reduction = "pca",resolution = 0.5)

DotPlot(object_scRNA,features = c("Ptprc",
                                  "Cd3d","Cd3g",
                                  "Klrb1c","Ncr1",
                                  "Gata3","Klrg1",
                                  "Cd19","Cd79a","Jchain",
                                  "Siglech","Ccr9",
                                  "Itgae","Xcr1","H2-DMb2",
                                  "Adgre1","C1qc","Csf1r",
                                  "G0s2","S100a8",
                                  "Mcpt2","Tpsb2",
                                  "Dcn","Cxcl12",
                                  "Acta2",
                                  "Pecam1","Lyve1",
                                  "Epcam","Cdh1"))+coord_flip()

#separating data and re-clustering 
object_wt <- subset(object_scRNA,group %in% "WT")
object_wt <- ScRNA_data_processing(object_wt,reduction = "pca",resolution = 0.1)
saveRDS(object_wt,file = "object_wt.rds")

object_tumor <- subset(object_scRNA,tissue %in% "Tumor")
object_tumor <- ScRNA_data_processing(object_tumor,reduction = "pca",resolution = 0.1)
saveRDS(object_tumor,file = "object_tumor.rds")



