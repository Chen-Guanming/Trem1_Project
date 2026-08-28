# R script for scRNA downstream analysis
library(Seurat)
library(harmony)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(plyr)
library(ggalluvial)
library(ComplexHeatmap)

source(file = "./Function.R")

object_wt <- readRDS("object_wt.rds")
object_tumor <- readRDS("object_tumor.rds")

object_wt$number <- 1
object_tumor$number <- 1

color_major <- c("#CC6666", "#FA8072", "#CC9966", "#CCCC66", 
                 "#99CC66", "#66CC99", "#66CCCC", "#6699CC", 
                 "#6666CC", "#CC66CC", "#B83D7A", "#CCB399", 
                 "#996680", "#CAAFCA", "#CC9999")
color_major_tumor <- c("#fdc086", "#f2991a", "#e6cc4c", "#cf6eb8", 
                       "#beaed4", "#1a734c", "#7fc97f", "#cc051a",
                       "#3387b5", "#b2d199", "#9c21b8", "#be71af", 
                       "#72daf2", "#3abca5", "#b7a39f")
color_subtype <- c("#cca6bf", "#cc7eb1", "#824880", "#5b86ab", 
                   "#d1edcb", "#47885e", "#d5c666")
#Figure 1J
Idents(object_wt) <- "majortype"
DimPlot(object_wt,raster = F,label = T,cols = color_major,pt.size = 0.1,repel = T)+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

#Figure 1K
Idents(object_wt) <- "majortype"
FeaturePlot(object_wt,raster = F,label = F,cols = c("#f1e9f1","#762a83"),
            pt.size = 0.1,features = c("Trem1"))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

#Figure 1I&M
vio_plot <- VlnPlot(subset(object_wt,majortype %in% "Macrophage"),
                    group.by = "tissue",features = "Trem1",pt.size = 0)[[1]]$data
colnames(vio_plot) <- c("gene","group")

a <- ggplot(vio_plot,aes(group,gene,fill=group))+
  geom_violin(color=NA)+
  geom_boxplot(aes(color=group),fill="white",width=0.1,outlier.shape = NA)+
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black",size = 0.5),
        legend.position = "none")+
  scale_fill_manual(values = c("#7c8fb2", "#daba53"))+
  scale_color_manual(values = c("#7c8fb2", "#daba53"))+
  stat_compare_means(aes(group=group),
                     label="p.signif",
                     method =c("t.test"),
                     hide.ns = F,
                     vjust = 1,
                     hjust = -1.5)+
  ylab("Trem1 Expression Level")+
  labs(title = "Macrophage")+
  xlab("Tissue")

vio_plot <- VlnPlot(subset(object_wt,majortype %in% "Neutrophil"),
                    group.by = "tissue",features = "Trem1",pt.size = 0)[[1]]$data
colnames(vio_plot) <- c("gene","group")

b <- ggplot(vio_plot,aes(group,gene,fill=group))+
  geom_violin(color=NA)+
  geom_boxplot(aes(color=group),fill="white",width=0.1,outlier.shape = NA)+
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black",size = 0.5),
        legend.position = "none")+
  scale_fill_manual(values = c("#7c8fb2", "#daba53"))+
  scale_color_manual(values = c("#7c8fb2", "#daba53"))+
  stat_compare_means(aes(group=group),
                     label="p.signif",
                     method =c("t.test"),
                     hide.ns = F,
                     vjust = 1,
                     hjust = -1.5)+
  ylab("Trem1 Expression Level")+
  labs(title = "Neutrophil")+
  xlab("Tissue")

a+b

#Figure 3A
Idents(object_tumor) <- "majortype"
DimPlot(object_tumor,raster = F,label = T,cols = color_major_tumor,pt.size = 0.1,repel = T)+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black",size = 1), 
        panel.grid = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  annotate("text",x=-10,y=20,label="133,294 cells",size = 5)

#Figure 3B
data <- ddply(object_tumor@meta.data,"group1",transform,percent=1/sum(number)*100)
data <- ddply(data,.(group,majortype),summarize,percent=sum(percent))

ggplot(data,aes(group1,percent,fill=majortype,alluvium=majortype))+
  geom_alluvium(width=0.6)+
  geom_bar(width=0.6,stat = "identity",position = "stack")+
  scale_fill_manual(values = color_major_tumor)+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black",size = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(colour = "black"))+
  ylab("Percentage of frequence(%)")

#Figure 3G
Tumor_majortype_DEGs <- NULL
temp <- SplitObject(object_tumor,split.by = "majortype")
for (i in c(1:length(temp))) {
  
  Idents(temp[[i]]) <- "group"

  marker <- FindMarkers(temp[[i]],ident.1 = "Trem1-/-",ident.2 ="WT",logfc.threshold = 0)
  
  marker$ident1 <- "Trem1-/-"
  marker$ident2 <- "WT"
  marker$majortype <- names(temp[i])
  marker$gene <- rownames(marker)
  
  Tumor_majortype_DEGs <- rbind(Tumor_majortype_DEGs,marker)
}
saveRDS(Tumor_majortype_DEGs,file = "Tumor_majortype_DEGs.rds")

library(ggrepel)
df <- subset(Tumor_majortype_DEGs, p_val_adj<0.05)
df$label <- ifelse(abs(df$avg_log2FC)>0.4 & -log10(df$p_val_adj)>25,yes = df$gene,no =NA)

ggplot(df,aes(avg_log2FC,-log10(p_val_adj),color=majortype))+
  geom_point(size=0.8)+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black",size = 0.5),
        axis.text = element_text(colour = "black"))+
  geom_vline(xintercept =c(0) ,col ="black",lty ="dashed")+
  geom_hline(yintercept =c(-log10(0.05)) ,col ="grey",lty ="dashed")+
  scale_color_manual(values = color_major_tumor)+
  geom_text_repel(aes(label=label),max.overlaps = getOption("ggrepel.max.overlaps", default = 10000))

#Figure 3H
object_majortype <- SplitObject(object_tumor,split.by = "majortype")
for (i in c(1:length(object_majortype))) {
  object_majortype[[i]] <- ScRNA_data_processing(object_majortype[[i]],
                                                  reduction = "harmony",
                                                  group.by = "Mix")
}

Idents(object_majortype$Macrophage) <- "subtype"
DimPlot(object_majortype$Macrophage,raster = F,label = T, 
                                    cols = color_subtype,pt.size = 0.1,repel = T)+
                                    theme(panel.background = element_blank(), 
                                          panel.grid = element_blank(), 
                                          panel.grid.minor = element_blank(),
                                          axis.text = element_blank(),
                                          axis.ticks = element_blank())

#Figure 3I
DotPlot(object_majortype$Macrophage,features = c("Ly6c2", "Chil3", "Plac8", 
                                                 "Cxcl9", "Cxcl10", "Nos2", 
                                                 "Arg1", "Inhba", "Ccl3", 
                                                 "Il1r2", "Nr4a2", "Ifitm1", 
                                                 "Vcam1", "Mmp13", "Apoe", "Ccl8",
                                                 "Fabp5", "Cd36", "Lpl", "Trem2", "Gpnmb","Spp1", 
                                                 "Top2a", "Stmn1", "Mki67"),
        cols = c("#f1e9f1","#8a5199"))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,face="italic"),
              axis.title = element_blank(),
              panel.border = element_rect(fill = NA, colour = "black",linewidth = rel(1.5)),
              axis.line = element_blank())

#Figure 3M
library(irGSEA)
macrophage_H <- irGSEA.score(object = object_majortype$Macrophage, assay = "RNA",ncores = 1,
                               species = "Mus musculus", category = "H")
macrophage_AUCell_Hallmark <- AverageExpression(macrophage_H,assays = "AUCell",group.by = "subtype")$AUCell
macrophage_AUCell_Hallmark <- as.data.frame(scale(t(macrophage_AUCell_Hallmark)))

Heatmap(t(macrophage_AUCell_Hallmark[,c(36,17,9,19)]),
        cluster_rows = F,cluster_columns = F,
        row_names_side = "left",name = "AUcell score",
        col = c("#006200","#a6dba0","#fdfdfd","#c2a5cf","#8a5199"),
        border = TRUE)

#Figure 3N
geneset <- list()
geneset$LAM <- c("Cd63","Cd68","Cstb","Ctsb","Ctsd","Fabp5","Gpnmb","Lgals3",
                  "Lipa","Lpl","Pld3","Plin2","Spp1","Trem2")
object_majortype$Macrophage <- AddModuleScore(object_majortype$Macrophage,features = geneset[1],name = "LAM")
object_majortype$Macrophage$LAM_score <- object_majortype$Macrophage$LAM1

FeaturePlot(object_majortype$Macrophage,raster = F,label = T,repel = T,
            cols = c("#DCEEFA","#CD0000"),min.cutoff="q2",
            pt.size = 0.1,features = c("LAM_score"))+
            theme(axis.text = element_blank(),
                  axis.ticks = element_blank())

#Figure 4A and B
library(tidyr)
library(stringr)

#load data
object_velo <- object_majortype@Macrophage
object_velo_WT <- subset(object_majortype@Macrophage,group %in% "WT")
object_velo_KO <- subset(object_majortype@Macrophage,group %in% "Trem1-/-")

#overlapping loom file and seurat object
DefaultAssay(object = object_velo) <- "RNA"
files = c("./UB_Combined_WTA_CRC1_MNX0D.loom",
          "./UB_Combined_WTA_CRC2_40GMT.loom",
          "./UB_Combined_WTA_CRC3_V2G91.loom",
          "./UB_Combined_WTA_CRC4_XL8XA.loom",
          "./UB_Combined_WTA_CRC5_ONRSS.loom")

sample_tag <-  levels(object_velo$Mix)

for (i in seq(1, length(files))) {
  
  file = files[i]

  x <- velocyto.R::read.loom.matrices(file = file, engine =  "hdf5r")
  
  s = x[["spliced"]][intersect(rownames(object_velo), rownames(x[["spliced"]])), ]
  u = x[["unspliced"]][intersect(rownames(object_velo), rownames(x[["unspliced"]])), ]
  a = x[["ambiguous"]][intersect(rownames(object_velo), rownames(x[["ambiguous"]])), ]
  
  colnames(s) <- str_split_fixed(str_split_fixed(colnames(s),pattern = ":",n=2)[,2],pattern = "x",n=2)[,1]
  colnames(u) <- str_split_fixed(str_split_fixed(colnames(u),pattern = ":",n=2)[,2],pattern = "x",n=2)[,1]
  colnames(a) <- str_split_fixed(str_split_fixed(colnames(a),pattern = ":",n=2)[,2],pattern = "x",n=2)[,1]
  
  
  colnames(s) <- paste(colnames(s),sample_tag[i],sep = "_")
  colnames(u) <- paste(colnames(u),sample_tag[i],sep = "_")
  colnames(a) <- paste(colnames(a),sample_tag[i],sep = "_")
  
  if (length(intersect(colnames(s), colnames(object_velo))) == 0){
    warning("No cells in", file, " used!")
  }
  
  if (i == 1){
    spliced = s
    unspliced = u
    ambiguous = a
  } else{
    spliced = cbind(spliced, s)
    unspliced = cbind(unspliced, u)
    ambiguous = cbind(ambiguous, a)
  }
}

object.cells <- intersect(colnames(spliced), colnames(object_velo))
object.cells_WT <- intersect(colnames(spliced), colnames(object_velo_WT))
object.cells_KO <- intersect(colnames(spliced), colnames(object_velo_KO))
length(object.cells) #check overlap cell number!
length(object.cells_WT) #check overlap cell number!
length(object.cells_KO) #check overlap cell number!

if (length(object.cells) != length(colnames(object_velo))) {
  object_velo <- subset(object_velo, cell = object.cells)
}
if (length(object.cells_WT) != length(colnames(object_velo_WT))) {
  object_velo_WT <- subset(object_velo_WT, cell = object.cells_WT)
}  
if (length(object.cells_KO) != length(colnames(object_velo_KO))) {
  object_velo_KO <- subset(object_velo_KO, cell = object.cells_KO)
}    

if (length(object.cells) != length(colnames(spliced))) {
  spliced <- spliced[, object.cells]
  unspliced <- unspliced[, object.cells]
  ambiguous <- ambiguous[, object.cells]
}
if (length(object.cells_WT) != length(colnames(spliced))) {
  spliced_WT <- spliced_WT[, object.cells_WT]
  unspliced_WT <- unspliced_WT[, object.cells_WT]
  ambiguous_WT <- ambiguous_WT[, object.cells_WT]
}
if (length(object.cells_KO) != length(colnames(spliced))) {
  spliced_KO <- spliced_KO[, object.cells_KO]
  unspliced_KO <- unspliced_KO[, object.cells_KO]
  ambiguous_KO <- ambiguous_KO[, object.cells_KO]
}

object_velo[["RNA"]] <- CreateAssayObject(counts = object_velo@assays$RNA@counts[rownames(spliced),])
object_velo[["spliced"]] <- CreateAssayObject(counts = spliced)
object_velo[["unspliced"]] <- CreateAssayObject(counts = unspliced)
object_velo[["ambiguous"]] <- CreateAssayObject(counts = ambiguous)
object_velo_WT[["RNA"]] <- CreateAssayObject(counts = object_velo_WT@assays$RNA@counts[rownames(spliced_WT),])
object_velo_WT[["spliced"]] <- CreateAssayObject(counts = spliced_WT)
object_velo_WT[["unspliced"]] <- CreateAssayObject(counts = unspliced_WT)
object_velo_WT[["ambiguous"]] <- CreateAssayObject(counts = ambiguous_WT)
object_velo_KO[["RNA"]] <- CreateAssayObject(counts = object_velo_KO@assays$RNA@counts[rownames(spliced_KO),])
object_velo_KO[["spliced"]] <- CreateAssayObject(counts = spliced_KO)
object_velo_KO[["unspliced"]] <- CreateAssayObject(counts = unspliced_KO)
object_velo_KO[["ambiguous"]] <- CreateAssayObject(counts = ambiguous_KO)

Idents(object_velo) <- "subtype"
object_velo$subtype <- object_velo@active.ident
Idents(object_velo_WT) <- "subtype"
object_velo_WT$subtype <- object_velo_WT@active.ident
Idents(object_velo_KO) <- "subtype"
object_velo_KO$subtype <- object_velo_KO@active.ident

seurat2anndata(obj = object_velo, outFile = "object_velo.h5ad")
seurat2anndata(obj = object_velo_WT, outFile = "object_velo_WT.h5ad")
seurat2anndata(obj = object_velo_KO, outFile = "object_velo_KO.h5ad")

#Figure 4C and 4D
metadata <- readr::read_csv("./Dynamo_metadata.csv")
data <- object_majortype$Macrophage@meta.data[metadata$cell_id,]
data$Dynamo_pseudotime <- metadata$umap_ddhodge_potential
data$Divergence_score <- metadata$divergence_pca

ggplot(data,aes(Dynamo_pseudotime,Divergence_score,color=group))+
        geom_smooth(method = loess)+
        scale_color_manual(values = c("#4169E1","#CD2626"))+
        theme_bw()+
        theme(axis.text = element_text(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))

ggplot(data,aes(Dynamo_pseudotime,color=group))+
        geom_density()+
        scale_color_manual(values =  c("#4169E1","#CD2626"))+
        theme_bw()+
        theme(axis.text = element_text(colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black"))+
        ggplot(data,aes(Dynamo_pseudotime,number,fill=subtype))+
        geom_bar(width=0.01,stat = "identity",position = "stack")+
        scale_fill_manual(values = color_subtype)+
        theme_bw()+
        theme(axis.text = element_text(colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black"))

#Figure 4G
geneset <- list()
geneset$Lipid <- c("Abca1", "Abca4", "Abca2", "Abcg1", "Abcd1", "Fabp4", "Apoa1", 
                   "Apoa2", "Apoa4", "Apoe", "Atp8a1", "Atp9a", "Atp10a", "Cd36", 
                   "Fabp3", "Fabp2", "Fabp1", "Gm2a", "Fabp5", "Abcc1", "Mttp", 
                   "Npc1", "Pctp", "Abcb1b", "Abcb4", "Abcb1a", "Pitpna", "Pitpnm1", 
                   "Plscr2", "Pltp", "Abcd3", "Abcd4", "Rbp4", "Pitpnm2", "Scp2", 
                   "Slc10a1", "Slc10a2", "Slc2a1", "Star", "Stra6", "Plscr1", 
                   "Slco2a1", "Ceacam1", "Ceacam2", "Slc27a1", "Slc27a2", "Slc27a5", 
                   "Slc27a3", "Slc27a4", "Abcd2", "Abca7", "Abca8b", "Abcg5", "Abca3", 
                   "Abcb11", "Slco1a1", "Slco1a4", "Slco1b2", "Slco1a6", "Ttpa", 
                   "Atp8a2", "Atp11a", "Atp9b", "Gramd1a", "Atp8b2", "Atp8b1", "Apom", 
                   "Pitpnb", "Gltp", "Slc43a3", "Slco1c1", "Stard3", "Osbpl1a", 
                   "Prelid3b", "Prelid1", "Atp8b3", "Abcg8", "Abca14", "Npc2", "Cert1", 
                   "Triap1", "Tmem30a", "Plscr3", "Ano9", "Osbpl3", "C2cd2l", "Pitpnc1", 
                   "Osbp2", "Osbpl10", "Abca12", "Slc10a6", "Vmp1", "Abca6", "Atp11b", 
                   "Osbp", "Abcc3", "Atg2b", "Mfsd2a", "Slc10a7", "Prelid2", 
                   "1700057G04Rik", "Osbpl5", "Cptp", "Plekha3", "Osbpl6", "Osbpl9", 
                   "Slco2b1", "Ano6", "Osbpl11", "Slc51a", "Slco1a5", "Stard4", 
                   "Stard5", "Slc22a27", "Abcg4", "Slc22a19", "Gramd1c", "Slc10a3", 
                   "Slc5a8", "Gltpd2", "Spns2", "Abca8a", "Abca9", "Abca5", "Slc27a6", 
                   "Prelid3a", "Ano3", "Osbpl2", "Atp10d", "Slc10a4", "Plekha8", 
                   "Tmem41b", "Abca16", "Gramd1b", "Plscr4", "Slc22a26", "Slc22a29", 
                   "Osbpl8", "Npc1l1", "Apob", "Tmem30b", "Abcc4", "Atp8b4", "Slc10a5", 
                   "Tnfaip8l3", "Atg9a", "Abca13", "Atp10b", "Slc22a30", "Ano4", 
                   "Atp8b5", "Abca15", "Atp11c", "Atg2a", "Slc51b", "Abca17", "Xkr9", 
                   "Xkr8", "Ano7", "Mfsd2b", "Slc22a28", "Gm5724", "Xkr4", "Slc10a4-ps", 
                   "Gm6614", "Plscr5")
geneset$Fatty <- c("Acacb", "Acly", "Acaca", "Pecr", "Acsl1",  
                   "Fasn", "Hadh", "Pcx", "Scd1", "Acsl6", "Mecr",  
                   "Acsl5", "Acsl4", "Ech1", "Echdc2", "Acaa2", "Echdc1", 
                   "Acss2", "Decr1", "Echdc3", "Acsl3", "Echs1")

object_majortype$Macrophage <- AddModuleScore(object_majortype$Macrophage,
                                              features = geneset,name = names(geneset))
object_majortype$Macrophage$`Lipid transporter activity` <- object_majortype$Macrophage$Lipid1
object_majortype$Macrophage$`Fatty acid biosynthesis` <- object_majortype$Macrophage$Fatty2

vio_plot <- subset(object_majortype$Macrophage@meta.data,subtype %in% "TAM_Spp1")
vio_plot$group <- factor(vio_plot$group,levels = c("WT","Trem1-/-"))
a <- ggplot(vio_plot,aes(group,`Lipid transporter activity`,fill=group))+
  geom_violin(color=NA,alpha = 0.5)+
  geom_boxplot(aes(colour=group),width=0.3,
               fill="white",outlier.shape = NA)+
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.title = element_blank())+
  scale_fill_manual(values = c("#4169E1","#CD2626"))+
  scale_color_manual(values = c("#4169E1","#CD2626"))+
  stat_compare_means(aes(group=group),
                     label="p.signif",
                     method =c("t.test"),
                     hide.ns = F,
                     vjust = 1)+
  ylab("Lipid transporter activity Score")

b <- ggplot(vio_plot,aes(group,`Fatty acid biosynthesis`,fill=group))+
  geom_violin(color=NA,alpha = 0.5)+
  geom_boxplot(aes(colour=group),width=0.3,
               fill="white",outlier.shape = NA)+
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.title = element_blank())+
  scale_fill_manual(values = c("#4169E1","#CD2626"))+
  scale_color_manual(values = c("#4169E1","#CD2626"))+
  stat_compare_means(aes(group=group),
                     label="p.signif",
                     method =c("t.test"),
                     hide.ns = F,
                     vjust = 1)+
  ylab("Fatty acid biosynthesis Score")

#Figure 5A
Tumor_subtype_DEGs <- NULL
temp <- SplitObject(object_tumor,split.by = "subtype")

for (i in c(1:length(temp))) {
  Idents(temp[[i]]) <- "group1"

  marker <- FindMarkers(temp[[i]],ident.1 = "Trem1-/-",ident.2 ="WT",logfc.threshold = 0)
  marker$ident1 <- "Trem1-/-"
  marker$ident2 <- "WT"
  marker$subtype <- names(temp[i])
  marker$gene <- rownames(marker)
  
  Tumor_subtype_DEGs <- rbind(Tumor_subtype_DEGs,marker)
}

df <- subset(Tumor_subtype_DEGs,subtype %in% "TAM_Spp1")
df$group <- "nosig"
df$group <- ifelse(df$avg_log2FC<0 & df$p_val_adj < 0.05, yes = "WT",no = df$group)
df$group <- ifelse(df$avg_log2FC>0 & df$p_val_adj < 0.05, yes = "Trem1-/-",no = df$group)
df$label <- ifelse(df$gene %in% c("Gpnmb", "Fabp5", "Dusp1", "Fos", "Lpl", "Ctsb", "Jun",
                                  "Ctsd", "Trem2", "Psap", "Ctss", "Apoe", "Cd36"),
                                  yes = df$gene,no =NA)
df$group <- factor(df$group,levels = c("WT","Trem1-/-","nosig"))

ggplot(df,aes(avg_log2FC,-log10(p_val_adj),color=group))+
  geom_point(size=0.4)+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black",size = 0.5),
        axis.text = element_text(colour = "black"))+
  geom_vline(xintercept =c(0) ,col ="black",lty ="dashed")+
  geom_hline(yintercept =c(-log10(0.05)) ,col ="grey",lty ="dashed")+
  scale_color_manual(values = c("#4169E1","#CD2626","grey"))+
  geom_text_repel(aes(label=label),max.overlaps = getOption("ggrepel.max.overlaps", default = 10000))


#Figure S2B
Idents(object_wt) <- "majortype"
DotPlot(object_tumor,features = c("Ptprc",
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
                                  "Epcam","Cdh1"),cols = c("#f1e9f1","#8a5199"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        axis.text.y = element_text(face="italic"),
        axis.title = element_blank())+coord_flip()

#Figure S4A
Idents(object_tumor) <- "majortype"
DimPlot(object_tumor,raster = F,label = F,cols = color_major_tumor,pt.size = 0.1,split.by = "sample",ncol = 5)+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

#Figure S4C
data <- ddply(object_tumor@meta.data,"sample",transform,percent=1/sum(number)*100)
data <- ddply(data,.(group,subtype,sample),summarize,percent=sum(percent))

ggplot(subset(data, subtype %in% c("mMDSC_Chil3", "mMDSC_Cxcl9", "mMDSC_Ccl3", "TAM_Il1r2",
                                   "TAM_Vcam1", "TAM_Spp1", "TAM_Mki67")),
       aes(group,percent,fill=group))+
  geom_boxplot(outlier.shape = NA,alpha = 0.5,width=0.5)+
  geom_point(aes(group,percent,fill=sample),fill="black")+
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black",size = 0.5),
        legend.position = "none",
        axis.title.x = element_blank())+
  scale_fill_manual(values = c("#4169E1","#CD2626"))+
  stat_compare_means(aes(group=group),
                     label="p.format",
                     method =c("t.test"),
                     vjust = 1.5)+
  facet_wrap(~subtype, ncol = 7, scales = "free_y")+
  ylab("Percentage of frequence(%)")

#Figure S7A
##SCORPION analysis
library(SCORPION)
library(pbapply)
library(Matrix)
library(Rfast)

#loading TF and PPI data which were download in SCORPION github (https://github.com/kuijjerlab/SCORPION/tree/main)
load("./mm10_TF.RData")
load("./mm_PPI.RData")

#for TAM_Spp1 subtype
temp <- subset(object_tumor,subtype %in% "TAM_Spp1")
Mat_exp <- temp@assays$RNA@counts
output_dir <- paste("./result","TAM_Spp1",sep = "_")
  
for (i in 1:length(levels(temp$sample))) {
    
    outFile <- paste0(output_dir,"/",levels(temp$sample)[i], '.RData')
    
    sample_exp <- subset(temp@meta.data, sample %in% levels(temp$sample)[i])$cell_id
    sample_exp <- Mat_exp[,sample_exp]
    
    N <- scorpion(tfMotifs = mmTF, gexMatrix = sample_exp, ppiNet = mmPPI)[[1]]
    N <- round(N,3)
    
    gc()
    save(N,file = outFile)
}

loadNetwork <- function(X){
  load(X)
  return(N)
}

makeComparable <- function(X, tfList, gList){
  X <- as.matrix(X)
  O <- matrix(data = 0, nrow = length(tfList), ncol = length(gList), dimnames = list(tfList, gList))
  O[rownames(X), colnames(X)] <- X
  O <- Matrix(O)
  gc()
  return(O)
}

fileList <- list.files('./result_TAM_Spp1/',  full.names = TRUE)
donorList <- gsub('.RData',"", basename(fileList))
fileContent <- lapply(fileList, loadNetwork)
names(fileContent) <- donorList
gList <- unique(unlist(lapply(fileContent, colnames)))
tfList <- unique(unlist(lapply(fileContent, rownames)))
fileContent <- pblapply(fileContent, function(X){makeComparable(X = X, tfList = tfList, gList = gList)})
fileContent <- sapply(fileContent, function(X){reshape2::melt(as.matrix(X))[,3]})                           
edgeList <- expand.grid(tfList, gList)

sideList <- as.data.frame(table(object_majortype$sample,object_majortype$group))
sideList <- subset(sideList,Freq >0)
colnames(sideList) <- c("sample","group","Freq")
rownames(sideList) <- sideList$sample
sideList <- sideList[,2,drop=F]
sideList <- sideList[donorList,]

fileContent <- data.frame(edgeList, fileContent)
R <- fileContent %>% group_by(Var1) %>% summarise(across(colnames(fileContent)[-c(1:2)], sum))
C <- as.matrix(R[,2:ncol(R)])
C <- C[,!is.na(sideList)]
sideList <- sideList[!is.na(sideList)]
rownames(C) <- R$Var1

O <- ttests(t(C), ina = (sideList == 'WT')+1)
O <- as.data.frame(O)
O$FDR <- p.adjust(O$pvalue, method = 'fdr')
rownames(O) <- R$Var1
O <- O[order(O$stat, decreasing = TRUE),]
write.csv(O,file = "TAM_SPP1_DEtf.csv")