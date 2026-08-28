#ATAC-seq downstream analysis

library(DiffBind)
library(ggplot2)
library(dplyr)
library(JASPAR2024)
library(seqLogo)
library(motifmatchr)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SummarizedExperiment)


#1. Find differential accessibility peaks
samples <- data.frame(
  SampleID = c("WT1", "WT2", "KO1", "KO2"),
  Condition=c(rep("WT", 2), rep("KO", 2)),
  bamReads = c("./WT1.rmChrM.markDup.filtered.bam",
               "./WT2.rmChrM.markDup.filtered.bam",
               "./KO1.rmChrM.markDup.filtered.bam",
               "./KO2.rmChrM.markDup.filtered.bam"),
  Peaks = c("./WT1_peaks.narrowPeak",
            "./WT2_peaks.narrowPeak",
            "./KO1_peaks.narrowPeak",
            "./KO2_peaks.narrowPeak"),
  PeakCaller = "narrow"
)

object_DBA <- dba(sampleSheet=samples)
counts <- dba.count(object_DBA,minOverlap = 2,summits=250,
                    bRemoveDuplicates=FALSE, bScaleControl=FALSE, bSubControl = FALSE)

dba.plotPCA(counts, attributes=DBA_CONDITION, label=DBA_ID)
dba.plotMA(counts, bNormalized=FALSE, sub="Non-Normalized",
           contrast=list(WT=counts$masks$WT,
                         KO=counts$masks$KO))

nor_count <- dba.normalize(counts,normalize = DBA_NORM_RLE) 
dba.plotMA(nor_count, sub="Normalized (Default)",
           contrast=list(WT=counts$masks$WT,
                         KO=counts$masks$KO))

# Establishing a contrast
DA <- dba.contrast(nor_count, categories=DBA_CONDITION,minMembers = 2)
DA <- dba.analyze(DA, method=DBA_DESEQ2)



comp1.deseq <- dba.report(DA, method=DBA_DESEQ2, contrast = 1, th=1)
result_deseq <- as.data.frame(comp1.deseq)
result_deseq$site <- paste(result_deseq$seqnames, ":", result_deseq$start, "-", result_deseq$end, sep = "")

#Create bed files for peaks
write.table(result_deseq[, c("seqnames", "start", "end", "strand", "Fold")],
            file="./deseq2_all.bed", sep="\t", quote=F, row.names=F, col.names=F)

#bash (Annotation using Homer)
##annotatePeaks.pl deseq2_all.bed mm10 > deseq2_all.annotation.txt

deseq.bed.annotation_all <- read.delim(file = "./deseq2_all.annotation.txt")
colnames(deseq.bed.annotation_all)[c(2:4)] <- c("seqnames", "start", "end")
deseq.bed.annotation_all$start <- deseq.bed.annotation_all$start - 1

deseq.bed.annotation_all <- deseq.bed.annotation_all[, c(2:4, 8, 10, 16, 19)]
deseq.bed.annotation_all <- unique(deseq.bed.annotation_all)
final_results <- merge(result_deseq, deseq.bed.annotation_all, by = c("seqnames", "start", "end"))

#2. motif enrichment
opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["collection"]] <- "CORE"
opts[["species"]] <- "Mus musculus"
opts[["all_versions"]] <- FALSE

jaspar <- JASPAR2024::JASPAR2024()
motifsToScan <- getMatrixSet(jaspar@db, opts)

temp <- subset(final_results, FDR <0.05 & combination == c("KO2,KO3,WT2,WT3"))
peaks <- GRanges(seqnames = temp$seqnames,
                 ranges = IRanges(start = temp$start,
                                  width = 500))

motifHits <- matchMotifs(motifsToScan, peaks, genome = "mm10",out = "matches") 
mmMatrix <- motifMatches(motifHits)

df <- cbind(as.data.frame(temp), as.matrix(mmMatrix))
df$direction <- ifelse(df$Fold>0,yes = "up",no="down")
motif_cols <- grep("^MA", colnames(df), value = TRUE)

# test each motif
fold_diff <- c()
p_values <- c()
for (motif in motif_cols) {

  up_count <- sum(df[[motif]] & df$direction == "up")
  up_total <- sum(df$direction == "up")
  up_prop <- up_count / up_total
  
  down_count <- sum(df[[motif]] & df$direction == "down")
  down_total <- sum(df$direction == "down")
  down_prop <- down_count / down_total
  
  fold_diff <- c(fold_diff, up_prop - down_prop)

  table <- table(df[[motif]], df$direction)
  test <- fisher.test(table)
  p_values <- c(p_values, test$p.value)
}

results <- data.frame(
  motif = motif_cols,
  fold_diff = fold_diff,
  p_value = p_values
)

# FDR adjust
results <- results %>% mutate(adj_p_value = p.adjust(p_value, method = "fdr"))

#visualization
results <- results %>%
  mutate(
    significance = case_when(
      adj_p_value < 0.05 & fold_diff > 0 ~ "up",
      adj_p_value < 0.05 & fold_diff < 0 ~ "down",
      TRUE ~ "not_significant"
    )
  )
names_list <- lapply(results$motif, function(idx) motifsToScan[[idx]]@name)
results$name <- unlist(names_list)
results$label <- ifelse(results$name %in% c("Atf3", "Foxf1", "Foxj2", "Foxj3", 
                                            "Foxl2", "Foxn1", "Foxo1", "Foxo3", 
                                            "Jun", "Mafg","Pparg::Rxra"), yes = results$name,no = "")

ggplot(results, aes(x = fold_diff, y = -log10(adj_p_value), color = significance)) +
        geom_point(size = 1.5) +
        scale_color_manual(values = c("up" = "#C22927", 
                                        "down" = "#4164AE", 
                                        "not_significant" = "gray")) +
        ggrepel::geom_text_repel(aes(label = label),
                                size = 4,
                                box.padding = 0.5,
                                max.overlaps = 1000)+
        labs(x = "Ratio Difference", y = "-log10(Adjusted p-value)") +
        theme_bw()


