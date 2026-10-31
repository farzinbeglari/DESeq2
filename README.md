# DESeq2
Differential expression analysis (GSE249584) COPD vs. Control
library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyverse)

# Data Import
data_directory <- getwd()
samples <- list.dirs(path = data_directory,
                     full.names = F,
                     recursive = F)
Data.list <- list()
for (sample in samples) {
  sample_directory <- file.path(data_directory, sample)
  
  counts <- ReadMtx(mtx = file.path(sample_directory, 'matrix.mtx.gz'),
                    cells = file.path(sample_directory, 'barcodes.tsv.gz'),
                    features = file.path(sample_directory, 'features.tsv.gz'),
                    feature.column = 2)
  
  Data.list[[sample]] <- CreateSeuratObject(counts = counts,
                                            project = sample)
}

# Merging Seurat Objects From 15 Samples
merged.data <- merge(x = Data.list[[1]],
                     y = Data.list[-1],
                     add.cell.ids = names(Data.list))
merged.data$sample <- merged.data$orig.ident
view(merged.data@meta.data)

Data <- AggregateExpression(merged.data,
assays = 'RNA',
return.seurat = F,
group.by = 'sample')

Data <- Data$RNA
Data <- as.data.frame(Data)

# Keep the genes with row sums counts >= 10
dds <- Data[rowSums(Data) >= 10,]

# colData annotation
colData <- as.data.frame(colnames(dds))
colData <- colData %>%
  mutate(condition = ifelse(grepl('Control', `colnames(dds)`),'Control', 'COPD')) %>%
  column_to_rownames(var = 'colnames(dds)')

colData$condition <- factor(colData$condition)

# Differential expression analysis (DEseq2)
dds <- DESeqDataSetFromMatrix(dds,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)

# Extract results: Genes with P adj value < 0.05 (alpha = 0.05) were considered as significant change
results <- results(dds, alpha = 0.05)
view(results)
plotMA(results)
