# scRNA
single cell Analysis

# 🔬 Single-Cell RNA-seq Analysis using Seurat (PBMC 10k Dataset)

This project is a complete single-cell RNA-seq (scRNA-seq) analysis workflow using **Seurat**, with cell type annotation using **SingleR** and reference datasets from **celldex**. We process a PBMC dataset from 10X Genomics, perform quality control, clustering, dimensionality reduction, and cell-type identification.

---

## 📁 Dataset

We use the [PBMC 10k v3 dataset from 10x Genomics](https://www.10xgenomics.com/resources/datasets/10-k-pbmcs-from-a-healthy-donor-single-indexed-3-v3-1-standard-3-0-0).

---

## 🧪 Step 1: Load Data & Create Seurat Object

```r
setwd("~/Documents/scRNA/")
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(celldex)

PBMCs <- Read10X_h5(filename = 'pbmc_10k_v3_raw_feature_bc_matrix.h5')
PBMCs.seurat.obj <- CreateSeuratObject(counts =PBMCs, project = "PBMCs", min.cells = 3, min.features = 200)

Step 2: Quality Control (QC)
We calculate mitochondrial and ribosomal gene percentages to identify poor-quality or stressed cells.

r
Copy
Edit
PBMCs.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(PBMCs.seurat.obj, pattern = "^MT-")
PBMCs.seurat.obj[["percent.rb"]] <- PercentageFeatureSet(PBMCs.seurat.obj, pattern = "^RP[SL]")

QC Violin Plots

VlnPlot(PBMCs.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
        ncol = 4, pt.size = 0.1, assay = "RNA", layer = "counts") 

![image](https://github.com/user-attachments/assets/2b3fa63c-051a-437f-b8db-6eb30aff5609)

Feature Scatter Plots

# Scatter plots for QC metric relationships
#X-axis: nCount_RNA = total number of UMIs (transcripts) per cell
#Y-axis: percent.mt = % of mitochondrial gene expression

plot1 <- FeatureScatter(PBMCs.seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  ggtitle("RNA Count vs Mitochondrial Percentage") + theme_minimal()+
  theme_minimal(base_size = 5)

#X-axis: nCount_RNA
#Y-axis: nFeature_RNA = number of unique genes detected

plot2 <- FeatureScatter(PBMCs.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  ggtitle("RNA Count vs Feature Count") + theme_minimal() +
  theme_minimal(base_size = 5)

#X-axis: nCount_RNA
#Y-axis: percent.rb = % expression of riboso
plot3 <- FeatureScatter(PBMCs.seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.rb") + 
  ggtitle("RNA Count vs Ribosomal Percentage") + theme_minimal()+
  theme_minimal(base_size = 5)


plot1 + plot2 + plot3

![image](https://github.com/user-attachments/assets/df4b12ee-d65f-4d8a-91ed-6163ec1f8aea)

Step 3: Filter Cells Based on QC
r
Copy
Edit
PBMCs.seurat.obj <- subset(PBMCs.seurat.obj, subset = 
  nCount_RNA < 75000 & 
  nFeature_RNA < 5000 & 
  percent.mt < 15 & 
  percent.rb < 50)

# Violin plots for QC metrics post-filtering
VlnPlot(PBMCs.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
        ncol = 4, pt.size = 0.1, assay = "RNA", layer = "counts") 

![image](https://github.com/user-attachments/assets/aa4570c9-ff2f-4f49-bf85-a5ec6b18110e)

Step 4: Normalization & Feature Selection
r
Copy
Edit
PBMCs.seurat.obj <- NormalizeData(PBMCs.seurat.obj)
PBMCs.seurat.obj <- FindVariableFeatures(PBMCs.seurat.obj, selection.method = "vst", nfeatures = 3000)

# Extract the top 10 most highly variable genes
top10 <- head(VariableFeatures(PBMCs.seurat.obj), 10)
print(top10)

#Visualize Highly Variable Features

# Plot variable features with log transformation to handle zeros safely
plot1 <- VariableFeaturePlot(PBMCs.seurat.obj) + 
  scale_x_continuous(trans = "log1p") + 
  ggtitle("Highly Variable Features") + 
  theme_minimal()

# Label the top 10 variable genes
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0, max.overlaps = 10)

# Display the plot
plot2

![image](https://github.com/user-attachments/assets/1c0d8c8b-462b-4456-9c69-f5d207669507)



