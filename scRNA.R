#path
setwd("~/Documents/scRNA/")

#library

library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(celldex)
# Load the PBMCs dataset

PBMCs <- Read10X_h5(filename = 'pbmc_10k_v3_raw_feature_bc_matrix.h5')
str(PBMCs)

# Initialize the Seurat object with the raw (non-normalized data).
#min.cells - features that are expressed in atleast 3 cells
#min.cells - keep cells with atleast 200 features/genes
PBMCs.seurat.obj <- CreateSeuratObject(counts =PBMCs , project = "PBMCs", min.cells = 3, min.features = 200)
str(PBMCs.seurat.obj)

# 1. QC -------
View(PBMCs.seurat.obj@meta.data)

#percentange of Mt  high mitochondrial RNA  means the cell is dying or stressed

PBMCs.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(PBMCs.seurat.obj, pattern = "^MT-")

#Ribosomal genes  can dominate RNA expression because they're highly abundant. 
#percentage of Ribosomal genes
PBMCs.seurat.obj[["percent.rb"]] <- PercentageFeatureSet(PBMCs.seurat.obj, pattern = "^RP[SL]")


# Violin plots for QC metrics to inspect distributions
VlnPlot(PBMCs.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
        ncol = 4, pt.size = 0.1, assay = "RNA", layer = "counts") 
 
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

#potential doublets or stressed cells with unusually high RNA content
#Cells with very high gene counts can be multiplets/doublets
# Filter cells based on observed QC thresholds
PBMCs.seurat.obj <- subset(PBMCs.seurat.obj, subset = nCount_RNA < 75000 & 
               nFeature_RNA < 5000 & 
               percent.mt < 15 & 
               percent.rb < 50)

# Violin plots for QC metrics post-filtering
VlnPlot(PBMCs.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
        ncol = 4, pt.size = 0.1, assay = "RNA", layer = "counts") 

#Normalize the Data

# Log-normalize the data
#Feature Counts per cell are divided by the total counts (UMIs) in that cel then multiplied by scale.factor then log1p of that to lower down the effect of 0 count
PBMCs.seurat.obj <- NormalizeData(PBMCs.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features using the 'vst' method
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

#Scale the Data

# Scale only the variable features to reduce noise and improve computation time
#scales the gene expression values for each gene across all cells
PBMCs.seurat.obj <- ScaleData(PBMCs.seurat.obj, features = VariableFeatures(object = PBMCs.seurat.obj))

#if we want to use all gene
#all.genes <- rownames(PBMCs.seurat.obj)
#PBMCs.seurat.obj <- ScaleData(PBMCs.seurat.obj, features = all.genes)

#Perform PCA

# Run PCA on the scaled variable features
PBMCs.seurat.obj <- RunPCA(PBMCs.seurat.obj, features = VariableFeatures(object = PBMCs.seurat.obj))

#View Top Genes Contributing to PCs

# Display the most important genes contributing to the top 6 principal components
print("Top genes contributing to the first 6 principal components:")
print(PBMCs.seurat.obj[["pca"]], dims = 1:6, nfeatures = 5)

# determine dimensionality of the data
ElbowPlot(PBMCs.seurat.obj)

#Visualize PCA Loadings

# Plot the top 15 genes contributing to each of the first 15 PCs
VizDimLoadings(PBMCs.seurat.obj, dims = 1:15, nfeatures = 15, reduction = "pca") +
  ggtitle("Top Genes in PCA Loadings") 

Visualize the PCA results in a basic PCA plot
DimPlot(PBMCs.seurat.obj, reduction = "pca") + 
  ggtitle("PCA Plot")

#Dimensional heatmap for the first 2 PCs with a subset of 500 cells for clarity
DimHeatmap(PBMCs.seurat.obj, dims = 1:2, cells = 500, balanced = TRUE)

#Clustering the Cells

#Find Neighbors and Identify Clusters

# Identify neighbors using the first 15 principal components
PBMCs.seurat.obj <- FindNeighbors(PBMCs.seurat.obj, dims = 1:15)

# understanding resolution
PBMCs.seurat.obj <- FindClusters(PBMCs.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(PBMCs.seurat.obj@meta.data)

DimPlot(PBMCs.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)

Idents(PBMCs.seurat.obj) <- "RNA_snn_res.0.1"

#Run Dimensional Reduction (UMAP and t-SNE)

# Run UMAP for non-linear dimensional reduction
PBMCs.seurat.obj <- RunUMAP(PBMCs.seurat.obj, dims = 1:15, verbose = FALSE)

# Run t-SNE for an alternative dimensional reduction
PBMCs.seurat.obj <- RunTSNE(PBMCs.seurat.obj, dims = 1:15, verbose = FALSE)

# UMAP plot with cluster labels
DimPlot(PBMCs.seurat.obj, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("UMAP of Clusters")

# t-SNE plot with cluster labels
DimPlot(PBMCs.seurat.obj, reduction = "tsne", label = TRUE, repel = TRUE) + 
  ggtitle("t-SNE of Clusters")

#Finding Cluster-Specific Marker Genes

#Find Markers for Cluster 1

# Identify markers for cluster 1 with a minimum expression percentage of 25%
cluster1.markers <- FindMarkers(PBMCs.seurat.obj, ident.1 = 1, min.pct = 0.25)

# Display the top 5 markers for cluster 1
head(cluster1.markers, n = 5)

#Visualize Top Markers for Cluster 1

# Violin plot for the top 2 markers of cluster 1
VlnPlot(PBMCs.seurat.obj, features = c(row.names(cluster1.markers)[1], row.names(cluster1.markers)[2]))

#Find Markers for All Clusters

# Identify markers for all clusters with a log-fold change threshold of 0.5
PBMCs.seurat.obj.marker <- FindAllMarkers(PBMCs.seurat.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5) %>%
  filter(p_val_adj < 0.05)  # Optional: Filter markers by adjusted p-value

#Identify Top Markers for Each Cluster

# Get the top 3 marker genes for each cluster based on average log-fold change
top_markers <- PBMCs.seurat.obj.marker %>% 
  group_by(cluster) %>% 
  top_n(n = 3, wt = avg_log2FC)

# Plotting Marker Genes

# Identify the top marker gene for each cluster
top_markers_1 <- PBMCs.seurat.obj.marker %>% 
  group_by(cluster) %>% 
  top_n(n = 1, wt = avg_log2FC)

# Plot the top marker genes for the first 12 clusters in groups of 4
FeaturePlot(PBMCs.seurat.obj, features = top_markers_1$gene[1:4])

# Cell Type Annotation

# Load Human Primary Cell Atlas reference
ref_human_primary <- celldex::HumanPrimaryCellAtlasData()

# Load Blueprint/ENCODE reference for additional coverage of immune cells
ref_blueprint <- celldex::BlueprintEncodeData()

# Inspect reference datasets
print(ref_human_primary)
print(ref_blueprint)

#Convert Seurat Object to SingleCellExperiment

# Convert Seurat object to SingleCellExperiment using DietSeurat to save memory
sce <- as.SingleCellExperiment(DietSeurat(PBMCs.seurat.obj))

# Run SingleR with Human Primary Cell Atlas reference
results_primary <- SingleR(test = sce, assay.type.test = 1, ref = ref_human_primary, labels = ref_human_primary$label.main)

# Run SingleR with Blueprint/ENCODE reference for comparison
results_blueprint <- SingleR(test = sce, assay.type.test = 1, ref = ref_blueprint, labels = ref_blueprint$label.main)

# Create a consensus annotation by comparing the two results
results_consensus <- ifelse(!is.na(results_primary$pruned.labels), 
                            results_primary$pruned.labels, 
                            results_blueprint$pruned.labels)

# Add annotations to the Seurat object metadata
PBMCs.seurat.obj@meta.data$celltype_primary <- results_primary$pruned.labels
PBMCs.seurat.obj@meta.data$celltype_blueprint <- results_blueprint$pruned.labels
PBMCs.seurat.obj@meta.data$celltype_consensus <- results_consensus

# Plot annotation scores to visualize confidence levels
plotScoreHeatmap(results_primary, main = "Annotation Scores - Human Primary Cell Atlas")

plotScoreHeatmap(results_blueprint, main = "Annotation Scores - Blueprint/ENCODE")

# Set Seurat object identities to consensus cell types and visualize with UMAP
PBMCs.seurat.obj <- SetIdent(PBMCs.seurat.obj, value = "celltype_consensus")
DimPlot(PBMCs.seurat.obj, label = TRUE, repel = TRUE, label.size = 3) + NoLegend() + 
  ggtitle("Consensus Cell Type Annotations")

# Display the number of cells in each annotation
table(PBMCs.seurat.obj@meta.data$celltype_consensus)




