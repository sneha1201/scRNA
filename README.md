# 🔬 Single-Cell RNA-seq Analysis with Seurat

This pipeline performs clustering, visualization, and annotation of PBMC 10k scRNA-seq data using Seurat and SingleR.

---

## 📦 Required Packages

```r
install.packages("Seurat")
install.packages("tidyverse")
install.packages("patchwork")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("SingleR")
BiocManager::install("celldex")
📁 Dataset
🔗 Download PBMC 10k v3 Raw Feature Matrix (.h5)


1️⃣ Quality Control (QC)
Purpose: Filter out low-quality cells/genes using metrics like mitochondrial/ribosomal gene content and total gene counts, retaining only high-quality cells for analysis.

R Code:

r
# Load data
data <- Read10X(data.dir = "/path/to/data")
rownames(data) <- gsub("_", "-", rownames(data))
df <- CreateSeuratObject(counts = data, project = "sc_project", min.cells = 3, min.features = 200)
rm(data)

# Calculate QC metrics
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
df[["percent.rb"]] <- PercentageFeatureSet(df, pattern = "^RP[SL]")

# Visualize metrics
VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1)

# Scatter plots
FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.rb")

# Filter based on thresholds
df <- subset(df, subset = nCount_RNA < 75000 & nFeature_RNA < 5000 & percent.mt < 15 & percent.rb < 50)

# Re-visualize post-filter
VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1)
Images:

Initial QC metrics: filtered_volcanoplot.png

Violin/Scatter plots: raw_scatter_plot.png, Elbowplot.png

2️⃣ Data Normalization & Highly Variable Genes
Purpose: Normalize gene expression and identify highly variable genes for downstream analyses.

R Code:

r
df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(df), 10)
VariableFeaturePlot(df) + scale_x_continuous(trans = "log1p")
LabelPoints(plot = VariableFeaturePlot(df), points = top10, repel = TRUE)
Images:

Variable gene plot: variable_gene_top10.png

3️⃣ Dimension Reduction: PCA Visualization
Purpose: Reduce data dimensionality with PCA and visualize sample separation and feature contributions.

R Code:

r
df <- ScaleData(df, features = VariableFeatures(object = df))
df <- RunPCA(df, features = VariableFeatures(object = df))
VizDimLoadings(df, dims = 1:6, nfeatures = 15, reduction = "pca")
DimPlot(df, reduction = "pca")
DimHeatmap(df, dims = 1:2, cells = 500, balanced = TRUE)
ElbowPlot(df)
Images:

PCA plot: pca_plot.png, 10_pc_plot

Dimensional heatmap: Dimensional__Heatmap.png

Elbow plot: Elbowplot.png

4️⃣ Clustering & Nonlinear Dimension Reduction (UMAP/t-SNE)
Purpose: Discover clusters and visualize them with UMAP and t-SNE.

R Code:

r
df <- FindNeighbors(df, dims = 1:10)
df <- FindClusters(df, resolution = 0.5)
df <- RunUMAP(df, dims = 1:10)
df <- RunTSNE(df, dims = 1:10)
DimPlot(df, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(df, reduction = "tsne", label = TRUE, repel = TRUE)
Images:

UMAP: umap.png, DIm_plot_groups.png

t-SNE: tsne_plot.png

5️⃣ Marker Gene Identification
Purpose: Identify cluster-specific marker genes.

R Code:

r
cluster1.markers <- FindMarkers(df, ident.1 = 1, min.pct = 0.25)
df.markers <- FindAllMarkers(df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5) %>% filter(p_val_adj < 0.05)
top_markers <- df.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
FeaturePlot(df, features = top_markers$gene[1:4])
Images:

Feature plots: feature_plot_top_markers.png, volcano_plot_2_markers.png

6️⃣ Cell Type Annotation
Purpose: Annotate clusters with reference atlases using SingleR and celldex.

R Code:

r
library(celldex)
ref_human_primary <- celldex::HumanPrimaryCellAtlasData()
ref_blueprint <- celldex::BlueprintEncodeData()
sce <- as.SingleCellExperiment(DietSeurat(df))
results_primary <- SingleR(test = sce, assay.type.test = 1, ref = ref_human_primary, labels = ref_human_primary$label.main)
results_blueprint <- SingleR(test = sce, assay.type.test = 1, ref = ref_blueprint, labels = ref_blueprint$label.main)
results_consensus <- ifelse(!is.na(results_primary$pruned.labels), results_primary$pruned.labels, results_blueprint$pruned.labels)
df@meta.data$celltype_consensus <- results_consensus
png("Cell_Type_Annotation.png", width = 2000, height = 1500, res = 300)
df <- SetIdent(df, value = "celltype_consensus")
DimPlot(df, label = TRUE, repel = TRUE, label.size = 3)
dev.off()
Images:

Cell type annotation: Cell_Type_Annotation.png, Annotate_score_Human_Atlas.png, Annotate_score_Blueprint.png

7️⃣ Subtype Identification (e.g., T/NK Cells)
Purpose: Subcluster specific cell types for finer annotation and marker gene discovery.

R Code:

r
T_df <- subset(df, idents = c("CD4+ T-cells", "CD8+ T-cells", "NK cells"))
T_df.markers <- FindAllMarkers(T_df, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
DimPlot(T_df, label = TRUE, repel = TRUE)
FeaturePlot(T_df, features = c("CD3E", "GNLY", "NKG7"))
Images:

UMAP for T/NK subclusters: T_NK_FineGrained_Annotations.png, feature_plot_top_markers.png

8️⃣ Differential Marker Gene Analysis for Myeloid Cells
Purpose: Find differential marker genes in specific cell subpopulations.

R Code:

r
myeloid <- subset(df, idents = c("Monocyte", "Neutrophils", "Macrophage", "Erythroblast", "Platelets"))
DefaultAssay(myeloid) <- "RNA"
myeloid <- SCTransform(myeloid, verbose = FALSE)
myeloid.markers <- FindAllMarkers(myeloid, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5, assay = "SCT")
top_markers <- myeloid.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(top_markers, file = "myeloid_top10_markers.csv", quote = FALSE, row.names = FALSE)
Images:

Violin plots for markers: filtered_volcanoplot.png

📁 File Structure
text
SingleCellRNAseq_Analysis/
├── README.md
├── Code/
│   ├── Install.R
│   └── analysis_script.R
├── images/
│   ├── umap.png
│   ├── pca_plot.png
│   ├── Elbowplot.png
│   ├── feature_plot_top_markers.png
│   ├── variable_gene_top10.png
│   └── ... (other .png outputs)
└── results/
    └── (tables, marker lists, cell type annotations)
🚀 Usage
Install packages:
Run the installation script as above.

Prepare data:
Place raw count matrices in the specified data directory.

Run analysis:
Follow step-by-step in the main R script, adjusting paths as needed.

Review output:
Visualize results in the /images folder and tabular output in /results.

🔍 Outputs and Interpretation
This pipeline produces:

High-quality filtered cell/gene matrices

QC, variable feature, PCA, UMAP/t-SNE, and clustering visualizations

Cluster marker identification tables and plots

Cell type annotation (atlas comparison and fine-grained reference)

Subtype and differential gene analyses

Publication-ready plots in images/
