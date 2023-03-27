### Reference: https://satijalab.org/seurat/articles/integration_introduction.html

install.packages("mailR")
install.packages("plyr")
install.packages('xlsx')
# install.packages("dplyr")
# install.packages("Seurat")
# install.packages("ggplot2")
# install.packages("cowplot")
# install.packages("patchwork")
# install packages for trajectory analysis
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
BiocManager::install("monocle")
install.packages("circlize")
install.packages("LSD")
install.packages("RColorBrewer")
BiocManager::install("Nebulosa")
# for GRNanalysis
BiocManager::install(c("AUCell", "RcisTarget"),ask = F,update = F) 
BiocManager::install(c("GENIE3"),ask = F,update = F)  # Optional. Can be replaced by GRNBoost
BiocManager::install(c("zoo", "mixtools", "rbokeh"),ask = F,update = F) 
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"),ask = F,update = F)
BiocManager::install(c("doMC", "doRNG"),ask = F,update = F)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")
install.packages('psych')
install.packages('qgraph')
install.packages('stringi')
install.packages('ggsci')
install.packages('igraph')
install.packages('glpk')
#先配置GLPK的环境，告诉R GLPK要到哪里去找
Sys.setenv(CPPFLAGS="-DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /lustre/home/acct-medzt/medzt/.conda/envs/cpdb/include")
Sys.setenv(LDFLAGS="-Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,/lustre/home/acct-medzt/medzt/.conda/envs/cpdb/lib -Wl,-rpath-link,/lustre/home/acct-medzt/medzt/.conda/envs/cpdb/lib -L/lustre/home/acct-medzt/medzt/.conda/envs/cpdb/lib")
remotes::install_github("igraph/rigraph@master")
if(!require(devtools)) install.packages("devtools");
devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
install.packages('randomcoloR')
install.packages('V8')
remove.packages('igraph')
devtools::install_version('igraph', version = "1.1.2",repos = "http://cran.us.r-project.org")
install.packages("rJava")
BiocManager::install("GeneOverlap")


options(connectionObserver = NULL)
library(dplyr)
library(plyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot())
library(ComplexHeatmap)
library(circlize)
library(LSD)
library(RColorBrewer)
library(Nebulosa)
library(xlsx)

library(psych)
library(qgraph)
library(igraph)
library(stringi)
library(ggsci)

library(monocle)
library(monocle3)


# load merged scRNA-seq data (a combined object encompassing several seuratObjects, every of which was quality controlled respectively ahead of combination)
load("scRAW data/sc.kNC_kJID_hNC.anchors2000ForIntegration.Rdata")

# split the dataset into a list of multiple seurat objects
sc.kNC_kJID_hNC.list <- SplitObject(sc.kNC_kJID_hNC.merge, split.by = "orig.ident")

# ______normalize and identify variable features for each dataset independently (*TaT: please define the "nfeatures" parameter to ensure that the gene of interest exists in the top "nfeatures" HVGs of every single sample/seuratObject via trial and error {NormalizeData() & FindVariableFeatures()} before combination)
sc.kNC_kJID_hNC.list_3000 <- lapply(X = sc.kNC_kJID_hNC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# ______select features that are repeatedly variable across datasets for integration (*TaT: please define the "nfeatures" parameter via trial and error to ensure that the gene of interest exists in the top "nfeatures" results)
features_2000 <- SelectIntegrationFeatures(object.list = sc.kNC_kJID_hNC.list_3000, nfeatures = 2000)
write.table(features_2000, file = "2000variableGenes_for_integration.xlsx", quote = F, sep = "\t", row.names = F, col.names = F)
# this step takes the longest time to execute
sc.anchors_2000 <- FindIntegrationAnchors(object.list = sc.kNC_kJID_hNC.list_3000, anchor.features = features_2000)
save(sc.anchors_2000, file = "sc.kNC_kJID_hNC.anchors2000ForIntegration.Rdata")
load("scRAW data/sc.kNC_kJID_hNC.anchors2000ForIntegration.Rdata")

# this command creates an 'integrated' data assay
sc.kNC_kJID_hNC.integrated_2000 <- IntegrateData(anchorset = sc.anchors_2000)
###

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000) <- "integrated"
save(sc.kNC_kJID_hNC.integrated_2000, file = "R/sc.kNC_kJID_hNC_3000HVG_integrated_2000.Rdata")

# Run the standard workflow for visualization and clustering
sc.kNC_kJID_hNC.integrated_2000 <- ScaleData(sc.kNC_kJID_hNC.integrated_2000, verbose = FALSE)
sc.kNC_kJID_hNC.integrated_2000 <- RunPCA(sc.kNC_kJID_hNC.integrated_2000, npcs = 30, verbose = FALSE)
sc.kNC_kJID_hNC.integrated_2000 <- RunUMAP(sc.kNC_kJID_hNC.integrated_2000, reduction = "pca", dims = 1:30)
Idents(sc.kNC_kJID_hNC.integrated_2000) <- "orig.ident"
# review integrated scRNA-seq cells in UMAP
DimPlot(sc.kNC_kJID_hNC.integrated_2000, reduction = "umap", shuffle = T)
# find neighbors for clustering
sc.kNC_kJID_hNC.integrated_2000 <- FindNeighbors(sc.kNC_kJID_hNC.integrated_2000, reduction = "pca", dims = 1:30)
save(sc.kNC_kJID_hNC.integrated_2000, file = "R/sc.kNC_kJID_hNC.integrated2000_NeighbourFound.Rdata")

load("R/sc.kNC_kJID_hNC.integrated2000_NeighbourFound.Rdata")


# ______Explore 10 clusters using 0.1 resolution (*TaT: please define the "resolution" parameter, increase of which leads to a greater number of clusters. Setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.)
sc.kNC_kJID_hNC.integrated_2000_0.1 <- FindClusters(sc.kNC_kJID_hNC.integrated_2000, resolution = 0.1)
DimPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, reduction = "umap", label = T, repel = T)
DimPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, reduction = "umap", label = T, repel = T, split.by = "zt.scarType")
save(sc.kNC_kJID_hNC.integrated_2000_0.1, file = "R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")

# ______Explore 23 clusters using 0.6 resolution (*TaT: please define the "resolution" parameter, increase of which leads to a greater number of clusters. Setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.)
sc.kNC_kJID_hNC.integrated_2000_0.6 <- FindClusters(sc.kNC_kJID_hNC.integrated_2000, resolution = 0.6)
DimPlot(sc.kNC_kJID_hNC.integrated_2000_0.6, reduction = "umap", label = T, repel = T)
save(sc.kNC_kJID_hNC.integrated_2000_0.6, file = "R/sc.kNC_kJID_hNC.integrated2000_0.6resolution_23clusters.Rdata")


#################### expression of well-established cell type markers #############################
load("R/sc.kNC_kJID_hNC.integrated2000_0.6resolution_23clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000_0.6) <- "RNA"

## cell types according to markers in three articles
# fibroblasts
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("COL1A1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("COL1A2"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("COL3A1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("PDGFRA"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("DCN"))
# SMC
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("TAGLN"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("ACTA2"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("TPM2"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("MYH11"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("CNN1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("RGS5"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("PDGFRB")) #muralCells
# endothelial cells
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("SELE"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("TM4SF1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("PECAM1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("ENG"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("CDH5"))
# keratinocyte
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("KRT14"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("KRT1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("KRT10"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("KRT5"))
# lymphatic endothelial cells
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("CCL21"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("LYVE1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("PROX1"))
# melanocyte
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("TYRP1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("PMEL"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("MLANA"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("DCT"))
# immune cells
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("LYZ"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("HLA-DRA"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("PTPRC"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("CXCR4"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("CD3D")) #Tcells
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("AIF1")) #Macrophages
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("ITGAX")) #dentriticCells
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("CD1A")) #dentriticCells
# sweat gland cells
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("SCGB1B2P"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("SCGB1D2"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("AQP5"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("MUCL1"))
# neural cells
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("NRXN1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("SCN7A"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("GLDN")) # Schwann cells
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("TJP1")) #widespread

#################################################


#################### identify DEGs of each cell type #############################
## identify DEGs of each one of 10 clusters
load("R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "RNA"
options(future.globals.maxSize= 3210612736)

# cluster 0 in 10 clusters
c0.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = 0, verbose = FALSE, logfc.threshold = 0.1)
write.table(c0.markers, file = "R/sc.KNC_kJID_hNC/c0.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 1 in 10 clusters
c1.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = 1, verbose = FALSE, logfc.threshold = 0.1)
write.table(c1.markers, file = "R/sc.KNC_kJID_hNC/c1.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 2 in 10 clusters
c2.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = 2, verbose = FALSE, logfc.threshold = 0.1)
write.table(c2.markers, file = "R/sc.KNC_kJID_hNC/c2.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 3 in 10 clusters
c3.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = 3, verbose = FALSE, logfc.threshold = 0.1)
write.table(c3.markers, file = "R/sc.KNC_kJID_hNC/c3.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 4 in 10 clusters
c4.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = 4, verbose = FALSE, logfc.threshold = 0.1)
write.table(c4.markers, file = "R/sc.KNC_kJID_hNC/c4.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 5 in 10 clusters
c5.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = 5, verbose = FALSE, logfc.threshold = 0.1)
write.table(c5.markers, file = "R/sc.KNC_kJID_hNC/c5.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 6 in 10 clusters
c6.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = 6, verbose = FALSE, logfc.threshold = 0.1)
write.table(c6.markers, file = "R/sc.KNC_kJID_hNC/c6.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 7 in 10 clusters
c7.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = 7, verbose = FALSE, logfc.threshold = 0.1)
write.table(c7.markers, file = "R/sc.KNC_kJID_hNC/c7.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 8 in 10 clusters
c8.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = 8, verbose = FALSE, logfc.threshold = 0.1)
write.table(c8.markers, file = "R/sc.KNC_kJID_hNC/c8.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 9 in 10 clusters
c9.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = 9, verbose = FALSE, logfc.threshold = 0.1)
write.table(c9.markers, file = "R/sc.KNC_kJID_hNC/c9.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

#################################################


#################### identify DEGs of 23 clusters #############################
## identify DEGs of each one in 23 clusters
options(future.globals.maxSize= 3210612736)
load("R/sc.kNC_kJID_hNC.integrated2000_0.6resolution_23clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000_0.6) <- "RNA"
# cluster 0 in 23 clusters
c0.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 0, verbose = FALSE, logfc.threshold = 0.1)
write.table(c0.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c0.markers", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 1 in 23 clusters
c1.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 1, verbose = FALSE, logfc.threshold = 0.1)
write.table(c1.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c1.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 2 in 23 clusters
c2.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 2, verbose = FALSE, logfc.threshold = 0.1)
write.table(c2.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c2.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 3 in 23 clusters
c3.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 3, verbose = FALSE, logfc.threshold = 0.1)
write.table(c3.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c3.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 4 in 23 clusters
c4.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 4, verbose = FALSE, logfc.threshold = 0.1)
write.table(c4.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c4.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 5 in 23 clusters
c5.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 5, verbose = FALSE, logfc.threshold = 0.1)
write.table(c5.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c5.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 6 in 23 clusters
c6.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 6, verbose = FALSE, logfc.threshold = 0.1)
write.table(c6.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c6.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 7 in 23 clusters
c7.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 7, verbose = FALSE, logfc.threshold = 0.1)
write.table(c7.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c7.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 8 in 23 clusters
c8.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 8, verbose = FALSE, logfc.threshold = 0.1)
write.table(c8.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c8.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 9 in 23 clusters
c9.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 9, verbose = FALSE, logfc.threshold = 0.1)
write.table(c9.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c9.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 10 in 23 clusters
c10.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 10, verbose = FALSE, logfc.threshold = 0.1)
write.table(c10.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c10.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 11 in 23 clusters
c11.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 11, verbose = FALSE, logfc.threshold = 0.1)
write.table(c11.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c11.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 12 in 23 clusters
c12.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 12, verbose = FALSE, logfc.threshold = 0.1)
write.table(c12.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c12.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 13 in 23 clusters
c13.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 13, verbose = FALSE, logfc.threshold = 0.1)
write.table(c13.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c13.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 14 in 23 clusters
c14.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 14, verbose = FALSE, logfc.threshold = 0.1)
write.table(c14.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c14.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 15 in 23 clusters
c15.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 15, verbose = FALSE, logfc.threshold = 0.1)
write.table(c15.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c15.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 16 in 23 clusters
c16.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 16, verbose = FALSE, logfc.threshold = 0.1)
write.table(c16.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c16.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 17 in 23 clusters
c17.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 17, verbose = FALSE, logfc.threshold = 0.1)
write.table(c17.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c17.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 18 in 23 clusters
c18.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 18, verbose = FALSE, logfc.threshold = 0.1)
write.table(c18.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c18.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 19 in 23 clusters
c19.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 19, verbose = FALSE, logfc.threshold = 0.1)
write.table(c19.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c19.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 20 in 23 clusters
c20.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 20, verbose = FALSE, logfc.threshold = 0.1)
write.table(c20.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c20.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 21 in 23 clusters
c21.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 21, verbose = FALSE, logfc.threshold = 0.1)
write.table(c21.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c21.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 22 in 23 clusters
c22.markers <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.6, ident.1 = 22, verbose = FALSE, logfc.threshold = 0.1)
write.table(c22.markers, file = "R/sc.KNC_kJID_hNC/23_clusters_marker_genes/c22.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

#################################################


#################### highlight a cluster in UMAP #############################
# highlight a specific cluster in UMAP
DimPlot(sc.kNC_kJID_hNC.integrated_2000_0.6, reduction = "umap", label = T, repel = T, cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.integrated_2000_0.6, idents = "0"), cols.highlight = "red")
DimPlot(sc.kNC_kJID_hNC.integrated_2000_0.6, reduction = "umap", label = T, repel = T, cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.integrated_2000_0.6, idents = "19"), cols.highlight = "red")
DimPlot(sc.kNC_kJID_hNC.integrated_2000_0.6, reduction = "umap", label = T, repel = T, split.by = "zt.scarType",
        cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.integrated_2000_0.6, idents = "19"), cols.highlight = "red")
DimPlot(sc.kNC_kJID_hNC.integrated_2000_0.6, reduction = "umap", label = T, repel = T, split.by = "zt.dataSource",
        cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.integrated_2000_0.6, idents = "19"), cols.highlight = "red")
DimPlot(sc.kNC_kJID_hNC.integrated_2000_0.6, reduction = "umap", label = T, repel = T, split.by = "orig.ident",
        cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.integrated_2000_0.6, idents = "19"), cols.highlight = "red")

# highlight a specific FBsubcluster in UMAP
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T)
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T, cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.FB_0.16, idents = c("0")), cols.highlight = "#e87d72")
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T, split.by = "zt.scarType",
        cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.FB_0.16, idents = c("0")), cols.highlight = "#e87d72")
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T, cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.FB_0.16, idents = c("1", "2")), cols.highlight = c("#53b64c", "#b49f33"))
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T, split.by = "zt.scarType",
        cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.FB_0.16, idents = c("1", "2")), cols.highlight = c("#53b64c", "#b49f33"))
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T, cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.FB_0.16, idents = c("3")), cols.highlight = "#54bcc2")
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T, split.by = "zt.scarType",
        cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.FB_0.16, idents = c("3")), cols.highlight = "#54bcc2")
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T, cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.FB_0.16, idents = c("4", "5")), cols.highlight = c("#e470dd", "#6c9cf8"))
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T, split.by = "zt.scarType",
        cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.FB_0.16, idents = c("4", "5")), cols.highlight = c("#e470dd", "#6c9cf8"))
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T, cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.FB_0.16, idents = c("2")), cols.highlight = c("#53b64c"))
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T, split.by = "zt.scarType",
        cells.highlight = CellsByIdentities(sc.kNC_kJID_hNC.FB_0.16, idents = c("2")), cols.highlight = c("#53b64c"))


#################################################


#################### highlight a cluster in LSD::Heatscatter #############################
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# backup
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
# heatscatter background
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = "SDC3", cols = c("lightgrey", "lightgrey")) + xlim(-8, 8) + ylim(-8, 8)

## FB.HTS/Keloid/Normal
# FB
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "FB")
# heatscatter of FB.HTS
Idents(sc.kNC_kJID_hNC.FB_0.16_test1) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, idents = "HTS")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.HTS@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.HTS@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "HTS")
# heatscatter of FB.Keloid
Idents(sc.kNC_kJID_hNC.FB_0.16_test1) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, idents = "Keloid")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.Keloid@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.Keloid@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Keloid")
# heatscatter of FB.Normal
Idents(sc.kNC_kJID_hNC.FB_0.16_test1) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.Normal <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, idents = "Normal")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.Normal@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.Normal@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Normal")





## FBc0
# heatscatter of FBc0
sc.kNC_kJID_hNC.FB_0.16_test1.c0 <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, idents = "0")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c0@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c0@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "FB")
# heatscatter of activated FBc0.HTS
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c0) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c0.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c0, idents = "HTS")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c0.HTS@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c0.HTS@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "HTS")
# heatscatter of activated FBc0.Keloid
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c0) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c0.Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c0, idents = "Keloid")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c0.Keloid@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c0.Keloid@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Keloid")
# heatscatter of activated FBc0.Normal
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c0) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c0.Normal <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c0, idents = "Normal")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c0.Normal@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c0.Normal@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Normal")

## FBc12
# heatscatter of FBc12
sc.kNC_kJID_hNC.FB_0.16_test1.c12 <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, idents = c("1", "2"))
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c12@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c12@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "FB")
# heatscatter of activated FBc12.HTS
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c12) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c12.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c12, idents = "HTS")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c12.HTS@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c12.HTS@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "HTS")
# heatscatter of activated FBc12.Keloid
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c12) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c12.Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c12, idents = "Keloid")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c12.Keloid@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c12.Keloid@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Keloid")
# heatscatter of activated FBc12.Normal
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c12) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c12.Normal <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c12, idents = "Normal")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c12.Normal@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c12.Normal@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Normal")

## FBc3
# heatscatter of FBc3
sc.kNC_kJID_hNC.FB_0.16_test1.c3 <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, idents = c("3"))
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c3@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c3@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "FB")
# heatscatter of activated FBc3.HTS
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c3) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c3.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c3, idents = "HTS")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c3.HTS@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c3.HTS@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "HTS")
# heatscatter of activated FBc3.Keloid
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c3) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c3.Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c3, idents = "Keloid")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c3.Keloid@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c3.Keloid@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Keloid")
# heatscatter of activated FBc3.Normal
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c3) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c3.Normal <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c3, idents = "Normal")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c3.Normal@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c3.Normal@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Normal")

## FBc45
# heatscatter of FBc45
sc.kNC_kJID_hNC.FB_0.16_test1.c45 <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, idents = c("4", "5"))
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c45@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c45@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "FB")
# heatscatter of activated FBc45.HTS
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c45) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c45.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c45, idents = "HTS")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c45.HTS@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c45.HTS@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "HTS")
# heatscatter of activated FBc45.Keloid
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c45) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c45.Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c45, idents = "Keloid")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c45.Keloid@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c45.Keloid@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Keloid")
# heatscatter of activated FBc45.Normal
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.c45) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.c45.Normal <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.c45, idents = "Normal")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.c45.Normal@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.c45.Normal@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Normal")


#################################################


#################### explore IL11RA #############################
### deep dive into IL11RA
## review IL11RA in all cell types
load("R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "RNA"
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("IL11RA"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("IL11RA"), split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("IL11RA"), combine = FALSE, group.by = "zt.scarType")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("IL11RA"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
# IL11RA in HTS
sc.kNC_kJID_hNC.integrated_2000_0.1_HTS <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "HTS")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1_HTS) <- "zt.cellType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1_HTS, features = c("IL11RA"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1_HTS, features = c("IL11RA"), combine = FALSE)
# IL11RA in Keloid
sc.kNC_kJID_hNC.integrated_2000_0.1_Keloid <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1_Keloid) <- "zt.cellType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1_Keloid, features = c("IL11RA"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1_Keloid, features = c("IL11RA"), combine = FALSE)
# IL11RA in NormalSkin
sc.kNC_kJID_hNC.integrated_2000_0.1_Normal <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Normal")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1_Normal) <- "zt.cellType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1_Normal, features = c("IL11RA"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1_Normal, features = c("IL11RA"), combine = FALSE)

## review IL11RA in FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.168resolution_8clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.168) <- "RNA"
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("IL11RA"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("IL11RA"), split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.168, features = c("IL11RA"), combine = FALSE, group.by = "zt.scarType")
Idents(sc.kNC_kJID_hNC.FB_0.168) <- "integrated_snn_res.0.168"
DotPlot(sc.kNC_kJID_hNC.FB_0.168, features = c("IL11RA"), cols = c("blue", "red", "blue", "red", "blue", "red"), dot.scale = 8) + RotatedAxis()
# IL11RA in HTS FB
Idents(sc.kNC_kJID_hNC.FB_0.168) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.168_HTS <- subset(sc.kNC_kJID_hNC.FB_0.168, idents = "HTS")
Idents(sc.kNC_kJID_hNC.FB_0.168_HTS) <- "integrated_snn_res.0.168"
DotPlot(sc.kNC_kJID_hNC.FB_0.168_HTS, features = c("IL11RA"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
VlnPlot(sc.kNC_kJID_hNC.FB_0.168_HTS, features = c("IL11RA"), combine = FALSE)
# IL11RA in Keloid FB
sc.kNC_kJID_hNC.FB_0.168_Keloid <- subset(sc.kNC_kJID_hNC.FB_0.168, idents = "Keloid")
Idents(sc.kNC_kJID_hNC.FB_0.168_Keloid) <- "integrated_snn_res.0.168"
DotPlot(sc.kNC_kJID_hNC.FB_0.168_Keloid, features = c("IL11RA"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
VlnPlot(sc.kNC_kJID_hNC.FB_0.168_Keloid, features = c("IL11RA"), combine = FALSE)
# IL11RA in Normal FB
sc.kNC_kJID_hNC.FB_0.168_Normal <- subset(sc.kNC_kJID_hNC.FB_0.168, idents = "Normal")
Idents(sc.kNC_kJID_hNC.FB_0.168_Normal) <- "integrated_snn_res.0.168"
DotPlot(sc.kNC_kJID_hNC.FB_0.168_Normal, features = c("IL11RA"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
VlnPlot(sc.kNC_kJID_hNC.FB_0.168_Normal, features = c("IL11RA"), combine = FALSE)


#################################################


#################### explore a particular gene: CD109 #############################
### deep dive into CD109
## review CD109 in all cell types
load("R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "RNA"
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("CD109"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("CD109"), split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("CD109"), combine = T, group.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("CD109"), combine = T, group.by = "zt.scarType", pt.size = 0)
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("CD109"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("CD109"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()+ scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType_scarType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("CD109"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("CD109"), cols = c("blue", "grey", "red"), dot.scale = 8, split.by = "zt.scarType") + RotatedAxis()

## deep dive into subseted FB
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
sc.kNC_kJID_hNC.integrated_2000_0.1.FB <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "FB")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1.FB) <- "zt.cellType_scarType"
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("CD109"), split.by = "zt.cellType_scarType")
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("CD109"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("CD109"), combine = T)
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("CD109"), combine = T, pt.size = 0)
# DEGs between FB_Keloid and FB_Normal
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1.FB) <- "zt.cellType_scarType"
DEG.FB.KeloidvsNormal <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, ident.1 = "FB_Keloid", ident.2 = "FB_Normal", verbose = FALSE, logfc.threshold = 0.01)
write.table(DEG.FB.KeloidvsNormal, file = "R/sc.KNC_kJID_hNC/CD109/DEG.FB.KeloidvsNormal.markers_FC0.01.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

## deep dive into subclustered FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# expression distribution
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("CD109"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("CD109"), split.by = "zt.scarType")
# expression in different scarType
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("CD109"), combine = T, group.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("CD109"), combine = T, group.by = "zt.scarType", pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("CD109"), dot.scale = 8, group.by = "zt.scarType") + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
# expression in different subclusters
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("CD109"), combine = T)
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("CD109"), combine = T, pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("CD109"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("CD109"), combine = T, split.by = "zt.scarType")
#################################################
#################### explore smc and pericyte markers #############################
## 
load("R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "RNA"
# expression distribution
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("RGS5"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("TAGLN"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("ACTA2"))

#################### explore a particular gene: DPP4 #############################
### deep dive into DPP4
## review DPP4 in all cell types
load("R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "RNA"
# expression distribution
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("DPP4"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("DPP4"), split.by = "zt.scarType")
# expression in different scarType
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("DPP4"), combine = T, group.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("DPP4"), combine = T, group.by = "zt.scarType", pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("DPP4"), cols = c("blue", "red"), dot.scale = 8, group.by = "zt.scarType") + RotatedAxis()
# expression in different cellType
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("DPP4"), combine = T)
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("DPP4"), combine = T, pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("DPP4"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("DPP4"), combine = T, split.by = "zt.scarType")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType_scarType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("DPP4"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')

## deep dive into subseted FB
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
sc.kNC_kJID_hNC.integrated_2000_0.1.FB <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "FB")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1.FB) <- "zt.cellType_scarType"
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("DPP4"), split.by = "zt.cellType_scarType")
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("DPP4"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("DPP4"), combine = T)
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("DPP4"), combine = T, pt.size = 0)

## deep dive into subclustered FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# expression distribution
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("DPP4"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("DPP4"), split.by = "zt.scarType")
# expression in different scarType
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("DPP4"), combine = T, group.by = "zt.cellType_scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("DPP4"), combine = T, group.by = "zt.cellType_scarType", pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("DPP4"), dot.scale = 8, group.by = "zt.cellType_scarType") + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
# expression in different subclusters
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("DPP4"), combine = T)
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("DPP4"), combine = T, pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("DPP4"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("DPP4"), combine = T, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("DPP4"), combine = T, split.by = "zt.scarType", pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("DPP4"), dot.scale = 8, group.by = "zt.FB6cluster_scarType") + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
#################################################


#################### expression of TOP3 DEGs of each cell type #############################
## cell type markers defined by TOP3 DEGs 
# cluster 0 in 10 clusters
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("COL1A1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("DCN"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("CFD"))
# cluster 1 in 10 clusters
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("SELE"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("TM4SF1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("ACKR1"))
# cluster 2 in 10 clusters
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("KRT1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("KRT10"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("DMKN"))
# cluster 3 in 10 clusters
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("TAGLN"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("RGS5"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("ACTA2"))
# cluster 4 in 10 clusters
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("KRT15"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("KRT14"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("S100A2"))
# cluster 5 in 10 clusters
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("CXCR4"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("PTPRC"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("SRGN"))
# cluster 6 in 10 clusters
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("CXCL8"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("HLA-DRA"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("LYZ"))
# cluster 7 in 10 clusters
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("CCL21"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("TFF3"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("MMRN1"))
# cluster 8 in 10 clusters
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("DCT"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("TYRP1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("PMEL"))
# cluster 9 in 10 clusters
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("DCD"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("MUCL1"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = c("SCGB2A2"))

## heatmap of cell markers defined by TOP3 DEGs
Idents(sc.kNC_kJID_hNC.integrated_2000_0.6) <- factor(Idents(sc.kNC_kJID_hNC.integrated_2000_0.6), levels = c("0","3","7","11","18","1","4","8","19","2","12","17","6","13","5","14","9","10","15","16","20","21","22"))
markers.to.plot <- c("COL1A1", "DCN", "CFD", "SELE", "TM4SF1", "ACKR1", "KRT1", "KRT10", "DMKN", "KRT15", "KRT14", "S100A2", "TAGLN", "RGS5", "ACTA2", "CXCR4", "PTPRC", "SRGN", "CXCL8", "HLA-DRA", "LYZ", "CCL21", "TFF3", "MMRN1", "DCT", "TYRP1", "PMEL", "DCD", "MUCL1", "SCGB2A2")
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = markers.to.plot,cols = c("lightgrey", "darkred"), dot.scale = 7) + RotatedAxis()
markers.to.plot <- c("IL11RA")
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "group.name") + RotatedAxis()
markers.to.plot <- c("ENTPD1")
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.6, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "group.name") + RotatedAxis()

DoHeatmap(sc.kNC_kJID_hNC.integrated_2000_0.6, downsample = 100, features = features, size = 3)


#################################################


#################### add new metadata relevant to 10 clusters #############################
## label 10 clusters with datatSource, scarType, cellType, cellType_sampleName, and cellType_scarType

# retrieve the id of all cells 
cells.id <- sc.kNC_kJID_hNC.integrated_2000_0.1@assays$RNA@data@Dimnames[2][[1]]
cells.id.head <- substr(cells.id, 1, 7)
cells.id.head <- gsub("NC.", "NC", cells.id.head)
# add a new column/tag(source of data: one of the three datasets) to the metadata
sc.kNC_kJID_hNC.integrated_2000_0.1$zt.dataSource <- cells.id.head
# check batch effects of 3 datasets: review integrated scRNA-seq cells in UMAP
DimPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, reduction = "umap", shuffle = T, group.by = "zt.dataSource")
DimPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, reduction = "umap", shuffle = T, split.by = "zt.dataSource", label = T)

## Define group names according to scarType (HTS, Keloid, or NormalSkin)
# retrieve the sampleName of all cells 
cells.sampleName <- sc.kNC_kJID_hNC.integrated_2000_0.1$orig.ident
cells.scarType <- gsub("HTS.", "HTS", cells.sampleName)
cells.scarType <- gsub("Normal.", "Normal", cells.scarType)
cells.scarType <- gsub(".*CASE", "Keloid", cells.scarType)
cells.scarType <- gsub(".*CTRL", "Normal", cells.scarType)
cells.scarType <- gsub("KF.", "Keloid", cells.scarType)
cells.scarType <- gsub("NF.", "Normal", cells.scarType)
# check the cell numbers in HTS, Keloid, or NormalSkin
table(cells.scarType)
# add a new column/tag(type of scar: HTS, Keloid, or Normal) to the metadata
sc.kNC_kJID_hNC.integrated_2000_0.1$zt.scarType <- cells.scarType
# review integrated scRNA-seq cells in UMAP: split by scarType
DimPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, reduction = "umap", shuffle = T, split.by = "zt.scarType", label = T)

# Name clusters
sc.kNC_kJID_hNC.integrated_2000_0.1 <- RenameIdents(sc.kNC_kJID_hNC.integrated_2000_0.1, "0" = "FB", "1" = "EC", "2" = "spKC", "3" = "SMC", "4" = "baKC", "5" = "LL", "6" = "ML", "7" = "LEC", "8" = "MELA/NEU", "9" = "SGC")
sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType <- sc.kNC_kJID_hNC.integrated_2000_0.1@active.ident

# add a new column/tag(type of scar: HTS, Keloid, or Normal) to the metadata
sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType_scarType <- paste(Idents(sc.kNC_kJID_hNC.integrated_2000_0.1), sc.kNC_kJID_hNC.integrated_2000_0.1$zt.scarType, sep = "_")
# calculate the numbers/percentage of different cellTypes in different scarTypes (HTS, Keloid, or Normal)
cellnumbers.cellType_scarType <- table(sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType_scarType)
write.table(cellnumbers.cellType_scarType, file = "R/sc.KNC_kJID_hNC/cellnumbers.cellType_scarType.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

# add a new column/tag(type of scar: HTS, Keloid, or Normal) to the metadata
sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType_sampleName <- paste(Idents(sc.kNC_kJID_hNC.integrated_2000_0.1), sc.kNC_kJID_hNC.integrated_2000_0.1$orig.ident, sep = "_")
# calculate the numbers/percentage of different cellTypes in different samples
cellnumbers.cellType_sampleName <- table(sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType_sampleName)
write.table(cellnumbers.cellType_sampleName, file = "R/sc.KNC_kJID_hNC/cellnumbers.cellType_sampleName.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

# add a new column to the metadata, in order to calculate cellNumbers in each cellType_scarType_dataSource
load("R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")
sc.kNC_kJID_hNC.integrated_2000_0.1$zt.scarType_dataSource <- paste(sc.kNC_kJID_hNC.integrated_2000_0.1$zt.scarType, sc.kNC_kJID_hNC.integrated_2000_0.1$zt.dataSource, sep = "_")
sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType_scarType_dataSource <- paste(sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType, sc.kNC_kJID_hNC.integrated_2000_0.1$zt.scarType_dataSource, sep = "_")
cellnumbers.cellType_scarType_dataSource <- data.frame(table(sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType_scarType_dataSource))
write.table(cellnumbers.cellType_scarType_dataSource, file = "R/sc.KNC_kJID_hNC/cellnumbers.cellType_scarType_dataSource.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

save(sc.kNC_kJID_hNC.integrated_2000_0.1, file = "R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")


#################################################


#################### add new metadata relevant to 23 clusters #############################
## label 23 clusters with datatSource, scarType, 23cluster_sampleName, and 23cluster_scarType
load("R/sc.kNC_kJID_hNC.integrated2000_0.6resolution_23clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000_0.6) <- "RNA"
# add a new column/tag(source of data: one of the three datasets) to the metadata
sc.kNC_kJID_hNC.integrated_2000_0.6$zt.dataSource <- sc.kNC_kJID_hNC.integrated_2000_0.1$zt.dataSource
# add a new column/tag(type of scar: HTS, Keloid, or Normal) to the metadata
sc.kNC_kJID_hNC.integrated_2000_0.6$zt.scarType <- sc.kNC_kJID_hNC.integrated_2000_0.1$zt.scarType
# add a new column/tag(23cluster_sampleName) to the metadata
sc.kNC_kJID_hNC.integrated_2000_0.6$zt.23cluster_sampleName <- paste(Idents(sc.kNC_kJID_hNC.integrated_2000_0.6), sc.kNC_kJID_hNC.integrated_2000_0.6$orig.ident, sep = "_")

# calculate the numbers/percentage of different clusters in different scarTypes (HTS, Keloid, or Normal)
cellnumbers.23cluster_scarType <- table(sc.kNC_kJID_hNC.integrated_2000_0.6$zt.23cluster_scarType)
write.table(cellnumbers.23cluster_scarType, file = "R/sc.KNC_kJID_hNC/cellnumbers.23cluster_scarType.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

# calculate the numbers/percentage of different clusters in different samples
cellnumbers.23cluster_sampleName <- table(sc.kNC_kJID_hNC.integrated_2000_0.6$zt.23cluster_sampleName)
write.table(cellnumbers.23cluster_sampleName, file = "R/sc.KNC_kJID_hNC/cellnumbers.23cluster_sampleName.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# add a new column/tag(23cluster_scarType) to the metadata
sc.kNC_kJID_hNC.integrated_2000_0.6$zt.23cluster_scarType <- paste(Idents(sc.kNC_kJID_hNC.integrated_2000_0.6), sc.kNC_kJID_hNC.integrated_2000_0.6$zt.scarType, sep = "_")

save(sc.kNC_kJID_hNC.integrated_2000_0.6, file = "R/sc.kNC_kJID_hNC.integrated2000_0.6resolution_23clusters.Rdata")


#################################################


#################### explore a specific cluster #############################
## Dig deep into a certain cluster (e.g. FB/fibroblasts)
load("R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")
sc.kNC_kJID_hNC.FB <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "FB")
Idents(sc.kNC_kJID_hNC.FB) <- "zt.scarType"
# create a new column in the meta.data dataframe
sc.kNC_kJID_hNC.integrated_2000_0.6$celltype.groupName <- paste(Idents(sc.kNC_kJID_hNC.integrated_2000_0.6), sc.kNC_kJID_hNC.integrated_2000_0.6$group.name, sep = "_")
sc.kNC_kJID_hNC.integrated_2000_0.6$celltype <- Idents(sc.kNC_kJID_hNC.integrated_2000_0.6)

# FB_HTSvsNormal: DEGs
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType_scarType"
DEG.FB_HTSvsNormal <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = "FB_HTS", ident.2 = "FB_Normal", verbose = FALSE, logfc.threshold = 0.1)
write.table(DEG.FB_HTSvsNormal, file = "R/sc.KNC_kJID_hNC/FB/DEG.FB_HTSvsNormal.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# FB_KeloidvsNormal: DEGs
DEG.FB_KeloidvsNormal <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = "FB_Keloid", ident.2 = "FB_Normal", verbose = FALSE, logfc.threshold = 0.1)
write.table(DEG.FB_KeloidvsNormal, file = "R/sc.KNC_kJID_hNC/FB/DEG.FB_KeloidvsNormal.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# FB_HTSvsKeloid: DEGs
DEG.FB_HTSvsKeloid <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1, ident.1 = "FB_HTS", ident.2 = "FB_Keloid", verbose = FALSE, logfc.threshold = 0.1)
write.table(DEG.FB_HTSvsKeloid, file = "R/sc.KNC_kJID_hNC/FB/DEG.FB_HTSvsKeloid.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

## recluster the particular cluster/cellType (e.g. FB/fibroblasts)
# specify that we will perform downstream analysis on the corrected data note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(sc.kNC_kJID_hNC.FB) <- "integrated"
# Run the standard workflow for visualization and clustering
sc.kNC_kJID_hNC.FB <- ScaleData(sc.kNC_kJID_hNC.FB, verbose = FALSE)
sc.kNC_kJID_hNC.FB <- RunPCA(sc.kNC_kJID_hNC.FB, npcs = 30, verbose = FALSE)
sc.kNC_kJID_hNC.FB <- RunUMAP(sc.kNC_kJID_hNC.FB, reduction = "pca", dims = 1:30)
# review cells in UMAP
DimPlot(sc.kNC_kJID_hNC.FB, reduction = "umap", shuffle = T)
DimPlot(sc.kNC_kJID_hNC.FB, reduction = "umap", shuffle = T, split.by = "zt.scarType")
# find neighbors for clustering
sc.kNC_kJID_hNC.FB <- FindNeighbors(sc.kNC_kJID_hNC.FB, reduction = "pca", dims = 1:30)

# review cells in t-SNE
sc.kNC_kJID_hNC.FB.tsne <- RunTSNE(sc.kNC_kJID_hNC.FB, reduction = "pca", dims = 1:30)
DimPlot(sc.kNC_kJID_hNC.FB.tsne, reduction = "tsne", shuffle = T)
DimPlot(sc.kNC_kJID_hNC.FB.tsne, reduction = "tsne", shuffle = T, group.by = "zt.scarType")
DimPlot(sc.kNC_kJID_hNC.FB.tsne, reduction = "tsne", shuffle = T, split.by = "zt.scarType")
# find neighbors for clustering
sc.kNC_kJID_hNC.FB.tsne <- FindNeighbors(sc.kNC_kJID_hNC.FB.tsne, reduction = "pca", dims = 1:30)
sc.kNC_kJID_hNC.FB.tsne_0.168 <- FindClusters(sc.kNC_kJID_hNC.FB.tsne, resolution = 0.168)
DimPlot(sc.kNC_kJID_hNC.FB.tsne_0.168, reduction = "tsne", shuffle = T)
DimPlot(sc.kNC_kJID_hNC.FB.tsne_0.168, reduction = "tsne", shuffle = T, split.by = "zt.scarType", label = T)

#################################################


#################### identify FB_6clusters #############################
# ______Explore 6 clusters using 0.16 resolution
DefaultAssay(sc.kNC_kJID_hNC.FB) <- "integrated"
sc.kNC_kJID_hNC.FB_0.16 <- FindClusters(sc.kNC_kJID_hNC.FB, resolution = 0.16)
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, repel = T)
DimPlot(sc.kNC_kJID_hNC.FB_0.16, reduction = "umap", label = T, split.by = "zt.scarType", shuffle = T)
save(sc.kNC_kJID_hNC.FB_0.16, file = "R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")

## identify the DEGs of each cluster in FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
options(future.globals.maxSize= 3210612736)
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# cluster 0 in 6 clusters
c0.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16, ident.1 = 0, verbose = FALSE, logfc.threshold = 0.1)
write.table(c0.markers, file = "R/sc.KNC_kJID_hNC/FB/6clusters/c0.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 1 in 6 clusters
c1.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16, ident.1 = 1, verbose = FALSE, logfc.threshold = 0.1)
write.table(c1.markers, file = "R/sc.KNC_kJID_hNC/FB/6clusters/c1.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 2 in 6 clusters
c2.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16, ident.1 = 2, verbose = FALSE, logfc.threshold = 0.1)
write.table(c2.markers, file = "R/sc.KNC_kJID_hNC/FB/6clusters/c2.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 3 in 6 clusters
c3.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16, ident.1 = 3, verbose = FALSE, logfc.threshold = 0.1)
write.table(c3.markers, file = "R/sc.KNC_kJID_hNC/FB/6clusters/c3.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 4 in 6 clusters
c4.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16, ident.1 = 4, verbose = FALSE, logfc.threshold = 0.1)
write.table(c4.markers, file = "R/sc.KNC_kJID_hNC/FB/6clusters/c4.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 5 in 6 clusters
c5.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16, ident.1 = 5, verbose = FALSE, logfc.threshold = 0.1)
write.table(c5.markers, file = "R/sc.KNC_kJID_hNC/FB/6clusters/c5.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)


#################################################


#################### add new metadata to FB_6clusters #############################
## label 6 FB clusters with zt.FB6cluster_sampleName and 23cluster_sampleName
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# add a new column/tag(zt.FB6cluster_sampleName) to the metadata
sc.kNC_kJID_hNC.FB_0.16$zt.FB6cluster_sampleName <- paste(Idents(sc.kNC_kJID_hNC.FB_0.16), sc.kNC_kJID_hNC.FB_0.16$orig.ident, sep = "_")
# add a new column/tag(zt.FB6cluster_scarType) to the metadata
sc.kNC_kJID_hNC.FB_0.16$zt.FB6cluster_scarType <- paste(Idents(sc.kNC_kJID_hNC.FB_0.16), sc.kNC_kJID_hNC.FB_0.16$zt.scarType, sep = "_")

# calculate the numbers/percentage of different clusters in different samples
cellnumbers.FB6cluster_sampleName <- table(sc.kNC_kJID_hNC.FB_0.16$zt.FB6cluster_sampleName)
write.table(cellnumbers.FB6cluster_sampleName, file = "R/sc.KNC_kJID_hNC/FB/6clusters/cellnumbers.FB6cluster_sampleName.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

# calculate the numbers/percentage of different clusters in different scarTypes (HTS, Keloid, or Normal)
cellnumbers.FB6cluster_scarType <- table(sc.kNC_kJID_hNC.FB_0.16$zt.FB6cluster_scarType)
write.table(cellnumbers.FB6cluster_scarType, file = "R/sc.KNC_kJID_hNC/FB/6clusters/cellnumbers.FB6cluster_scarType.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

# add a new column to the metadata, in order to calculate cellNumbers in each FB6cluster_scarType_dataSource
sc.kNC_kJID_hNC.FB_0.16$zt.FB6cluster_scarType_dataSource <- paste(sc.kNC_kJID_hNC.FB_0.16$integrated_snn_res.0.16, sc.kNC_kJID_hNC.FB_0.16$zt.scarType_dataSource, sep = "_")
cellnumbers.FB6cluster_scarType_dataSource <- data.frame(table(sc.kNC_kJID_hNC.FB_0.16$zt.FB6cluster_scarType_dataSource))
write.table(cellnumbers.FB6cluster_scarType_dataSource, file = "R/sc.KNC_kJID_hNC/FB/6clusters/cellnumbers.FB6cluster_scarType_dataSource.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

save(sc.kNC_kJID_hNC.FB_0.16, file = "R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")


#################################################


################### expression of papillary signature ##############################
## papillary signature in PMID: 32327715
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
zt.combine_data.1 <- sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["DCN",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["CSPG4",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COL18A1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["TGFB2",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["APCDD1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["PLXNC1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["HSPB3",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COLEC12",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["NPTX2",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["LOXL3",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["ROBO2",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COL10A1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COL23A1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["S100A8",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["AXIN2",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["IL15",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["PTGDS",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["INHBB",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["PDPN",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["CTSC",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["NTF3",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COL7A1",]
sc.kNC_kJID_hNC.FB_0.16_test1.subset <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("SDC1"))
sc.kNC_kJID_hNC.FB_0.16_test1.subset@assays$RNA@data["SDC1",] <- zt.combine_data.1
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC1"), cols = c("lightblue", "red"), min.cutoff = 6, max.cutoff = 20)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC1"), cols = c("lightblue", "red"), min.cutoff = 6, max.cutoff = 20, split.by = "zt.scarType")

FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("DCN", "CSPG4", "COL18A1", "TGFB2", "APCDD1", "PLXNC1", "HSPB3", "COLEC12", "NPTX2", "LOXL3", "ROBO2", "COL10A1", "COL23A1", "S100A8", "AXIN2", "IL15", "PTGDS", "INHBB", "PDPN", "CTSC", "NTF3", "COL7A1"), 
            cols = c("lightblue", "red"), ncol = 8)
# test delete DCN
zt.combine_data.1 <- sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["CSPG4",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COL18A1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["TGFB2",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["APCDD1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["PLXNC1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["HSPB3",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COLEC12",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["NPTX2",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["LOXL3",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["ROBO2",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COL10A1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COL23A1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["S100A8",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["AXIN2",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["IL15",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["PTGDS",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["INHBB",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["PDPN",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["CTSC",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["NTF3",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COL7A1",]
sc.kNC_kJID_hNC.FB_0.16_test1.subset <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("SDC1"))
sc.kNC_kJID_hNC.FB_0.16_test1.subset@assays$RNA@data["SDC1",] <- zt.combine_data.1
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC1"), cols = c("lightblue", "red"), min.cutoff = 6, max.cutoff = 20)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC1"), cols = c("lightblue", "red"), min.cutoff = 6, max.cutoff = 20, split.by = "zt.scarType")
# 

VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC1"), pt.size = 0, combine = FALSE)
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC1"), pt.size = 0.1, combine = FALSE)
# DEGS o c1 vs c0 c3 c4
c1vs0.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 1, ident.2 = 0, verbose = FALSE, logfc.threshold = 0.1)
c1vs3.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 1, ident.2 = 3, verbose = FALSE, logfc.threshold = 0.1)
c1vs4.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 1, ident.2 = 4, verbose = FALSE, logfc.threshold = 0.1)
c1vs5.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 1, ident.2 = 5, verbose = FALSE, logfc.threshold = 0.1)
# 
c2vs0.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 2, ident.2 = 0, verbose = FALSE, logfc.threshold = 0.1)
c2vs3.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 2, ident.2 = 3, verbose = FALSE, logfc.threshold = 0.1)
c2vs4.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 2, ident.2 = 4, verbose = FALSE, logfc.threshold = 0.1)
c2vs5.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 2, ident.2 = 5, verbose = FALSE, logfc.threshold = 0.1)
# 
markers.bind <- rbind("c1vs0.markers" = c1vs0.markers, "c1vs3.markers" = c1vs3.markers, "c1vs4.markers" = c1vs4.markers, "c1vs5.markers" = c1vs5.markers, "c2vs0.markers" = c2vs0.markers, "c2vs3.markers" = c2vs3.markers, "c2vs4.markers" = c2vs4.markers, "c2vs5.markers" = c2vs5.markers)
write.table(markers.bind, file = "R/sc.KNC_kJID_hNC/FB/6clusters/papillarySignature_statistics.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)


#################################################


################### expression of reticular signature ##############################
## reticular signature in PMID: 32327715
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
zt.combine_data.1 <- sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["MFAP5",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COL11A1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["TGM2",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["IGF1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["MGP",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["PCOLCE2",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["SFRP4",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["EFEMP1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["FGF7",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["GPC4",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["ADAMTSL1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["CRLF1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["ACAN",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["PCSK5",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["BMP6",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["TAGLN",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COL14A1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["MAP1B",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["A2M",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["PPP1R14A",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["ANGPTL1",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["CDH2",] + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data["COMP",]
sc.kNC_kJID_hNC.FB_0.16_test1.subset <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("SDC1"))
sc.kNC_kJID_hNC.FB_0.16_test1.subset@assays$RNA@data["SDC1",] <- zt.combine_data.1
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC1"), cols = c("lightblue", "red"), min.cutoff = 10, max.cutoff = 25)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC1"), cols = c("lightblue", "red"), min.cutoff = 6, max.cutoff = 20, split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("MFAP5", "COL11A1", "TGM2", "IGF1", "MGP", "PCOLCE2", "SFRP4", "EFEMP1", "FGF7", "GPC4", "ADAMTSL1", "CRLF1", "ACAN", "PCSK5", "BMP6", "TAGLN", "COL14A1", "MAP1B", "A2M", "PPP1R14A", "ANGPTL1", "CDH2", "COMP"), 
            cols = c("lightblue", "red"), ncol = 8)
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC1"), pt.size = 0, combine = FALSE)
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC1"), pt.size = 0.1, combine = FALSE)
# c0 vs c1 c2 c3 c5
c0vs1.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 0, ident.2 = 1, verbose = FALSE, logfc.threshold = 0.1)
c0vs2.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 0, ident.2 = 2, verbose = FALSE, logfc.threshold = 0.1)
c0vs3.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 0, ident.2 = 3, verbose = FALSE, logfc.threshold = 0.1)
c0vs5.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 0, ident.2 = 5, verbose = FALSE, logfc.threshold = 0.1)
# c1 vs c2 c3 c5
c1vs2.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 1, ident.2 = 2, verbose = FALSE, logfc.threshold = 0.1)
c1vs3.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 1, ident.2 = 3, verbose = FALSE, logfc.threshold = 0.1)
c1vs5.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 1, ident.2 = 5, verbose = FALSE, logfc.threshold = 0.1)
# c4 vs c1 c2 c3 c5
c4vs1.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 4, ident.2 = 1, verbose = FALSE, logfc.threshold = 0.1)
c4vs2.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 4, ident.2 = 2, verbose = FALSE, logfc.threshold = 0.1)
c4vs3.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 4, ident.2 = 3, verbose = FALSE, logfc.threshold = 0.1)
c4vs5.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = 4, ident.2 = 5, verbose = FALSE, logfc.threshold = 0.1)
markers.bind <- rbind("c0vs1.markers" = c0vs1.markers, "c0vs2.markers" = c0vs2.markers, "c0vs3.markers" = c0vs3.markers, "c0vs5.markers" = c0vs5.markers, "c1vs2.markers" = c1vs2.markers, "c1vs3.markers" = c1vs3.markers, "c1vs5.markers" = c1vs5.markers, "c4vs1.markers" = c4vs1.markers, "c4vs2.markers" = c4vs2.markers, "c4vs3.markers" = c4vs3.markers, "c4vs5.markers" = c4vs5.markers)
write.table(markers.bind, file = "R/sc.KNC_kJID_hNC/FB/6clusters/reticularSignature_statistics.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

#################################################


##################### heatmap: FB_HTSorKeloidorNormal ############################
## heatmap using Seurat
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
# backup
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_forHeatmap))

# FB_KeloidvsNormal: DEGs
DEG.FB_HTSvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_HTSvsNormal.xlsx")
# pick out DEGs with p<0.05&FC>0.1
DEG.FB_HTSvsNormal_p0.05FC0.1 <- DEG.FB_HTSvsNormal[which(DEG.FB_HTSvsNormal["p_val_adj"]<0.05),]
features = top_n(DEG.FB_HTSvsNormal_p0.05FC0.1, n = 100, wt = avg_log2FC)
DoHeatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = rownames(features), slot = "scale.data", size = 3)

# FB_KeloidvsNormal: DEGs
DEG.FB_KeloidvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_KeloidvsNormal.xlsx")
# pick out DEGs with p<0.05&FC>0.1
DEG.FB_KeloidvsNormal_p0.05FC0.1 <- DEG.FB_KeloidvsNormal[which(DEG.FB_KeloidvsNormal["p_val_adj"]<0.05),]
features = top_n(DEG.FB_KeloidvsNormal_p0.05FC0.1, n = 100, wt = avg_log2FC)
DoHeatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = rownames(features), slot = "scale.data", size = 3)


## heatmap using ComplexHeatmap (Reference: COMPLEXHEATMAP展示单细胞聚类: https://www.freesion.com/article/3766164158/)
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# backup
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_forHeatmap$integrated_snn_res.0.16)

# FB_HTSvsNormal: DEGs
DEG.FB_HTSvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_HTSvsNormal.xlsx")
# pick out DEGs with p<0.05&FC>0.1
DEG.FB_HTSvsNormal_p0.05FC0.1 <- DEG.FB_HTSvsNormal[which(DEG.FB_HTSvsNormal["p_val_adj"]<0.05),]
# ___________define features, which will be subject to heatmap
features <- DEG.FB_HTSvsNormal_p0.05FC0.1 #all DEGs with p<0.05&FC>0.1
features = top_n(DEG.FB_HTSvsNormal_p0.05FC0.1, n = 100, wt = avg_log2FC) # top n DEGs with p<0.05&FC>0.1
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[rownames(features), names(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c", "#64bdc2", "#6c9cf8", "#e470dd"))
rownames(col) <- levels(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = F, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation)


## union of top50 DEGs: heatmap using ComplexHeatmap in HTS FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_HTS <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = "HTS")
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "integrated_snn_res.0.168"
Idents(sc.kNC_kJID_hNC.FB_0.16_HTS) <- "integrated_snn_res.0.168"
# backup
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16_HTS
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap$integrated_snn_res.0.16)
# define features as the uninon of top50 DEGs of HTS.FBvs.Normal.FB and Keloid.FBvsNormal.FB
union.top50 <- read.csv("R/sc.KNC_kJID_hNC/FB/union of FBHTSvsFBNormal_top50 and FBKeloidvsFBNormal_top50.txt", header = F)
features <- as.vector(union.top50[,1])
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c", "#64bdc2", "#6c9cf8", "#e470dd"))
rownames(col) <- levels(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info)
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

## union of down50 DEGs: heatmap using ComplexHeatmap in HTS FB
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap$integrated_snn_res.0.16)
# define features as the uninon of top50 DEGs of HTS.FBvs.Normal.FB and Keloid.FBvsNormal.FB
union.down50 <- read.csv("R/sc.KNC_kJID_hNC/FB/union of FBHTSvsFBNormal_down50 and FBKeloidvsFBNormal_down50.txt", header = F)
features <- as.vector(union.down50[,1])
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c", "#64bdc2", "#6c9cf8", "#e470dd"))
rownames(col) <- levels(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info)
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


## union of top50 DEGs: heatmap using ComplexHeatmap in Keloid FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = "Keloid")
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "integrated_snn_res.0.168"
Idents(sc.kNC_kJID_hNC.FB_0.16_Keloid) <- "integrated_snn_res.0.168"
# backup
sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16_Keloid
sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap$integrated_snn_res.0.16)
# define features as the uninon of top50 DEGs of HTS.FBvs.Normal.FB and Keloid.FBvsNormal.FB
union.top50 <- read.csv("R/sc.KNC_kJID_hNC/FB/union of FBHTSvsFBNormal_top50 and FBKeloidvsFBNormal_top50.txt", header = F)
features <- as.vector(union.top50[,1])
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.cluster_info))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c", "#64bdc2", "#6c9cf8", "#e470dd"))
rownames(col) <- levels(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.cluster_info)
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

## union of down50 DEGs: heatmap using ComplexHeatmap in Keloid FB
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap$integrated_snn_res.0.16)
# define features as the uninon of top50 DEGs of HTS.FBvs.Normal.FB and Keloid.FBvsNormal.FB
union.down50 <- read.csv("R/sc.KNC_kJID_hNC/FB/union of FBHTSvsFBNormal_down50 and FBKeloidvsFBNormal_down50.txt", header = F)
features <- as.vector(union.down50[,1])
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.cluster_info))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c", "#64bdc2", "#6c9cf8", "#e470dd"))
rownames(col) <- levels(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.cluster_info)
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_Keloid_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


## union of top50 DEGs: heatmap using ComplexHeatmap in Normal FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_Normal <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = "Normal")
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "integrated_snn_res.0.168"
Idents(sc.kNC_kJID_hNC.FB_0.16_Normal) <- "integrated_snn_res.0.168"
# backup
sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16_Normal
sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap$integrated_snn_res.0.16)
# define features as the uninon of top50 DEGs of HTS.FBvs.Normal.FB and Keloid.FBvsNormal.FB
union.top50 <- read.csv("R/sc.KNC_kJID_hNC/FB/union of FBHTSvsFBNormal_top50 and FBKeloidvsFBNormal_top50.txt", header = F)
features <- as.vector(union.top50[,1])
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.cluster_info))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c", "#64bdc2", "#6c9cf8", "#e470dd"))
rownames(col) <- levels(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.cluster_info)
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

## union of down50 DEGs: heatmap using ComplexHeatmap in Normal FB
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap$integrated_snn_res.0.16)
# define features as the uninon of top50 DEGs of HTS.FBvs.Normal.FB and Keloid.FBvsNormal.FB
union.down50 <- read.csv("R/sc.KNC_kJID_hNC/FB/union of FBHTSvsFBNormal_down50 and FBKeloidvsFBNormal_down50.txt", header = F)
features <- as.vector(union.down50[,1])
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.cluster_info))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c", "#64bdc2", "#6c9cf8", "#e470dd"))
rownames(col) <- levels(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.cluster_info)
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_Normal_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


#################################################


#################### Heatmap of FB_HTSorKeloidorNormal: selected DEGs in HTS and/or Keloid #############################
## union.top50
# union.top50_HTSunique
union.top50_HTSunique <- c("HBB", "HBA2", "C1QTNF3", "IFI27", "RPL13A", "RPS20", "TM4SF1", "TMSB4X", "HLA-DRA", "C12orf75", "IFI6", "RPLP2", "CRABP2", "MMP2", "CTGF", "PPIC", "TPM1", "KDELR3", "LGALS1", "ATP5F1E", "CD24", "KDELR2", "SERPINH1", "CAVIN3")
# union.top50_common
union.top50_common <- c("POSTN", "ASPN", "COMP", "SFRP4", "COL3A1", "BGN", "SPARC", "MDK", "COL5A2", "LUM", "COL1A1", "CTHRC1", "FN1", "COL1A2", "ADAM12", "THBS4", "SFRP2", "ELN", "TAGLN", "COL14A1", "SPON2", "MMP23B", "CERCAM")
# union.top50_Keloidunique
union.top50_Keloidunique <- c("COL5A1", "PRSS23", "HTRA1", "NREP", "COL6A3", "OGN", "CPXM1", "ARL4C", "COPZ2", "FAP", "COL6A1", "TNC", "VCAN", "C1QTNF6", "OLFML2B", "COL6A2", "COL11A1", "BASP1", "SULF2", "MARCKS", "THY1", "GOLM1", "EDIL3", "LRRC15", "TPM2", "CHPF", "DCD")

## union.down50
# union.down50_HTSunique
union.down50_HTSunique <- c("GEM", "CXCL3", "CCNL1", "SLC38A2", "IER5", "TXNIP", "IGFBP7", "IER3", "CXCL2", "TUBB4B", "WTAP", "GADD45A", "RASD1", "RRAD", "DNAJA1", "HSP90AB1", "ABCA8", "NFIA", "TNFAIP6", "UBE2S", "NFKBIA", "APOE", "GADD45B", "SOD2", "JUND", "SFRP1", "DNAJB1", "EIF4A3", "HSPA6", "CCL19", "IER5L", "TWIST1", "HSP90AA1", "CD81", "HSPA1B", "ID1", "PTGDS", "C3", "HSPA1A", "PLCG2", "APOD", "C11orf96")
# union.down50_common
union.down50_common <- c("JUNB", "HMOX1", "SELM", "CXCL1", "GPX3", "PLA2G2A", "SEPP1", "GNB2L1")
# union.down50_Keloidunique
union.down50_Keloidunique <- c("CHP1", "GPNMB", "HLA-B", "FCGRT", "HLA-C", "MFGE8", "ADIRF", "ANXA2", "MYC", "TPPP3", "ZFP36L2", "CLDN11", "ZFP36", "LEPR", "FBLN1", "CES1", "RPS17", "CLU", "PCSK1N", "ALDOA", "MT1M", "TSPAN8", "PLPP3", "DCN", "F10", "KRT14", "CEBPD", "TSC22D3", "TNXB", "SOD3", "RPS26", "SLPI", "CD9", "CCL2", "PI16", "CFD", "ADH1B", "AXL", "GSN")
# union_HTSupKeloiddown
union_HTSupKeloiddown <- c("RPL7", "PCOLCE2", "RPL31")

## union.top50_down50
union.top50_down50 <- c(union.top50_HTSunique, union.top50_common, union.top50_Keloidunique, union.down50_HTSunique, union.down50_common, union.down50_Keloidunique, union_HTSupKeloiddown)


## union of top50 DEGs: heatmap using ComplexHeatmap in FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
# backup
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_forHeatmap$zt.scarType)
# define features
features <- union.top50_down50
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


# FB_HTSvsNormal: DEGs_p<0.05
DEG.FB_HTSvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_HTSvsNormal.xlsx")
DEG.FB_HTSvsNormal_p0.05 <- subset(DEG.FB_HTSvsNormal, p_val_adj < 0.05)
DEG.FB_HTSvsNormal_p0.05 <- DEG.FB_HTSvsNormal_p0.05[order(DEG.FB_HTSvsNormal_p0.05$avg_log2FC, decreasing = T), ]
DEG.FB_HTSvsNormal_p0.05 <- rbind(head(DEG.FB_HTSvsNormal_p0.05, 50), tail(DEG.FB_HTSvsNormal_p0.05, 50))
# FB_KeloidvsNormal: DEGs_p<0.05
DEG.FB_KeloidvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_KeloidvsNormal.xlsx")
DEG.FB_KeloidvsNormal_p0.05 <- subset(DEG.FB_KeloidvsNormal, p_val_adj < 0.05)
DEG.FB_KeloidvsNormal_p0.05 <- DEG.FB_KeloidvsNormal_p0.05[order(DEG.FB_KeloidvsNormal_p0.05$avg_log2FC, decreasing = T), ]
DEG.FB_KeloidvsNormal_p0.05 <- rbind(head(DEG.FB_KeloidvsNormal_p0.05, 50), tail(DEG.FB_KeloidvsNormal_p0.05, 50))
# all DEGs: the union of DEG.FB_HTSvsNormal_p0.05 and DEG.FB_KeloidvsNormal_p0.05
DEG.FB_HTSandKeloidvsNormal_p0.05 <- union(rownames(DEG.FB_HTSvsNormal_p0.05), rownames(DEG.FB_KeloidvsNormal_p0.05))
## heatmap of all DEGs: heatmap using ComplexHeatmap in FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
# backup
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_forHeatmap$zt.scarType)
# define features
features <- DEG.FB_HTSandKeloidvsNormal_p0.05
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = F, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))



# FB_HTSvsNormal: DEGs_p<0.05
DEG.FB_HTSvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_HTSvsNormal.xlsx")
DEG.FB_HTSvsNormal_p0.05 <- subset(DEG.FB_HTSvsNormal, p_val_adj < 0.05)
DEG.FB_HTSvsNormal_p0.05 <- DEG.FB_HTSvsNormal_p0.05[order(DEG.FB_HTSvsNormal_p0.05$avg_log2FC, decreasing = T), ]
DEG.FB_HTSvsNormal_p0.05 <- rbind(head(DEG.FB_HTSvsNormal_p0.05, 1000), tail(DEG.FB_HTSvsNormal_p0.05, 1000))
# all DEGs: the union of DEG.FB_HTSvsNormal_p0.05 and DEG.FB_KeloidvsNormal_p0.05
DEG.FB_HTSandKeloidvsNormal_p0.05 <- rownames(head(DEG.FB_HTSvsNormal_p0.05, 100))
DEG.FB_HTSandKeloidvsNormal_p0.05 <- rownames(DEG.FB_HTSvsNormal_p0.05)
## heatmap of all DEGs: heatmap using ComplexHeatmap in FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
# backup
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = c("HTS", "Normal"))
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_forHeatmap$zt.scarType)
# define features
features <- DEG.FB_HTSandKeloidvsNormal_p0.05
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix, 
        cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = F, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))





# FB_HTSvsNormal: DEGs_p<0.05
DEG.FB_HTSvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_HTSvsNormal.xlsx")
DEG.FB_HTSvsNormal_p0.05 <- subset(DEG.FB_HTSvsNormal, p_val_adj < 0.05)
DEG.FB_HTSvsNormal_p0.05 <- DEG.FB_HTSvsNormal_p0.05[order(DEG.FB_HTSvsNormal_p0.05$avg_log2FC, decreasing = T), ]
DEG.FB_HTSvsNormal_p0.05 <- head(DEG.FB_HTSvsNormal_p0.05, 100)
# FB_KeloidvsNormal: DEGs_p<0.05
DEG.FB_KeloidvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_KeloidvsNormal.xlsx")
DEG.FB_KeloidvsNormal_p0.05 <- subset(DEG.FB_KeloidvsNormal, p_val_adj < 0.05)
DEG.FB_KeloidvsNormal_p0.05 <- DEG.FB_KeloidvsNormal_p0.05[order(DEG.FB_KeloidvsNormal_p0.05$avg_log2FC, decreasing = T), ]
DEG.FB_KeloidvsNormal_p0.05 <- head(DEG.FB_KeloidvsNormal_p0.05, 100)
# all DEGs: the union of DEG.FB_HTSvsNormal_p0.05 and DEG.FB_KeloidvsNormal_p0.05
DEG.FB_HTSandKeloidvsNormal_p0.05 <- union(rownames(DEG.FB_HTSvsNormal_p0.05), rownames(DEG.FB_KeloidvsNormal_p0.05))
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
# backup
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- subset(sc.kNC_kJID_hNC.FB_0.16)
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_forHeatmap))

## DoHeatmap
# define features
features <- DEG.FB_HTSandKeloidvsNormal_p0.05
DoHeatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = features, slot = "scale.data", size = 3)

## ComplexHeatmap
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_forHeatmap$zt.scarType)
# define features
features <- DEG.FB_HTSandKeloidvsNormal_p0.05
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = F, show_column_names = F, show_row_names = F, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5), 
        col = colorRamp2(c(-2, 0, 2), c("purple","black","yellow")))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = F, show_column_names = F, show_row_names = F, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5), 
        col = colorRamp2(seq(min(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix), max(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix), length=3), c("purple","black","yellow")))


#################################################


#################### For heatmap of FB_6 clusters: selected DEGs in HTS and/or Keloid #############################
## Define features for heatmap
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# FB_HTSvsNormal: DEGs
DEG.FB_HTSvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_HTSvsNormal.xlsx")
# FB_KeloidvsNormal: DEGs
DEG.FB_KeloidvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_KeloidvsNormal.xlsx")

## union.top50
# union.top50_HTSunique
union.top50_HTSunique <- c("HBB", "HBA2", "C1QTNF3", "IFI27", "RPL13A", "RPS20", "TM4SF1", "TMSB4X", "HLA-DRA", "C12orf75", "IFI6", "RPLP2", "CRABP2", "MMP2", "CTGF", "PPIC", "TPM1", "KDELR3", "LGALS1", "ATP5F1E", "CD24", "KDELR2", "SERPINH1", "CAVIN3")
# union.top50_common
union.top50_common <- c("POSTN", "ASPN", "COMP", "SFRP4", "COL3A1", "BGN", "SPARC", "MDK", "COL5A2", "LUM", "COL1A1", "CTHRC1", "FN1", "COL1A2", "ADAM12", "THBS4", "SFRP2", "ELN", "TAGLN", "COL14A1", "SPON2", "MMP23B", "CERCAM")
# union.top50_Keloidunique
union.top50_Keloidunique <- c("COL5A1", "PRSS23", "HTRA1", "NREP", "COL6A3", "OGN", "CPXM1", "ARL4C", "COPZ2", "FAP", "COL6A1", "TNC", "VCAN", "C1QTNF6", "OLFML2B", "COL6A2", "COL11A1", "BASP1", "SULF2", "MARCKS", "THY1", "GOLM1", "EDIL3", "LRRC15", "TPM2", "CHPF", "DCD")

## union.down50
# union.down50_HTSunique
union.down50_HTSunique <- c("GEM", "CXCL3", "CCNL1", "SLC38A2", "IER5", "TXNIP", "IGFBP7", "IER3", "CXCL2", "TUBB4B", "WTAP", "GADD45A", "RASD1", "RRAD", "DNAJA1", "HSP90AB1", "ABCA8", "NFIA", "TNFAIP6", "UBE2S", "NFKBIA", "APOE", "GADD45B", "SOD2", "JUND", "SFRP1", "DNAJB1", "EIF4A3", "HSPA6", "CCL19", "IER5L", "TWIST1", "HSP90AA1", "CD81", "HSPA1B", "ID1", "PTGDS", "C3", "HSPA1A", "PLCG2", "APOD", "C11orf96")
# union.down50_common
union.down50_common <- c("JUNB", "HMOX1", "SELM", "CXCL1", "GPX3", "PLA2G2A", "SEPP1", "GNB2L1")
# union.down50_Keloidunique
union.down50_Keloidunique <- c("CHP1", "GPNMB", "HLA-B", "FCGRT", "HLA-C", "MFGE8", "ADIRF", "ANXA2", "MYC", "TPPP3", "ZFP36L2", "CLDN11", "ZFP36", "LEPR", "FBLN1", "CES1", "RPS17", "CLU", "PCSK1N", "ALDOA", "MT1M", "TSPAN8", "PLPP3", "DCN", "F10", "KRT14", "CEBPD", "TSC22D3", "TNXB", "SOD3", "RPS26", "SLPI", "CD9", "CCL2", "PI16", "CFD", "ADH1B", "AXL", "GSN")
# union_HTSupKeloiddown
union_HTSupKeloiddown <- c("RPL7", "PCOLCE2", "RPL31")

## union.top50_down50
union.top50_down50 <- c(union.top50_HTSunique, union.top50_common, union.top50_Keloidunique, union.down50_HTSunique, union.down50_common, union.down50_Keloidunique, union_HTSupKeloiddown)


#################################################


##################### ScarsignatureGenes_FB6clusters: GeneOverlap ############################
library(GeneOverlap)
data(GeneOverlap)
# Gene signature for ordering: the union of DEGs in HTSvsNormal and KeloidvsNormal
# FB_HTSvsNormal: DEGs_p<0.05
DEG.FB_HTSvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_HTSvsNormal.xlsx")
DEG.FB_HTSvsNormal_p0.05 <- subset(DEG.FB_HTSvsNormal, p_val_adj < 0.05) %>% rownames()
# FB_KeloidvsNormal: DEGs_p<0.05
DEG.FB_KeloidvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_KeloidvsNormal.xlsx")
DEG.FB_KeloidvsNormal_p0.05 <- subset(DEG.FB_KeloidvsNormal, p_val_adj < 0.05) %>% rownames() %>% as.vector()
# Define the backgroundGenome
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
backgroundGenome <- sc.kNC_kJID_hNC.FB_0.16@assays$RNA@data@Dimnames[[1]]
# test
testOverlap <- newGeneOverlap(DEG.FB_HTSvsNormal_p0.05, DEG.FB_KeloidvsNormal_p0.05, genome.size = gs.RNASeq)
testOverlap <- testGeneOverlap(testOverlap)
print(testOverlap)


#################################################


##################### ScarsignatureGenes_FB6clusters: violinPlot ############################
## scarSignature
scarSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/scarSignature.txt",
                      header = T, sep = "\t")
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
# construct an expressionMatrix for scarSignature as a single gene
zt.combine_data.1 <- sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data[scarSig[1,],]
for (i in 2:dim(scarSig)[1]) {
  zt.combine_data.1 <- zt.combine_data.1 + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data[scarSig[i,],]
}
sc.kNC_kJID_hNC.FB_0.16_test1.subset <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("SDC3"))
sc.kNC_kJID_hNC.FB_0.16_test1.subset@assays$RNA@data["SDC3",] <- zt.combine_data.1
# featurePlot of scarSignature
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("lightblue", "red"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("lightblue", "red"), split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("lightblue", "red"), split.by = "zt.scarType_dataSource")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), pt.size = 0, combine = FALSE, group.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), pt.size = 1, combine = FALSE, group.by = "zt.scarType_dataSource")
DotPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("blue", "red"), group.by = "zt.scarType") + RotatedAxis()
DotPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("blue", "red"), group.by = "zt.scarType_dataSource") + RotatedAxis()
# violin plots of scarSignature in FBc0-5, splitbyScarType
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), pt.size = 0, combine = FALSE, split.by = "zt.scarType")
# or (separated)
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "0", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "1", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "2", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "3", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "4", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "5", pt.size = 0, combine = FALSE, split.by = "zt.scarType")

# statistical analysis between scarTypes
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.scarType"
HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "HTS", ident.2 = "Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(HTSvsNormal.markers)
KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Keloid", ident.2 = "Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(KeloidvsNormal.markers)
HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "HTS", ident.2 = "Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(HTSvsKeloid.markers)
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.scarType_dataSource"
Normal_sc.hNCvsNormal_sc.kJID.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Normal_sc.hNC", ident.2 = "Normal_sc.kJID", verbose = FALSE, logfc.threshold = 0.0001)
print(Normal_sc.hNCvsNormal_sc.kJID.markers)
Normal_sc.hNCvsNormal_sc.kNC.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Normal_sc.hNC", ident.2 = "Normal_sc.kNC", verbose = FALSE, logfc.threshold = 0.0001)
print(Normal_sc.hNCvsNormal_sc.kNC.markers)
Normal_sc.kJIDvsNormal_sc.kNC.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Normal_sc.kJID", ident.2 = "Normal_sc.kNC", verbose = FALSE, logfc.threshold = 0.0001)
print(Normal_sc.kJIDvsNormal_sc.kNC.markers)

# statistical analysis between scarTypes in FB6clusters 
# FBc0
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc0.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "0_HTS", ident.2 = "0_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc0.HTSvsNormal.markers)
FBc0.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "0_Keloid", ident.2 = "0_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc0.KeloidvsNormal.markers)
FBc0.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "0_HTS", ident.2 = "0_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc0.HTSvsKeloid.markers)
# FBc1
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc1.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "1_HTS", ident.2 = "1_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc1.HTSvsNormal.markers)
FBc1.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "1_Keloid", ident.2 = "1_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc1.KeloidvsNormal.markers)
FBc1.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "1_HTS", ident.2 = "1_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc1.HTSvsKeloid.markers)
# FBc2
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc2.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "2_HTS", ident.2 = "2_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc2.HTSvsNormal.markers)
FBc2.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "2_Keloid", ident.2 = "2_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc2.KeloidvsNormal.markers)
FBc2.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "2_HTS", ident.2 = "2_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc2.HTSvsKeloid.markers)
# FBc3
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc3.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "3_HTS", ident.2 = "3_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc3.HTSvsNormal.markers)
FBc3.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "3_Keloid", ident.2 = "3_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc3.KeloidvsNormal.markers)
FBc3.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "3_HTS", ident.2 = "3_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc3.HTSvsKeloid.markers)
# FBc4
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc4.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "4_HTS", ident.2 = "4_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc4.HTSvsNormal.markers)
FBc4.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "4_Keloid", ident.2 = "4_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc4.KeloidvsNormal.markers)
FBc4.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "4_HTS", ident.2 = "4_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc4.HTSvsKeloid.markers)
# FBc5
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc5.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "5_HTS", ident.2 = "5_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc5.HTSvsNormal.markers)
FBc5.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "5_Keloid", ident.2 = "5_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc5.KeloidvsNormal.markers)
FBc5.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "5_HTS", ident.2 = "5_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc5.HTSvsKeloid.markers)



## htsSignature
htsSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/htsSignature.txt",
                     header = T, sep = "\t")
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
# construct an expressionMatrix for htsSignature as a single gene
zt.combine_data.1 <- sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data[htsSig[1,],]
for (i in 2:dim(htsSig)[1]) {
  zt.combine_data.1 <- zt.combine_data.1 + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data[htsSig[i,],]
}
sc.kNC_kJID_hNC.FB_0.16_test1.subset <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("SDC3"))
sc.kNC_kJID_hNC.FB_0.16_test1.subset@assays$RNA@data["SDC3",] <- zt.combine_data.1
# featurePlot of scarSignature
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("lightblue", "red"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("lightblue", "red"), split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("lightblue", "red"), split.by = "zt.scarType_dataSource")
DotPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("blue", "red"), group.by = "zt.scarType") + RotatedAxis()
DotPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("blue", "red"), group.by = "zt.scarType_dataSource") + RotatedAxis()
# violin plots of scarSignature in FBc0-5, splitbyScarType
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), pt.size = 0, combine = FALSE, group.by = "zt.scarType_dataSource")
# or (separated)
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "0", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "1", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "2", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "3", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "4", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "5", pt.size = 0, combine = FALSE, split.by = "zt.scarType")

# statistical analysis between scarTypes
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.scarType"
HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "HTS", ident.2 = "Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(HTSvsNormal.markers)
KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Keloid", ident.2 = "Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(KeloidvsNormal.markers)
HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "HTS", ident.2 = "Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(HTSvsKeloid.markers)
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.scarType_dataSource"
Normal_sc.hNCvsNormal_sc.kJID.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Normal_sc.hNC", ident.2 = "Normal_sc.kJID", verbose = FALSE, logfc.threshold = 0.0001)
print(Normal_sc.hNCvsNormal_sc.kJID.markers)
Normal_sc.hNCvsNormal_sc.kNC.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Normal_sc.hNC", ident.2 = "Normal_sc.kNC", verbose = FALSE, logfc.threshold = 0.0001)
print(Normal_sc.hNCvsNormal_sc.kNC.markers)
Normal_sc.kJIDvsNormal_sc.kNC.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Normal_sc.kJID", ident.2 = "Normal_sc.kNC", verbose = FALSE, logfc.threshold = 0.0001)
print(Normal_sc.kJIDvsNormal_sc.kNC.markers)

# statistical analysis between scarTypes in FB6clusters 
# FBc0
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc0.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "0_HTS", ident.2 = "0_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc0.HTSvsNormal.markers)
FBc0.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "0_Keloid", ident.2 = "0_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc0.KeloidvsNormal.markers)
FBc0.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "0_HTS", ident.2 = "0_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc0.HTSvsKeloid.markers)
# FBc1
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc1.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "1_HTS", ident.2 = "1_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc1.HTSvsNormal.markers)
FBc1.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "1_Keloid", ident.2 = "1_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc1.KeloidvsNormal.markers)
FBc1.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "1_HTS", ident.2 = "1_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc1.HTSvsKeloid.markers)
# FBc2
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc2.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "2_HTS", ident.2 = "2_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc2.HTSvsNormal.markers)
FBc2.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "2_Keloid", ident.2 = "2_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc2.KeloidvsNormal.markers)
FBc2.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "2_HTS", ident.2 = "2_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc2.HTSvsKeloid.markers)
# FBc3
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc3.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "3_HTS", ident.2 = "3_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc3.HTSvsNormal.markers)
FBc3.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "3_Keloid", ident.2 = "3_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc3.KeloidvsNormal.markers)
FBc3.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "3_HTS", ident.2 = "3_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc3.HTSvsKeloid.markers)
# FBc4
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc4.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "4_HTS", ident.2 = "4_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc4.HTSvsNormal.markers)
FBc4.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "4_Keloid", ident.2 = "4_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc4.KeloidvsNormal.markers)
FBc4.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "4_HTS", ident.2 = "4_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc4.HTSvsKeloid.markers)
# FBc5
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc5.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "5_HTS", ident.2 = "5_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc5.HTSvsNormal.markers)
FBc5.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "5_Keloid", ident.2 = "5_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc5.KeloidvsNormal.markers)
FBc5.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "5_HTS", ident.2 = "5_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc5.HTSvsKeloid.markers)



## keloidSignature
keloidSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/keloidSignature.txt",
                     header = T, sep = "\t")
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
# construct an expressionMatrix for htsSignature as a single gene
zt.combine_data.1 <- sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data[keloidSig[1,],]
for (i in 2:dim(keloidSig)[1]) {
  zt.combine_data.1 <- zt.combine_data.1 + sc.kNC_kJID_hNC.FB_0.16_test1@assays$RNA@data[keloidSig[i,],]
}
sc.kNC_kJID_hNC.FB_0.16_test1.subset <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("SDC3"))
sc.kNC_kJID_hNC.FB_0.16_test1.subset@assays$RNA@data["SDC3",] <- zt.combine_data.1
# featurePlot of scarSignature
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("lightblue", "red"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("lightblue", "red"), split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("lightblue", "red"), split.by = "zt.scarType_dataSource")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), pt.size = 0, combine = FALSE, group.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), pt.size = 0, combine = FALSE, group.by = "zt.scarType_dataSource")
DotPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("blue", "red"), group.by = "zt.scarType") + RotatedAxis()
DotPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), cols = c("blue", "red"), group.by = "zt.scarType_dataSource") + RotatedAxis()
# violin plots of scarSignature in FBc0-5, splitbyScarType
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), pt.size = 0, combine = FALSE, split.by = "zt.scarType")
# or (separated)
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "0", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "1", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "2", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "3", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "4", pt.size = 0, combine = FALSE, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), idents = "5", pt.size = 0, combine = FALSE, split.by = "zt.scarType")

# statistical analysis between scarTypes
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.scarType"
HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "HTS", ident.2 = "Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(HTSvsNormal.markers)
KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Keloid", ident.2 = "Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(KeloidvsNormal.markers)
HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "HTS", ident.2 = "Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(HTSvsKeloid.markers)
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.scarType_dataSource"
Normal_sc.hNCvsNormal_sc.kJID.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Normal_sc.hNC", ident.2 = "Normal_sc.kJID", verbose = FALSE, logfc.threshold = 0.0001)
print(Normal_sc.hNCvsNormal_sc.kJID.markers)
Normal_sc.hNCvsNormal_sc.kNC.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Normal_sc.hNC", ident.2 = "Normal_sc.kNC", verbose = FALSE, logfc.threshold = 0.0001)
print(Normal_sc.hNCvsNormal_sc.kNC.markers)
Normal_sc.kJIDvsNormal_sc.kNC.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "Normal_sc.kJID", ident.2 = "Normal_sc.kNC", verbose = FALSE, logfc.threshold = 0.0001)
print(Normal_sc.kJIDvsNormal_sc.kNC.markers)

# statistical analysis between scarTypes in FB6clusters 
# FBc0
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc0.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "0_HTS", ident.2 = "0_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc0.HTSvsNormal.markers)
FBc0.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "0_Keloid", ident.2 = "0_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc0.KeloidvsNormal.markers)
FBc0.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "0_HTS", ident.2 = "0_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc0.HTSvsKeloid.markers)
# FBc1
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc1.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "1_HTS", ident.2 = "1_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc1.HTSvsNormal.markers)
FBc1.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "1_Keloid", ident.2 = "1_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc1.KeloidvsNormal.markers)
FBc1.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "1_HTS", ident.2 = "1_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc1.HTSvsKeloid.markers)
# FBc2
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc2.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "2_HTS", ident.2 = "2_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc2.HTSvsNormal.markers)
FBc2.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "2_Keloid", ident.2 = "2_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc2.KeloidvsNormal.markers)
FBc2.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "2_HTS", ident.2 = "2_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc2.HTSvsKeloid.markers)
# FBc3
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc3.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "3_HTS", ident.2 = "3_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc3.HTSvsNormal.markers)
FBc3.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "3_Keloid", ident.2 = "3_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc3.KeloidvsNormal.markers)
FBc3.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "3_HTS", ident.2 = "3_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc3.HTSvsKeloid.markers)
# FBc4
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc4.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "4_HTS", ident.2 = "4_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc4.HTSvsNormal.markers)
FBc4.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "4_Keloid", ident.2 = "4_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc4.KeloidvsNormal.markers)
FBc4.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "4_HTS", ident.2 = "4_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc4.HTSvsKeloid.markers)
# FBc5
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.FB6cluster_scarType"
FBc5.HTSvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "5_HTS", ident.2 = "5_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc5.HTSvsNormal.markers)
FBc5.KeloidvsNormal.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "5_Keloid", ident.2 = "5_Normal", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc5.KeloidvsNormal.markers)
FBc5.HTSvsKeloid.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.16_test1.subset, ident.1 = "5_HTS", ident.2 = "5_Keloid", verbose = FALSE, logfc.threshold = 0.0001)
print(FBc5.HTSvsKeloid.markers)


#################################################


####################### ScarsignatureGenes_FB: AddModuleScore in LSD::Heatscatter ##########################
#* AddModuleScore is a function in the Seurat R package AND is used to calculate the average expression level of a pre-determined gene set (reference: PMID_35476447)
## scarSignature
scarSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/scarSignature.txt",
                      header = T, sep = "\t")
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# backup
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
# calculate ModuleScore of every cell
sc.kNC_kJID_hNC.FB_0.16_test1 <- AddModuleScore(sc.kNC_kJID_hNC.FB_0.16_test1, features = scarSig, name = "moduleScore.scarSignature")
colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data)[length(colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data))] <- 'moduleScore.scarSignature' 
# heatscatter background
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = "SDC3", cols = c("lightgrey", "lightgrey")) + xlim(-8, 8) + ylim(-8, 8)
# add a new column (in/activatied scarSignature) to the metadata
sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.scarSignature <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature)
threshold.activated <- quantile(sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature)[4]
for (i in 1:length(sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.scarSignature)) {
  if (sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.scarSignature[i] > threshold.activated) {
    sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.scarSignature[i] <-  "activated"
  } else {
    sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.scarSignature[i] <-  "inactivated"
  }
}
sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature <- sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.scarSignature
# heatscatter of activated FB
sc.kNC_kJID_hNC.FB_0.16_test1.activated <- sc.kNC_kJID_hNC.FB_0.16_test1
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "moduleScore.scarSignature"
sc.kNC_kJID_hNC.FB_0.16_test1.activated <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "activated")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "FB")
# heatscatter of activated FB.HTS
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.activated.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "HTS")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated.HTS@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated.HTS@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "HTS")
# heatscatter of activated FB.Keloid
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.activated.Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "Keloid")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated.Keloid@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated.Keloid@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Keloid")
# heatscatter of activated FB.Normal
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.activated.Normal <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "Normal")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated.Normal@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated.Normal@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Normal")


## htsSignature
htsSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/htsSignature.txt",
                        header = T, sep = "\t")
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# backup
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
# calculate ModuleScore of every cell
sc.kNC_kJID_hNC.FB_0.16_test1 <- AddModuleScore(sc.kNC_kJID_hNC.FB_0.16_test1, features = htsSig, name = "moduleScore.htsSignature")
colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data)[length(colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data))] <- 'moduleScore.htsSignature' 
# heatscatter background
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = "SDC3", cols = c("lightgrey", "lightgrey"))
# add a new column (in/activatied htsSignature) to the metadata
sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.htsSignature <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.htsSignature)
threshold.activated <- quantile(sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.htsSignature)[4]
for (i in 1:length(sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.htsSignature)) {
  if (sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.htsSignature[i] > threshold.activated) {
    sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.htsSignature[i] <-  "activated"
  } else {
    sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.htsSignature[i] <-  "inactivated"
  }
}
sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.htsSignature <- sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.htsSignature
# heatscatter of activated FB
sc.kNC_kJID_hNC.FB_0.16_test1.activated <- sc.kNC_kJID_hNC.FB_0.16_test1
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "moduleScore.htsSignature"
sc.kNC_kJID_hNC.FB_0.16_test1.activated <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "activated")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "FB")
# heatscatter of activated FB.HTS
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.activated.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "HTS")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated.HTS@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated.HTS@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "HTS")
# heatscatter of activated FB.Keloid
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.activated.Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "Keloid")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated.Keloid@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated.Keloid@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Keloid")
# heatscatter of activated FB.Normal
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.activated.Normal <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "Normal")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated.Normal@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated.Normal@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Normal")


## keloidSignature
keloidSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/keloidSignature.txt",
                      header = T, sep = "\t")
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# backup
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
# calculate ModuleScore of every cell
sc.kNC_kJID_hNC.FB_0.16_test1 <- AddModuleScore(sc.kNC_kJID_hNC.FB_0.16_test1, features = keloidSig, name = "moduleScore.keloidSignature")
colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data)[length(colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data))] <- 'moduleScore.keloidSignature' 
# heatscatter background
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = "SDC3", cols = c("lightgrey", "lightgrey"))
# add a new column (in/activatied keloidSignature) to the metadata
sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.keloidSignature <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.keloidSignature)
threshold.activated <- quantile(sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.keloidSignature)[4]
for (i in 1:length(sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.keloidSignature)) {
  if (sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.keloidSignature[i] > threshold.activated) {
    sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.keloidSignature[i] <-  "activated"
  } else {
    sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.keloidSignature[i] <-  "inactivated"
  }
}
sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.keloidSignature <- sc.kNC_kJID_hNC.FB_0.16_test1.moduleScore.keloidSignature
# heatscatter of activated FB
sc.kNC_kJID_hNC.FB_0.16_test1.activated <- sc.kNC_kJID_hNC.FB_0.16_test1
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "moduleScore.keloidSignature"
sc.kNC_kJID_hNC.FB_0.16_test1.activated <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "activated")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "FB")
# heatscatter of activated FB.HTS
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.activated.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "HTS")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated.HTS@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated.HTS@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "HTS")
# heatscatter of activated FB.Keloid
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.activated.Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "Keloid")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated.Keloid@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated.Keloid@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Keloid")
# heatscatter of activated FB.Normal
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.activated) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_test1.activated.Normal <- subset(sc.kNC_kJID_hNC.FB_0.16_test1.activated, idents = "Normal")
heatscatter(sc.kNC_kJID_hNC.FB_0.16_test1.activated.Normal@reductions$umap@cell.embeddings[, 1], sc.kNC_kJID_hNC.FB_0.16_test1.activated.Normal@reductions$umap@cell.embeddings[, 2],
            xlab = "UMAP1", ylab= "UMAP2", xlim = c(-8, 8), ylim = c(-8, 8),
            colpal= "bl2gr2rd",
            add.contour = T, main = "Normal")


#################################################


####################### ScarsignatureGenes_FB: expression distribution by kernel density estimation: Nebulosa ##########################
## scarSignature
scarSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/scarSignature.txt",
                      header = T, sep = "\t")
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# backup
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
# calculate ModuleScore of every cell
sc.kNC_kJID_hNC.FB_0.16_test1 <- AddModuleScore(sc.kNC_kJID_hNC.FB_0.16_test1, features = scarSig, name = "moduleScore.scarSignature")
colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data)[length(colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data))] <- 'moduleScore.scarSignature' 
sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature <- sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature - min(sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("moduleScore.scarSignature"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("moduleScore.scarSignature"), split.by = "zt.scarType")
# construct an expressionMatrix for scarSignature as a single gene
sc.kNC_kJID_hNC.FB_0.16_test1.subset <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("SDC3"))
sc.kNC_kJID_hNC.FB_0.16_test1.subset@assays$RNA@data["SDC3",] <- sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature
# expression distribution by kernel density estimation
plot_density(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + xlim(-8, 8) + ylim(-10, 10)
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.scarType"
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16_test1.subset, idents = "HTS"), features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + labs(tag = "HTS", title = "scarSignature") + xlim(-8, 8) + ylim(-10, 10)
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16_test1.subset, idents = "Keloid"), features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + labs(tag = "Keloid", title = "scarSignature") + xlim(-8, 8) + ylim(-10, 10)
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16_test1.subset, idents = "Normal"), features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + labs(tag = "Normal", title = "scarSignature") + xlim(-8, 8) + ylim(-10, 10)


## htsSignature
htsSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/htsSignature.txt",
                     header = T, sep = "\t")
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# backup
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
# calculate ModuleScore of every cell
sc.kNC_kJID_hNC.FB_0.16_test1 <- AddModuleScore(sc.kNC_kJID_hNC.FB_0.16_test1, features = htsSig, name = "moduleScore.scarSignature")
colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data)[length(colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data))] <- 'moduleScore.scarSignature'
sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature <- sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature - min(sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("moduleScore.scarSignature"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("moduleScore.scarSignature"), split.by = "zt.scarType")
# construct an expressionMatrix for htsSignature as a single gene
sc.kNC_kJID_hNC.FB_0.16_test1.subset <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("SDC3"))
sc.kNC_kJID_hNC.FB_0.16_test1.subset@assays$RNA@data["SDC3",] <- sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature
# expression distribution by kernel density estimation
plot_density(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + xlim(-8, 8) + ylim(-10, 10)
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.scarType"
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16_test1.subset, idents = "HTS"), features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + labs(tag = "HTS", title = "htsSignature") + xlim(-8, 8) + ylim(-10, 10)
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16_test1.subset, idents = "Keloid"), features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + labs(tag = "Keloid", title = "htsSignature") + xlim(-8, 8) + ylim(-10, 10)
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16_test1.subset, idents = "Normal"), features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + labs(tag = "Normal", title = "htsSignature") + xlim(-8, 8) + ylim(-10, 10)


## keloidSignature
keloidSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/keloidSignature.txt",
                        header = T, sep = "\t")
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# backup
sc.kNC_kJID_hNC.FB_0.16_test1 <- sc.kNC_kJID_hNC.FB_0.16
# calculate ModuleScore of every cell
sc.kNC_kJID_hNC.FB_0.16_test1 <- AddModuleScore(sc.kNC_kJID_hNC.FB_0.16_test1, features = keloidSig, name = "moduleScore.scarSignature")
colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data)[length(colnames(sc.kNC_kJID_hNC.FB_0.16_test1@meta.data))] <- 'moduleScore.scarSignature'
sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature <- sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature - min(sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("moduleScore.scarSignature"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("moduleScore.scarSignature"), split.by = "zt.scarType")
# construct an expressionMatrix for scarSignature as a single gene
sc.kNC_kJID_hNC.FB_0.16_test1.subset <- subset(sc.kNC_kJID_hNC.FB_0.16_test1, features = c("SDC3"))
sc.kNC_kJID_hNC.FB_0.16_test1.subset@assays$RNA@data["SDC3",] <- sc.kNC_kJID_hNC.FB_0.16_test1$moduleScore.scarSignature
# expression distribution by kernel density estimation
plot_density(sc.kNC_kJID_hNC.FB_0.16_test1.subset, features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + xlim(-8, 8) + ylim(-10, 10)
Idents(sc.kNC_kJID_hNC.FB_0.16_test1.subset) <- "zt.scarType"
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16_test1.subset, idents = "HTS"), features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + labs(tag = "HTS", title = "keloidSignature") + xlim(-8, 8) + ylim(-10, 10)
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16_test1.subset, idents = "Keloid"), features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + labs(tag = "Keloid", title = "keloidSignature") + xlim(-8, 8) + ylim(-10, 10)
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16_test1.subset, idents = "Normal"), features = c("SDC3"), pal = "viridis", limits=c(0, 0.06)) + labs(tag = "Normal", title = "keloidSignature") + xlim(-8, 8) + ylim(-10, 10)


#################################################


####################### ScarsignatureGenes_FB6clusters: heatmap ##########################
## scarSignature: heatmap using ComplexHeatmap in FB6clusters of HTS
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16_HTS <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = "HTS")
# backup
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16_HTS
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap$integrated_snn_res.0.16)
# define features
scarSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/scarSignature.txt",
                      header = T, sep = "\t")
features <- unlist(scarSig)
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c", "#53bbc2", "#6b9cf8", "#e46fdc"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.exprMatrix, 
        cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_HTS_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


#################################################


####################### ScarsignatureGenes_FB: heatmap ##########################
## scarSignature: heatmap using ComplexHeatmap in FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# backup
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_forHeatmap))

## DoHeatmap
# define features
scarSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/scarSignature.txt",
                      header = T, sep = "\t")
features <- unlist(scarSig)
Idents(sc.kNC_kJID_hNC.FB_0.16_forHeatmap) <- "zt.scarType"
DoHeatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = features, slot = "scale.data", size = 3)

## ComplexHeatmap
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_forHeatmap$zt.scarType)
# define features
scarSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/scarSignature.txt",
                      header = T, sep = "\t")
features <- unlist(scarSig)
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix, col = colorRamp2(c(-2, 0, 2), c("grey", "white", "red")),
        cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


## htsSignature: heatmap using ComplexHeatmap in FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# backup
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_forHeatmap))

## DoHeatmap
# define features
htsSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/htsSignature.txt",
                      header = T, sep = "\t")
features <- unlist(htsSig)
Idents(sc.kNC_kJID_hNC.FB_0.16_forHeatmap) <- "zt.scarType"
DoHeatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = features, slot = "scale.data", size = 3)

## ComplexHeatmap
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_forHeatmap$zt.scarType)
# define features
htsSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/htsSignature.txt",
                      header = T, sep = "\t")
features <- unlist(htsSig)
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix, 
        cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


## keloidSignature: heatmap using ComplexHeatmap in FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# backup
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_forHeatmap))

## DoHeatmap
# define features
keloidSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/keloidSignature.txt",
                      header = T, sep = "\t")
features <- unlist(keloidSig)
Idents(sc.kNC_kJID_hNC.FB_0.16_forHeatmap) <- "zt.scarType"
DoHeatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, features = features, slot = "scale.data", size = 3)

## ComplexHeatmap
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_forHeatmap$zt.scarType)
# define features
keloidSig <- read.table(file = "R/sc.KNC_kJID_hNC/FB/SignatureGenes-scar_keloid_HTS/newSignature_20220401/keloidSignature.txt",
                     header = T, sep = "\t")
features <- unlist(keloidSig)
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


#################################################


####################### FB6clusters: PCA ##########################
## Combine FBCells into 3scarTypes
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
# FB.HTS
sc.kNC_kJID_hNC.FB_0.16.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'HTS')
sc.kNC_kJID_hNC.FB_0.16.HTS.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.HTS@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.HTS.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.HTS.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.HTS.exprMatr_sum) <- 'HTS'
# FB.Keloid
sc.kNC_kJID_hNC.FB_0.16.Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'Keloid')
sc.kNC_kJID_hNC.FB_0.16.Keloid.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.Keloid@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.Keloid.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.Keloid.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.Keloid.exprMatr_sum) <- 'Keloid'
# FB.Normal
sc.kNC_kJID_hNC.FB_0.16.Normal <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'Normal')
sc.kNC_kJID_hNC.FB_0.16.Normal.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.Normal@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.Normal.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.Normal.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.Normal.exprMatr_sum) <- 'Normal'
# put 3scarTypes in one table
sc.kNC_kJID_hNC.FB_0.16.3scarTypes.exprMatr_sum <- cbind(sc.kNC_kJID_hNC.FB_0.16.HTS.exprMatr_sum, sc.kNC_kJID_hNC.FB_0.16.Keloid.exprMatr_sum, sc.kNC_kJID_hNC.FB_0.16.Normal.exprMatr_sum)
write.csv(sc.kNC_kJID_hNC.FB_0.16.3scarTypes.exprMatr_sum, file = 'R/sc.KNC_kJID_hNC/FB/pseudoBulkPCA/exprMatr_sum.Splitby3scarTypes.csv')


## Combine FBCells into 3dataSources
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.dataSource"
# FB.sc.hNC
sc.kNC_kJID_hNC.FB_0.16.sc.hNC <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'sc.hNC')
sc.kNC_kJID_hNC.FB_0.16.sc.hNC.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.sc.hNC@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.sc.hNC.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.sc.hNC.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.sc.hNC.exprMatr_sum) <- 'sc.hNC'
# FB.sc.kJID
sc.kNC_kJID_hNC.FB_0.16.sc.kJID <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'sc.kJID')
sc.kNC_kJID_hNC.FB_0.16.sc.kJID.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.sc.kJID@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.sc.kJID.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.sc.kJID.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.sc.kJID.exprMatr_sum) <- 'sc.kJID'
# FB.sc.kNC
sc.kNC_kJID_hNC.FB_0.16.sc.kNC <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'sc.kNC')
sc.kNC_kJID_hNC.FB_0.16.sc.kNC.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.sc.kNC@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.sc.kNC.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.sc.kNC.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.sc.kNC.exprMatr_sum) <- 'sc.kNC'
# put 3dataSources in one table
sc.kNC_kJID_hNC.FB_0.16.3dataSources.exprMatr_sum <- cbind(sc.kNC_kJID_hNC.FB_0.16.sc.hNC.exprMatr_sum, sc.kNC_kJID_hNC.FB_0.16.sc.kJID.exprMatr_sum, sc.kNC_kJID_hNC.FB_0.16.sc.kNC.exprMatr_sum)
write.csv(sc.kNC_kJID_hNC.FB_0.16.3dataSources.exprMatr_sum, file = 'R/sc.KNC_kJID_hNC/FB/pseudoBulkPCA/exprMatr_sum.Splitby3dataSources.csv')


## Combine FBCells into 20sampleNames
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "orig.ident"
# FB.HTS1
sc.kNC_kJID_hNC.FB_0.16.HTS1 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'HTS1')
sc.kNC_kJID_hNC.FB_0.16.HTS1.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.HTS1@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.HTS1.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.HTS1.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.HTS1.exprMatr_sum) <- 'HTS1'
# FB.HTS2
sc.kNC_kJID_hNC.FB_0.16.HTS2 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'HTS2')
sc.kNC_kJID_hNC.FB_0.16.HTS2.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.HTS2@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.HTS2.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.HTS2.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.HTS2.exprMatr_sum) <- 'HTS2'
# FB.HTS3
sc.kNC_kJID_hNC.FB_0.16.HTS3 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'HTS3')
sc.kNC_kJID_hNC.FB_0.16.HTS3.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.HTS3@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.HTS3.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.HTS3.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.HTS3.exprMatr_sum) <- 'HTS3'
# FB.K007CASE
sc.kNC_kJID_hNC.FB_0.16.K007CASE <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'K007CASE')
sc.kNC_kJID_hNC.FB_0.16.K007CASE.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.K007CASE@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.K007CASE.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.K007CASE.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.K007CASE.exprMatr_sum) <- 'K007CASE'
# FB.K007CTRL
sc.kNC_kJID_hNC.FB_0.16.K007CTRL <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'K007CTRL')
sc.kNC_kJID_hNC.FB_0.16.K007CTRL.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.K007CTRL@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.K007CTRL.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.K007CTRL.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.K007CTRL.exprMatr_sum) <- 'K007CTRL'
# FB.K009CASE
sc.kNC_kJID_hNC.FB_0.16.K009CASE <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'K009CASE')
sc.kNC_kJID_hNC.FB_0.16.K009CASE.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.K009CASE@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.K009CASE.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.K009CASE.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.K009CASE.exprMatr_sum) <- 'K009CASE'
# FB.K009CTRL
sc.kNC_kJID_hNC.FB_0.16.K009CTRL <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'K009CTRL')
sc.kNC_kJID_hNC.FB_0.16.K009CTRL.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.K009CTRL@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.K009CTRL.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.K009CTRL.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.K009CTRL.exprMatr_sum) <- 'K009CTRL'
# FB.K012CASE
sc.kNC_kJID_hNC.FB_0.16.K012CASE <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'K012CASE')
sc.kNC_kJID_hNC.FB_0.16.K012CASE.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.K012CASE@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.K012CASE.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.K012CASE.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.K012CASE.exprMatr_sum) <- 'K012CASE'
# FB.K012CTRL
sc.kNC_kJID_hNC.FB_0.16.K012CTRL <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'K012CTRL')
sc.kNC_kJID_hNC.FB_0.16.K012CTRL.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.K012CTRL@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.K012CTRL.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.K012CTRL.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.K012CTRL.exprMatr_sum) <- 'K012CTRL'
# FB.K013CASE
sc.kNC_kJID_hNC.FB_0.16.K013CASE <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'K013CASE')
sc.kNC_kJID_hNC.FB_0.16.K013CASE.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.K013CASE@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.K013CASE.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.K013CASE.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.K013CASE.exprMatr_sum) <- 'K013CASE'
# FB.K013CTRL
sc.kNC_kJID_hNC.FB_0.16.K013CTRL <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'K013CTRL')
sc.kNC_kJID_hNC.FB_0.16.K013CTRL.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.K013CTRL@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.K013CTRL.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.K013CTRL.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.K013CTRL.exprMatr_sum) <- 'K013CTRL'
# FB.KF1
sc.kNC_kJID_hNC.FB_0.16.KF1 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'KF1')
sc.kNC_kJID_hNC.FB_0.16.KF1.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.KF1@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.KF1.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.KF1.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.KF1.exprMatr_sum) <- 'KF1'
# FB.KF2
sc.kNC_kJID_hNC.FB_0.16.KF2 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'KF2')
sc.kNC_kJID_hNC.FB_0.16.KF2.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.KF2@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.KF2.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.KF2.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.KF2.exprMatr_sum) <- 'KF2'
# FB.KF3
sc.kNC_kJID_hNC.FB_0.16.KF3 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'KF3')
sc.kNC_kJID_hNC.FB_0.16.KF3.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.KF3@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.KF3.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.KF3.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.KF3.exprMatr_sum) <- 'KF3'
# FB.NF1
sc.kNC_kJID_hNC.FB_0.16.NF1 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'NF1')
sc.kNC_kJID_hNC.FB_0.16.NF1.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.NF1@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.NF1.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.NF1.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.NF1.exprMatr_sum) <- 'NF1'
# FB.NF2
sc.kNC_kJID_hNC.FB_0.16.NF2 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'NF2')
sc.kNC_kJID_hNC.FB_0.16.NF2.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.NF2@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.NF2.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.NF2.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.NF2.exprMatr_sum) <- 'NF2'
# FB.NF3
sc.kNC_kJID_hNC.FB_0.16.NF3 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'NF3')
sc.kNC_kJID_hNC.FB_0.16.NF3.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.NF3@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.NF3.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.NF3.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.NF3.exprMatr_sum) <- 'NF3'
# FB.Normal1
sc.kNC_kJID_hNC.FB_0.16.Normal1 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'Normal1')
sc.kNC_kJID_hNC.FB_0.16.Normal1.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.Normal1@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.Normal1.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.Normal1.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.Normal1.exprMatr_sum) <- 'Normal1'
# FB.Normal2
sc.kNC_kJID_hNC.FB_0.16.Normal2 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'Normal2')
sc.kNC_kJID_hNC.FB_0.16.Normal2.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.Normal2@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.Normal2.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.Normal2.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.Normal2.exprMatr_sum) <- 'Normal2'
# FB.Normal3
sc.kNC_kJID_hNC.FB_0.16.Normal3 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'Normal3')
sc.kNC_kJID_hNC.FB_0.16.Normal3.exprMatr <- sc.kNC_kJID_hNC.FB_0.16.Normal3@assays$RNA@counts
sc.kNC_kJID_hNC.FB_0.16.Normal3.exprMatr_sum <- apply(sc.kNC_kJID_hNC.FB_0.16.Normal3.exprMatr, 1, sum) %>% data.frame()
colnames(sc.kNC_kJID_hNC.FB_0.16.Normal3.exprMatr_sum) <- 'Normal3'
# put 20sampleNames in one table
sc.kNC_kJID_hNC.FB_0.16.20sampleNames.exprMatr_sum <- cbind(sc.kNC_kJID_hNC.FB_0.16.HTS1.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.HTS2.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.HTS3.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.K007CASE.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.K007CTRL.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.K009CASE.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.K009CTRL.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.K012CASE.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.K012CTRL.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.K013CASE.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.K013CTRL.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.KF1.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.KF2.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.KF3.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.NF1.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.NF2.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.NF3.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.Normal1.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.Normal2.exprMatr_sum,
                                                            sc.kNC_kJID_hNC.FB_0.16.Normal3.exprMatr_sum)
write.csv(sc.kNC_kJID_hNC.FB_0.16.20sampleNames.exprMatr_sum, file = 'R/sc.KNC_kJID_hNC/FB/pseudoBulkPCA/exprMatr_sum.Splitby20sampleNames.csv')


#################################################


####################### heatmap: FB_6 clusters ##########################
## union of top50 DEGs: heatmap using ComplexHeatmap in FBc0
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
sc.kNC_kJID_hNC.FB_0.16_c0 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = "0")
Idents(sc.kNC_kJID_hNC.FB_0.16_c0) <- "zt.scarType"
# backup
sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16_c0
sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap$zt.scarType)
# define features
features <- union.top50_down50
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap.exprMatrix, 
        cluster_rows = F, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_c0_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

## union of top50 DEGs: heatmap using ComplexHeatmap in FBc1
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
sc.kNC_kJID_hNC.FB_0.16_c1 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = "1")
Idents(sc.kNC_kJID_hNC.FB_0.16_c1) <- "zt.scarType"
# backup
sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16_c1
sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap$zt.scarType)
# define features
features <- union.top50_down50
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap.exprMatrix, 
        cluster_rows = F, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_c1_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

## union of top50 DEGs: heatmap using ComplexHeatmap in FBc2
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
sc.kNC_kJID_hNC.FB_0.16_c2 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = "2")
Idents(sc.kNC_kJID_hNC.FB_0.16_c2) <- "zt.scarType"
# backup
sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16_c2
sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap$zt.scarType)
# define features
features <- union.top50_down50
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap.exprMatrix, 
        cluster_rows = F, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_c2_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

## union of top50 DEGs: heatmap using ComplexHeatmap in FBc3
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
sc.kNC_kJID_hNC.FB_0.16_c3 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = "3")
Idents(sc.kNC_kJID_hNC.FB_0.16_c3) <- "zt.scarType"
# backup
sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16_c3
sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap$zt.scarType)
# define features
features <- union.top50_down50
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap.exprMatrix, 
        cluster_rows = F, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_c3_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

## union of top50 DEGs: heatmap using ComplexHeatmap in FBc4
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
sc.kNC_kJID_hNC.FB_0.16_c4 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = "4")
Idents(sc.kNC_kJID_hNC.FB_0.16_c4) <- "zt.scarType"
# backup
sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16_c4
sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap$zt.scarType)
# define features
features <- union.top50_down50
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap.exprMatrix, 
        cluster_rows = F, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_c4_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

## union of top50 DEGs: heatmap using ComplexHeatmap in FBc5
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
sc.kNC_kJID_hNC.FB_0.16_c5 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = "5")
Idents(sc.kNC_kJID_hNC.FB_0.16_c5) <- "zt.scarType"
# backup
sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap <- sc.kNC_kJID_hNC.FB_0.16_c5
sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap))
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap, slot = "scale.data")
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap$zt.scarType)
# define features
features <- union.top50_down50
# talor the expression matrix for heatmap (not all genes, but all pre-defined features)
sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap.exprMatrix <- as.matrix(sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap.exprMatrix[features, names(sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap.cluster_info)])
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e77d71", "#b49f33", "#64b74c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap.exprMatrix, 
        cluster_rows = F, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_c5_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


#################################################


####################### FB(DEG.union of HTSandKeloidvsNormal).trajectoryAnalysis:monocle2 ##########################
# 
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# subtract raw data to establish a cds object
FB.counts <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16, assay = 'RNA', slot = 'counts')
FB.cell_metadata <- sc.kNC_kJID_hNC.FB_0.16@meta.data
FB.gene_annotation <- data.frame(gene_short_name = rownames(FB.counts), row.names = rownames(FB.counts))
pd <- new('AnnotatedDataFrame', data = FB.cell_metadata) 
fd <- new('AnnotatedDataFrame', data = FB.gene_annotation)
cds <- newCellDataSet(FB.counts,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Gene signature for ordering: the union of DEGs in HTSvsNormal and KeloidvsNormal
# FB_HTSvsNormal: DEGs_p<0.05
DEG.FB_HTSvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_HTSvsNormal.xlsx")
DEG.FB_HTSvsNormal_p0.05 <- subset(DEG.FB_HTSvsNormal, p_val_adj < 0.05)
# FB_KeloidvsNormal: DEGs_p<0.05
DEG.FB_KeloidvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_KeloidvsNormal.xlsx")
DEG.FB_KeloidvsNormal_p0.05 <- subset(DEG.FB_KeloidvsNormal, p_val_adj < 0.05)
# all DEGs: the union of DEG.FB_HTSvsNormal_p0.05 and DEG.FB_KeloidvsNormal_p0.05
union_DEG.FB_HTSandKeloidvsNormal_p0.05 <- union(rownames(DEG.FB_HTSvsNormal_p0.05), rownames(DEG.FB_KeloidvsNormal_p0.05))

ordergene <- union_DEG.FB_HTSandKeloidvsNormal_p0.05
cds <- setOrderingFilter(cds, ordergene)

#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
plot_ordering_genes(cds)
#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

cds <- orderCells(cds)
#⚠️使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
#cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点
save(cds, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/Union/cds.Rdata")

# check the trajectory by Pseudotime
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)
# check the trajectory by FB_6clusters
plot_cell_trajectory(cds,color_by="integrated_snn_res.0.16", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "integrated_snn_res.0.16") + facet_wrap("~integrated_snn_res.0.16", nrow = 1)
# check the trajectory by HTS/Keloid/Normal
plot_cell_trajectory(cds,color_by="zt.scarType", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "zt.scarType") + facet_wrap("~zt.scarType", nrow = 1)
# check the trajectory by dataSource
plot_cell_trajectory(cds,color_by="zt.dataSource", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "zt.dataSource") + facet_wrap("~zt.dataSource", nrow = 1)
# check the trajectory by Source
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)

## check the expression of a certain gene
cds_subset <- cds[c("COL1A1", "ACTA2"),]
# the distribution of gene expression in different State/Pseudotime
plot_genes_in_pseudotime(cds_subset, color_by = "State")
plot_genes_in_pseudotime(cds_subset, color_by = "zt.scarType")
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
# the distribution of gene expression in the trajectory tree
pData(cds)$COL1A1 = log2(exprs(cds)['COL1A1',]+1)
plot_cell_trajectory(cds, color_by = "COL1A1", cell_size=1) + scale_alpha(range = c(0, 1))
pData(cds)$ACTA2 = log2(exprs(cds)['ACTA2',]+1)
plot_cell_trajectory(cds, color_by = "ACTA2", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$PCNA = log2(exprs(cds)['PCNA',]+1)
plot_cell_trajectory(cds, color_by = "PCNA", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$MKI67 = log2(exprs(cds)['MKI67',]+1)
plot_cell_trajectory(cds, color_by = "MKI67", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red') + scale_alpha(range = c(0, 1))


#################################################


####################### FB(DEG.intersection of HTSandKeloidvsNormal).trajectoryAnalysis:monocle2 ##########################
# 
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# subtract raw data to establish a cds object
FB.counts <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16, assay = 'RNA', slot = 'counts')
FB.cell_metadata <- sc.kNC_kJID_hNC.FB_0.16@meta.data
FB.gene_annotation <- data.frame(gene_short_name = rownames(FB.counts), row.names = rownames(FB.counts))
pd <- new('AnnotatedDataFrame', data = FB.cell_metadata) 
fd <- new('AnnotatedDataFrame', data = FB.gene_annotation)
cds <- newCellDataSet(FB.counts,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Gene signature for ordering: the union of DEGs in HTSvsNormal and KeloidvsNormal
# FB_HTSvsNormal: DEGs_p<0.05
DEG.FB_HTSvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_HTSvsNormal.xlsx")
DEG.FB_HTSvsNormal_p0.05 <- subset(DEG.FB_HTSvsNormal, p_val_adj < 0.05)
# FB_KeloidvsNormal: DEGs_p<0.05
DEG.FB_KeloidvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_KeloidvsNormal.xlsx")
DEG.FB_KeloidvsNormal_p0.05 <- subset(DEG.FB_KeloidvsNormal, p_val_adj < 0.05)
# all DEGs: the union of DEG.FB_HTSvsNormal_p0.05 and DEG.FB_KeloidvsNormal_p0.05
intersect_DEG.FB_HTSandKeloidvsNormal_p0.05 <- intersect(rownames(DEG.FB_HTSvsNormal_p0.05), rownames(DEG.FB_KeloidvsNormal_p0.05))

ordergene <- intersect_DEG.FB_HTSandKeloidvsNormal_p0.05
load( file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/cds.Rdata")
cds <- setOrderingFilter(cds, ordergene)

#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
plot_ordering_genes(cds)
#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

cds <- orderCells(cds)
#⚠️使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
#cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点
save(cds, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/Intersection/cds.Rdata")

load(file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/Intersection/cds.Rdata")
# check the trajectory by Pseudotime
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)
# check the trajectory by FB_6clusters
plot_cell_trajectory(cds,color_by="integrated_snn_res.0.16", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "integrated_snn_res.0.16") + facet_wrap("~integrated_snn_res.0.16", nrow = 1)
# check the trajectory by HTS/Keloid/Normal
plot_cell_trajectory(cds,color_by="zt.scarType", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "zt.scarType") + facet_wrap("~zt.scarType", nrow = 1)
# check the trajectory by dataSource
plot_cell_trajectory(cds,color_by="zt.dataSource", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "zt.dataSource") + facet_wrap("~zt.dataSource", nrow = 1)
# check the trajectory by Source
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)

## check the expression of a certain gene
load("R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/unionofAllDEGofHTSandKeloid/cds.ordered.Rdata")
cds_subset <- cds[c("COL1A1", "ACTA2"),]
# the distribution of gene expression in different State/Pseudotime
plot_genes_in_pseudotime(cds_subset, color_by = "State")
plot_genes_in_pseudotime(cds_subset, color_by = "zt.scarType")
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
# the distribution of gene expression in the trajectory tree
pData(cds)$COL1A1 = log2(exprs(cds)['COL1A1',]+1)
plot_cell_trajectory(cds, color_by = "COL1A1", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$ACTA2 = log2(exprs(cds)['ACTA2',]+1)
plot_cell_trajectory(cds, color_by = "ACTA2", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$PCNA = log2(exprs(cds)['PCNA',]+1)
plot_cell_trajectory(cds, color_by = "PCNA", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$MKI67 = log2(exprs(cds)['MKI67',]+1)
plot_cell_trajectory(cds, color_by = "MKI67", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red') + scale_alpha(range = c(0, 1))
pData(cds)$POSTN = log2(exprs(cds)['POSTN',]+1)
plot_cell_trajectory(cds, color_by = "POSTN", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red') + scale_alpha(range = c(0, 1))


#################################################


####################### FB_6clusters(KeloidvsNormal.DEG).trajectoryAnalysis:monocle2 ##########################
# 
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# subtract raw data to establish a cds object
FB.counts <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16, assay = 'RNA', slot = 'counts')
FB.cell_metadata <- sc.kNC_kJID_hNC.FB_0.16@meta.data
FB.gene_annotation <- data.frame(gene_short_name = rownames(FB.counts), row.names = rownames(FB.counts))
pd <- new('AnnotatedDataFrame', data = FB.cell_metadata) 
fd <- new('AnnotatedDataFrame', data = FB.gene_annotation)
cds <- newCellDataSet(FB.counts,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

diff <- FindVariableFeatures(sc.kNC_kJID_hNC.FB_0.16, selection.method = "vst", nfeatures = 2000)
diff <- diff@assays$integrated@var.features
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)

# Gene signature for ordering: the union of DEGs in HTSvsNormal and KeloidvsNormal
# FB_KeloidvsNormal: DEGs_p<0.05
DEG.FB_KeloidvsNormal <- read.table("R/sc.KNC_kJID_hNC/FB/DEG.FB_KeloidvsNormal.xlsx")
DEG.FB_KeloidvsNormal_p0.05 <- subset(DEG.FB_KeloidvsNormal, p_val_adj < 0.05)

deg <- DEG.FB_KeloidvsNormal_p0.05
head(deg)

## 轨迹构建基因可视化
load( file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/cds.Rdata")
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
#这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
plot_ordering_genes(cds)
#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

cds <- orderCells(cds)
#⚠️使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
#cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点
save(cds, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/HVG2000/cds.Rdata")

# check the trajectory by Pseudotime
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)
# check the trajectory by FB_6clusters
plot_cell_trajectory(cds,color_by="integrated_snn_res.0.16", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "integrated_snn_res.0.16") + facet_wrap("~integrated_snn_res.0.16", nrow = 1)
# check the trajectory by HTS/Keloid/Normal
plot_cell_trajectory(cds,color_by="zt.scarType", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "zt.scarType") + facet_wrap("~zt.scarType", nrow = 1)
# check the trajectory by dataSource
plot_cell_trajectory(cds,color_by="zt.dataSource", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "zt.dataSource") + facet_wrap("~zt.dataSource", nrow = 1)
# check the trajectory by Source
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)

## check the expression of a certain gene
load("R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/HVG2000/cds.Rdata")
cds_subset <- cds[c("COL1A1", "ACTA2"),]
# the distribution of gene expression in different State/Pseudotime
plot_genes_in_pseudotime(cds_subset, color_by = "State")
plot_genes_in_pseudotime(cds_subset, color_by = "zt.scarType")
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
# the distribution of gene expression in the trajectory tree
pData(cds)$COL1A1 = log2(exprs(cds)['COL1A1',]+1)
plot_cell_trajectory(cds, color_by = "COL1A1", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$ACTA2 = log2(exprs(cds)['ACTA2',]+1)
plot_cell_trajectory(cds, color_by = "ACTA2", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$PCNA = log2(exprs(cds)['PCNA',]+1)
plot_cell_trajectory(cds, color_by = "PCNA", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$MKI67 = log2(exprs(cds)['MKI67',]+1)
plot_cell_trajectory(cds, color_by = "MKI67", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red') + scale_alpha(range = c(0, 1))


#################################################


####################### FB_6clusters(dFeature.DEGs).trajectoryAnalysis:monocle2 ##########################
# 
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# subtract raw data to establish a cds object
FB.counts <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16, assay = 'RNA', slot = 'counts')
FB.cell_metadata <- sc.kNC_kJID_hNC.FB_0.16@meta.data
FB.gene_annotation <- data.frame(gene_short_name = rownames(FB.counts), row.names = rownames(FB.counts))
pd <- new('AnnotatedDataFrame', data = FB.cell_metadata) 
fd <- new('AnnotatedDataFrame', data = FB.gene_annotation)
cds <- newCellDataSet(FB.counts,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

diff <- differentialGeneTest(cds,fullModelFormulaStr="~integrated_snn_res.0.16",cores=1) 
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)
save(diff, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/allDEGinFB6cluster_monocle2.Rdata")
save(cds, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/cds.Rdata")

load(file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/allDEGinFB6cluster_monocle2.Rdata")
##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
deg <- subset(diff, pval < 0.05 & qval < 0.01) 
deg <- deg[order(deg$qval,decreasing=F),]
deg <- head(deg, 2000) #选出top2000个基因
head(deg)


## 轨迹构建基因可视化
load(file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/cds.Rdata")
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)  
#这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
plot_ordering_genes(cds)
#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

cds <- orderCells(cds)
#⚠️使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
#cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点
cds <- orderCells(cds, root_state = 4)
save(cds, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/cds.Rdata")

load(file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/cds.Rdata")
# check the trajectory by Pseudotime
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)  + facet_wrap("~zt.scarType", nrow = 1)
# check the trajectory by FB_6clusters
plot_cell_trajectory(cds,color_by="integrated_snn_res.0.16", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "integrated_snn_res.0.16") + facet_wrap("~integrated_snn_res.0.16", nrow = 1)
# check the trajectory by HTS/Keloid/Normal
plot_cell_trajectory(cds,color_by="zt.scarType", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "zt.scarType") + facet_wrap("~zt.scarType", nrow = 1)
# check the trajectory by dataSource
plot_cell_trajectory(cds,color_by="zt.dataSource", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "zt.dataSource") + facet_wrap("~zt.dataSource", nrow = 1)
# check the trajectory by Source
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
# add a new column to the metadata, in order to calculate cellNumbers in each state
cds$state_scarType <- paste(cds$State, cds$zt.scarType, sep = "_")
cellNumbers <- as.data.frame(table(cds$state_scarType))
write.table(cellNumbers, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/cellNumbersOfStates.xlsx", sep = "\t")
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
cds.State <- data.frame(cds$State, row.names = cds@assayData$exprs@Dimnames[[2]])
cds.State$cds.State <- gsub('4', 'preBranch', cds.State$cds.State)
cds.State$cds.State <- gsub('1|2|5', 'cellFate1', cds.State$cds.State)
cds.State$cds.State <- gsub('3', 'cellFate2', cds.State$cds.State)
sc.kNC_kJID_hNC.FB_0.16$dFeatureState <- cds.State
save(sc.kNC_kJID_hNC.FB_0.16, file = "R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
# add a new column to the metadata, in order to calculate cellNumbers of every cluster in each state
cds$state_scarType_FB6cluster <- paste(cds$State, cds$zt.scarType, cds$integrated_snn_res.0.16,sep = "_")
cellNumbers <- as.data.frame(table(cds$state_scarType_FB6cluster))
colnames(cellNumbers) <- c("state_scarType_FB6cluster", "Freq")
write.table(cellNumbers, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/cellNumbersOfFB6clustersInStates.xlsx", sep = "\t")


## check the expression of a certain gene
load("R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/cds.Rdata")
cds_subset <- cds[c("COL1A1", "ACTA2"),]
# the distribution of gene expression in different State/Pseudotime
plot_genes_in_pseudotime(cds_subset, color_by = "State")
plot_genes_in_pseudotime(cds_subset, color_by = "zt.scarType")
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
# the distribution of gene expression in the trajectory tree
pData(cds)$COL1A1 = log2(exprs(cds)['COL1A1',]+1)
plot_cell_trajectory(cds, color_by = "COL1A1", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$ACTA2 = log2(exprs(cds)['ACTA2',]+1)
plot_cell_trajectory(cds, color_by = "ACTA2", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$PCNA = log2(exprs(cds)['PCNA',]+1)
plot_cell_trajectory(cds, color_by = "PCNA", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$MKI67 = log2(exprs(cds)['MKI67',]+1)
plot_cell_trajectory(cds, color_by = "MKI67", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red') + scale_alpha(range = c(0, 1))
pData(cds)$POSTN = log2(exprs(cds)['POSTN',]+1)
plot_cell_trajectory(cds, color_by = "POSTN", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red') + scale_alpha(range = c(0, 1))
pData(cds)$STAT3 = log2(exprs(cds)['STAT3',]+1)
plot_cell_trajectory(cds, color_by = "STAT3", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red') + scale_alpha(range = c(0, 1))
pData(cds)$MYC = log2(exprs(cds)['MYC',]+1)
plot_cell_trajectory(cds, color_by = "MYC", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red') + scale_alpha(range = c(0, 1))


## plot_genes_branched_heatmap
BEAM_res <- BEAM(cds[ordergene,], branch_point = 1, cores = 2) 
save(BEAM_res, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/plot_genes_branched_heatmap.BEAM_res.Rdata")
load("R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/plot_genes_branched_heatmap.BEAM_res.Rdata")
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #对2829个基因进行排序，运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_genes <- row.names(subset(BEAM_res, qval < 1e-4))
p <- plot_genes_branched_heatmap(cds[BEAM_genes,],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T,
                            return_heatmap = T)
# extract allDEGs in the heatmap
hp.genes <- p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- BEAM_res[hp.genes, c("gene_short_name", "pval", "qval")]
write.table(BEAM_sig, "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/FB.branched_heatmap.allDEGs.xlsx", sep = "\t", row.names = F)
# top100
BEAM_genes.c1.top25 <- c("MMP11", "GXYLT2", "UQCR11", "GPX7", "S100A6", "CREB3L1", "MAP2K2", "HSPB3", "WFDC1", "HSPB6", "IGF1", "CPE", "PMEPA1", "GPX4", "TCEAL9", "CDC42EP3", "PLPP5", "SERTAD4-AS1", "CD109", "MMP16", "NRN1", "AK5", "PDGFD", "RCAN1", "RBP4")
BEAM_genes.c2.top25 <- c("IER5", "BTG3", "PLEKHH2", "TPT1", "BTF3", "HLA-F", "PLSCR1", "PFKFB3", "RPL6", "PGF", "CBX4", "ZC3HAV1", "CSRNP1", "HES4", "CSF1", "KMT2E", "AL118516.1", "RPS5", "DOK6", "WNT5A", "CFI", "CSGALNACT1", "HIC1", "RPL28", "SMIM3")
BEAM_genes.c3.top25 <- c("ABCA10", "EVA1B", "IL24", "THSD4", "PLXDC1", "PAPPA", "GDF10", "CMKLR1", "UNC5B", "DAB2", "ITM2B", "FRMD6", "PODNL1", "COL7A1", "SEPT9", "CISD1", "SBSPON", "TMEM52", "STC2", "SDC2", "TPBG", "APELA", "EBF2", "PHLDA3", "MMP1")
BEAM_genes.c4.top25 <- c("SEMA3D", "ARHGAP29", "DSTN", "RIN2", "TYROBP", "PTPRE", "SELE", "RASGEF1B", "RNASE1", "ECSCR", "CCL3", "ADAMTS9", "MATN2", "CCDC3", "IL1B", "TOMM7", "DAAM1", "LYZ", "COL13A1", "MCTP1", "MTSS1", "NFATC2", "C1orf56", "CCL3L1", "CRISPLD2")
BEAM_genes <- c(BEAM_genes.c1.top25, BEAM_genes.c2.top25, BEAM_genes.c3.top25, BEAM_genes.c4.top25)
BEAM_genes <- BEAM_genes.c1.top25
p.100 <- plot_genes_branched_heatmap(cds[BEAM_genes,],
                                 branch_point = 1, #绘制的是哪个分支
                                 num_clusters = 4, #分成几个cluster，根据需要调整
                                 cores = 1,
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T)

# selectedGenes (important TFs)
selectedGenes <- c('COL1A1', 'TWIST1', 'ALX3', 'USF2', 'HOXC5', 'ATF2', 'ZNF76', 'ZFX', 'NRF1', 'FOXM1', 'HOXA3', 'CRTC2', 'CARF', 'ZNF319', 'MAX', 'NPAS2', 'FOXA1', 'WT1', 'FOXL1', 'ZNF562', 'ATF1', 'ZBTB33', 'TFF3', 'FOXD3', 'DLX5', 'ZNF823', 'USF1', 'UBE2V1', 'SOX12', 'BARX2', 'POU2F2', 'ZNF12', 'ZNF689', 'RFX3', 'ATF6B', 'ZNF639', 'TBP', 'VDR', 'HOXA4', 'ARID3A')
p.selectedGenes <- plot_genes_branched_heatmap(cds[selectedGenes,],
                                     branch_point = 1, #绘制的是哪个分支
                                     num_clusters = 4, #分成几个cluster，根据需要调整
                                     cores = 1,
                                     use_gene_short_name = T,
                                     show_rownames = T,
                                     return_heatmap = T)
# extract genes in the heatmap
hp.genes <- p.100$ph_res$tree_row$labels[p.100$ph_res$tree_row$order]
BEAM_sig <- BEAM_res[hp.genes, c("gene_short_name", "pval", "qval")]
write.table(BEAM_sig, "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/FB.branched_heatmap.genes.xlsx", sep = "\t")
# plot_genes_branched_heatmap of TF
TF <- c("FOXP1", "ZNF469", "SCX", "CREB5", "TWIST2", "ZFHX4", "KLF6", "MXD4", "BNC2", "CARHSP1", "XBP1", "CREB3L1", "CEBPD", "NR4A1", "ZFP36L2", "ZFP36L1", "HES1", "FOS", "SOX4", "STAT3", "YBX3", "MYC", "NFE2L2", "NFIA", "ARID5B", "REL", "MSX1", "TWIST1", "ATF3", "IRF8", "MAFF", "ZNF331", "IRF1", "KLF5", "KLF9", "ID1", "ID4", "JUNB", "JUN", "SNAI1", "NR4A2", "FOSB", "ID3", "NFIL3", "NR2F2", "BCL6", "ID2", "CEBPB", "ETS2", "EGR3", "NFKB1", "EGR1", "HBP1", "NR4A3", "MAFB", "BHLHE40", "THAP2", "TBX18", "SFPQ", "NFKB2", "KLF10", "PRRX2", "JUND", "IRX1", "ZBTB7A", "ZEB2", "EBF1", "GPBP1", "TSC22D1", "LITAF", "FOXO3", "AHR", "ERF", "EPAS1", "ARID4B", "MAFG", "KLF2", "ETS1", "SNAI2", "RELB", "IRX2", "CREM", "FOXD1", "MSC", "TCF7L1", "CREBRF", "CNBP", "JARID2", "ZNF267", "CEBPA", "NR3C1", "PRDM1", "POU3F1", "DDIT3", "KLF13", "BHLHE22", "YY1", "GATAD2B", "MEIS2", "FOXC1", "ARID5A", "ZBTB1", "TFAM", "TCF7L2", "SON", "TBX3", "PRDM8", "RUNX3", "ZBTB20", "HES5", "EGR2", "BBX", "CEBPZ", "PGR", "TRPS1", "ELF1", "PPARG", "PURA", "TSHZ2", "RUNX1", "TSHZ3", "FOXC2", "ZFHX3", "HOPX", "TFAP2A", "ETV5", "HIC1", "HES4", "CSRNP1", "PLSCR1", "TCF4", "MEF2C", "ALX4", "HMGB1", "SOX18", "MEOX2", "GLI1", "PRRX1", "ETV1", "SOX11", "FOXS1", "EBF2", "TSC22D3", "NFIB", "FOSL2", "TBX15", "ELMSAN1", "TFAP2C", "MAF", "MEF2A", "LRRFIP1", "NFATC2")
BEAM_genes <- TF
p <- plot_genes_branched_heatmap(cds[BEAM_genes,],
                                 branch_point = 1, #绘制的是哪个分支
                                 num_clusters = 4, #分成几个cluster，根据需要调整
                                 cores = 1,
                                 use_gene_short_name = T,
                                 show_rownames = T,
                                 return_heatmap = T)
# plot_genes_branched_heatmap of TF, which is the intersection of 10 iterations of GeneRegulatoryNetwork
TF <- c('MEF2C', 'CEBPD', 'FOSL2', 'CREB3L1', 'FOS', 'TBX15', 'MXD4', 'TFAP2A', 'CREB5', 'EGR2', 'XBP1', 'STAT3', 'PRDM1', 'FOXO3', 'EGR1', 'MYC', 'NFE2L2', 'ETS2', 'CREM', 'REL', 'NFKB1', 'MAFG', 'ELF1', 'RELB', 'YY1', 'ATF3', 'IRF8', 'ZBTB7A', 'MAFF', 'ETS1', 'IRF1', 'KLF5', 'BHLHE40', 'JUNB', 'JUN', 'JUND', 'FOSB', 'NFIL3', 'FOXD1', 'NR2F2', 'EGR3', 'NFKB2', 'CEBPB')
BEAM_genes <- TF
p <- plot_genes_branched_heatmap(cds[BEAM_genes,],
                                 branch_point = 1, #绘制的是哪个分支
                                 num_clusters = 4, #分成几个cluster，根据需要调整
                                 cores = 1,
                                 use_gene_short_name = T,
                                 show_rownames = T,
                                 return_heatmap = T)


## plot_pseudotime_heatmap
# 这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, 
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
save(Time_diff, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/plot_pseudotime_heatmap.Time_diff.Rdata")
load(file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/plot_pseudotime_heatmap.Time_diff.Rdata")
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
Time_genes.100 <- head(top_n(Time_diff, n = 100, desc(qval)), 100) %>% pull(gene_short_name) %>% as.character()
p.100 = plot_pseudotime_heatmap(cds[Time_genes.100,], num_clusters=4, show_rownames=T, return_heatmap=T)
plot_pseudotime_heatmap(cds[Time_genes.100,], num_clusters=3, show_rownames=T, return_heatmap=T)
ggsave("R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/FB_heatmap.pdf", p.100, width = 5, height = 10)
# extract genes in the heatmap
clustering <- data.frame(cutree(p.100$tree_row, k = 4))
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
write.table(clustering, "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/FB.pseudotime_heatmap.genesTop100.xlsx", sep = "\t")


## plot_pseudotime_heatmap: FB.HTS
# subset of HTS
cds_subset.HTS <- cds[, cds@phenoData@data$zt.scarType == "HTS"]
# 这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
Time_diff.HTS <- differentialGeneTest(cds_subset.HTS[ordergene,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes.HTS <- row.names(subset(Time_diff.HTS, qval < 1e-4))
p.HTS <- plot_pseudotime_heatmap(cds_subset.HTS[Time_genes.HTS,], num_clusters=4, show_rownames=F, return_heatmap=T)
# plot_genes_heatmap of TF, which is the intersection of 10 iterations of GeneRegulatoryNetwork
TF <- c('BHLHE40', 'CEBPB', 'CHD2', 'ATF3', 'CREB3L1', 'FOSB', 'ETS1', 'CREM', 'JUND', 'MAFG', 'FOS', 'ELF1', 'MXD4', 'NFIL3', 'IRF1', 'ETS2', 'TRIM28', 'NFKB1', 'JUNB', 'EZH2', 'XBP1', 'MAFF', 'FOSL2', 'NFE2L2', 'MEF2C', 'NR2F2', 'REL', 'PRDM1', 'YY1', 'TBX15')
p.HTS <- plot_pseudotime_heatmap(cds_subset.HTS[TF,], num_clusters=4, show_rownames=T, return_heatmap=T)


## plot_pseudotime_heatmap: FB.Keloid
# subset of Keloid
cds_subset.Keloid <- cds[, cds@phenoData@data$zt.scarType == "Keloid"]
# 这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
Time_diff.Keloid <- differentialGeneTest(cds_subset.Keloid[ordergene,], cores = 1, 
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes.Keloid <- row.names(subset(Time_diff.Keloid, qval < 1e-4))
p.Keloid <- plot_pseudotime_heatmap(cds_subset.Keloid[Time_genes.Keloid,], num_clusters=4, show_rownames=F, return_heatmap=T)
# plot_genes_heatmap of TF, which is the intersection of 10 iterations of GeneRegulatoryNetwork
TF <- c('ATF3', 'CREB3L1', 'CREB5', 'BHLHE40', 'MXD4', 'IRF8', 'CEBPB', 'REL', 'CEBPD', 'XBP1', 'CHD2', 'EGR1', 'EGR2', 'EGR3', 'FOS', 'FOSB', 'FOXD1', 'JUN', 'JUNB', 'JUND', 'MEF2C', 'NR2F2', 'STAT3', 'TBX15', 'TFAP2A', 'TRIM28')
p.Keloid <- plot_pseudotime_heatmap(cds_subset.Keloid[TF,], num_clusters=4, show_rownames=T, return_heatmap=T)

## check the expression of a specific gene in branches

genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c( "STAT3", "MYC")))

plot_genes_branched_pseudotime(cds[genes,],
                               branch_point = 1,
                               color_by = NULL,
                               ncol = 1)









## heatscatter of FB: density distribution of cells in plot_cell_trajectory
load("R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/cds.Rdata")
# FB
heatscatter(cds@reducedDimS[1,], cds@reducedDimS[2,],
            xlab = "Component 1", ylab= "Component 2", 
            colpal= "bl2gr2rd",
            add.contour = T, main = "FB")
# subset of HTS
cds_subset.HTS <- cds[, cds@phenoData@data$zt.scarType == "HTS"]
heatscatter(cds@reducedDimS[1, cds_subset.HTS@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.HTS@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", 
            colpal= "heat",
            add.contour = F, main = "HTS")
# subset of Keloid
cds_subset.Keloid <- cds[, cds@phenoData@data$zt.scarType == "Keloid"]
heatscatter(cds@reducedDimS[1, cds_subset.Keloid@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.Keloid@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", 
            colpal= "heat",
            add.contour = F, main = "Keloid")
# subset of Normal
cds_subset.Normal <- cds[, cds@phenoData@data$zt.scarType == "Normal"]
heatscatter(cds@reducedDimS[1, cds_subset.Normal@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.Normal@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", 
            colpal= "heat",
            add.contour = F, main = "Normal")
# subset of Normal_sc.hNC
cds$zt.scarType_dataSource <- paste(cds$zt.scarType, cds$zt.dataSource, sep = "_")
cds_subset.Normal_sc.hNC <- cds[, cds@phenoData@data$zt.scarType_dataSource == "Normal_sc.hNC"]
heatscatter(cds@reducedDimS[1, cds_subset.Normal_sc.hNC@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.Normal_sc.hNC@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", xlim = c(-10, 10), ylim = c(-5, 10),
            colpal= "heat",
            add.contour = F, main = "Normal_sc.hNC")
# subset of Normal_sc.kJID
cds_subset.Normal_sc.kJID <- cds[, cds@phenoData@data$zt.scarType_dataSource == "Normal_sc.kJID"]
heatscatter(cds@reducedDimS[1, cds_subset.Normal_sc.kJID@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.Normal_sc.kJID@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", 
            colpal= "heat",
            add.contour = F, main = "Normal_sc.kJID")
# subset of Normal_sc.kNC
cds_subset.Normal_sc.kNC <- cds[, cds@phenoData@data$zt.scarType_dataSource == "Normal_sc.kNC"]
heatscatter(cds@reducedDimS[1, cds_subset.Normal_sc.kNC@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.Normal_sc.kNC@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", 
            colpal= "heat",
            add.contour = F, main = "Normal_sc.kNC")
# add a new column to the metadata, in order to calculate cellNumbers of every cluster in each state, split by scarType
cds$state_zt.scarType_dataSource_FB6cluster <- paste(cds$State, cds$zt.scarType_dataSource, cds$integrated_snn_res.0.16,sep = "_")
cellNumbers <- as.data.frame(table(cds$state_zt.scarType_dataSource_FB6cluster))
colnames(cellNumbers) <- c("state_zt.scarType_dataSource_FB6cluster", "Freq")
write.table(cellNumbers, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/dFeature2000/cellNumbersOfFB6clustersInStates_splitByScarType.xlsx", sep = "\t")
# subset of FB6cluster.c0
cds_subset.c0 <- cds[, cds@phenoData@data$integrated_snn_res.0.16 == "0"]
heatscatter(cds@reducedDimS[1, cds_subset.c0@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.c0@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", xlim = c(-10, 10), ylim = c(-5, 10),
            colpal= "heat",
            add.contour = F, main = "c0")
# subset of FB6cluster.c1
cds_subset.c1 <- cds[, cds@phenoData@data$integrated_snn_res.0.16 == "1"]
heatscatter(cds@reducedDimS[1, cds_subset.c1@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.c1@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", xlim = c(-10, 10), ylim = c(-5, 10),
            colpal= "heat",
            add.contour = F, main = "c1")
# subset of FB6cluster.c2
cds_subset.c2 <- cds[, cds@phenoData@data$integrated_snn_res.0.16 == "2"]
heatscatter(cds@reducedDimS[1, cds_subset.c2@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.c2@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", xlim = c(-10, 10), ylim = c(-5, 10),
            colpal= "heat",
            add.contour = F, main = "c2")
# subset of FB6cluster.c3
cds_subset.c3 <- cds[, cds@phenoData@data$integrated_snn_res.0.16 == "3"]
heatscatter(cds@reducedDimS[1, cds_subset.c3@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.c3@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", xlim = c(-10, 10), ylim = c(-5, 10),
            colpal= "heat",
            add.contour = F, main = "c3")
# subset of FB6cluster.c4
cds_subset.c4 <- cds[, cds@phenoData@data$integrated_snn_res.0.16 == "4"]
heatscatter(cds@reducedDimS[1, cds_subset.c4@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.c4@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", xlim = c(-10, 10), ylim = c(-5, 10),
            colpal= "heat",
            add.contour = F, main = "c4")
# subset of FB6cluster.c5
cds_subset.c5 <- cds[, cds@phenoData@data$integrated_snn_res.0.16 == "5"]
heatscatter(cds@reducedDimS[1, cds_subset.c5@assayData$exprs@Dimnames[[2]]], cds@reducedDimS[2, cds_subset.c5@assayData$exprs@Dimnames[[2]]],
            xlab = "Component 1", ylab= "Component 2", xlim = c(-10, 10), ylim = c(-5, 10),
            colpal= "heat",
            add.contour = F, main = "c5")


#################################################


####################### FB_6clusters(FB.HVG).trajectoryAnalysis:monocle2 ##########################
# 
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# subtract raw data to establish a cds object
FB.counts <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16, assay = 'RNA', slot = 'counts')
FB.cell_metadata <- sc.kNC_kJID_hNC.FB_0.16@meta.data
FB.gene_annotation <- data.frame(gene_short_name = rownames(FB.counts), row.names = rownames(FB.counts))
pd <- new('AnnotatedDataFrame', data = FB.cell_metadata) 
fd <- new('AnnotatedDataFrame', data = FB.gene_annotation)
cds <- newCellDataSet(FB.counts,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

diff <- FindVariableFeatures(sc.kNC_kJID_hNC.FB_0.16, selection.method = "vst", nfeatures = 2000)
diff <- diff@assays$integrated@var.features
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)

## FB的HighVariableGenes作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
deg <- diff
head(deg)

## 轨迹构建基因可视化
load( file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/cds.Rdata")
ordergene <- deg
cds <- setOrderingFilter(cds, ordergene)  
#这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
plot_ordering_genes(cds)
#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

cds <- orderCells(cds)
#⚠️使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
#cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点
save(cds, file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/HVG2000/cds.Rdata")

load(file = "R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/HVG2000/cds.Rdata")
# check the trajectory by Pseudotime
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)
# check the trajectory by FB_6clusters
plot_cell_trajectory(cds,color_by="integrated_snn_res.0.16", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "integrated_snn_res.0.16") + facet_wrap("~integrated_snn_res.0.16", nrow = 1)
# check the trajectory by HTS/Keloid/Normal
plot_cell_trajectory(cds,color_by="zt.scarType", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "zt.scarType") + facet_wrap("~zt.scarType", nrow = 1)
# check the trajectory by dataSource
plot_cell_trajectory(cds,color_by="zt.dataSource", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "zt.dataSource") + facet_wrap("~zt.dataSource", nrow = 1)
# check the trajectory by Source
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)

## check the expression of a certain gene
load("R/sc.KNC_kJID_hNC/FB/trajectoryAnalysis/monocle2/unionofAllDEGofHTSandKeloid/cds.ordered.Rdata")
cds_subset <- cds[c("COL1A1", "ACTA2"),]
# the distribution of gene expression in different State/Pseudotime
plot_genes_in_pseudotime(cds_subset, color_by = "State")
plot_genes_in_pseudotime(cds_subset, color_by = "zt.scarType")
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
# the distribution of gene expression in the trajectory tree
pData(cds)$COL1A1 = log2(exprs(cds)['COL1A1',]+1)
plot_cell_trajectory(cds, color_by = "COL1A1", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$ACTA2 = log2(exprs(cds)['ACTA2',]+1)
plot_cell_trajectory(cds, color_by = "ACTA2", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$PCNA = log2(exprs(cds)['PCNA',]+1)
plot_cell_trajectory(cds, color_by = "PCNA", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red')
pData(cds)$MKI67 = log2(exprs(cds)['MKI67',]+1)
plot_cell_trajectory(cds, color_by = "MKI67", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red') + scale_alpha(range = c(0, 1))
pData(cds)$POSTN = log2(exprs(cds)['POSTN',]+1)
plot_cell_trajectory(cds, color_by = "POSTN", cell_size=0.5) + scale_color_gradient2(mid='grey', high = 'red') + scale_alpha(range = c(0, 1))


#################################################


####################### FB_6clusters.trajectoryAnalysis:monocle3 ##########################
# 
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# subtract raw data to establish a cds object
FB.counts <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16, assay = 'RNA', slot = 'counts')
FB.cell_metadata <- sc.kNC_kJID_hNC.FB_0.16@meta.data
FB.gene_annotation <- data.frame(gene_short_name = rownames(FB.counts), row.names = rownames(FB.counts))
cds <- new_cell_data_set(FB.counts,
                         cell_metadata = FB.cell_metadata,
                         gene_metadata = FB.gene_annotation)
# reduce dimensionality using PCA
cds <- preprocess_cds(cds, num_dim = 50)
# check the efficiency of the dimension number
plot_pc_variance_explained(cds)
# check batch effects
plot_cells(cds, color_cells_by="zt.dataSource", label_cell_groups=FALSE)
# remove batch effects by sampleName
cds = align_cds(cds, num_dim = 50, alignment_group = "orig.ident")
cds = reduce_dimension(cds)
plot_cells(cds, color_cells_by="zt.dataSource", label_cell_groups=FALSE)
# check FB.6clusters
plot_cells(cds, color_cells_by="integrated_snn_res.0.16", group_label_size = 5)
# check the expression of a particular gene
plot_cells(cds, color_cells_by="integrated_snn_res.0.16", group_label_size = 5, genes = c("POSTN")) + scale_color_gradient2(low = "blue", mid = "white", high = "red")
# clustering
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")
# learn trajectory
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "integrated_snn_res.0.16",
           label_groups_by_cluster=F, group_label_size = 3, 
           label_leaves=FALSE,
           label_branch_points=F)
# define the root/startpoint of the trajectory
cds <- order_cells(cds)
# show trajectory by pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=2)
# check the expression dynamics of some genes
test_genes <- c("POSTN", "COL1A1", "IL11RA")
test_lineage_cds <- cds[rowData(cds)$gene_short_name %in% test_genes,]
plot_genes_in_pseudotime(test_lineage_cds,
                         color_cells_by="pseudotime",
                         min_expr=0.5)


#################################################


#################### FB_6clusters.GeneRegulatoryNetwork #############################
# 提取表达矩阵，并命名为scRNA.csv，注意矩阵一定要转置，不然会报错 serving pySCENIC_linux workflow
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
write.csv(t(as.matrix(sc.kNC_kJID_hNC.FB_0.16@assays$RNA@counts)), file = "R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/pySCENIC_Linux_jobtest/scRNA.FB6clusters.csv")


library(SCENIC)
library(SCopeLoomR)
library(AUCell)


## Read information from loom file
loom <- open_loom('R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/HTS.pySCENIC_Linux/sample_SCENIC.loom')
# a list of regulons (TF and its targets)
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
write.table(names(regulons), "R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/HTS.pySCENIC_Linux/ROutput/regulons.TFname.xlsx", sep = "\t")
# regulon activity (& the corresponding threshold) in each cell
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
regulonAucThresholds <- data.frame(names(regulonAucThresholds), row.names = regulonAucThresholds)
embeddings <- get_embeddings(loom)
close_loom(loom)

## ComplexHeatmap of regulon activity
# Split complexHeatmap: dFeatureState
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
sc.kNC_kJID_hNC.FB_0.16.forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.expreMat <- t(scale(t(regulonAUC@assays@data$AUC)))
Idents(sc.kNC_kJID_hNC.FB_0.16.forHeatmap) <- "dFeatureState"
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16.forHeatmap$dFeatureState)
cluster_annotation = HeatmapAnnotation(dFeatureState = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info,
                                       col = list(dFeatureState = c("preBranch" = "#979797", "cellFate1" = "#df6166", "cellFate2" = "#7d90c4")))
Heatmap(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.expreMat,
        col = colorRamp2(c(-2, 0, 2), c("grey", "white", "red")),
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


## Binary heatmap of regulon activity
# Binarize regulonActivity
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds
  regulonsCells <- setNames(lapply(rownames(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x,]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),rownames(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  return(binaryRegulonActivity)
}
binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)
class(binaryRegulonActivity) <- "matrix"
save(binaryRegulonActivity, file = 'R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/pySCENIC_Linux/ROutput/binaryRegulonActivity.Rdata')
load('R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/pySCENIC_Linux/ROutput/binaryRegulonActivity.Rdata')

# Non-split complexHeatmap: zt.scarType
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
sc.kNC_kJID_hNC.FB_0.16.forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
Idents(sc.kNC_kJID_hNC.FB_0.16.forHeatmap) <- "zt.scarType"
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16.forHeatmap$zt.scarType)
cluster_annotation = HeatmapAnnotation(cluster = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info,
                                       col = list(cluster = c("HTS" = "#e97d72", "Keloid" = "#b49f33", "Normal" = "#53b64c")))
Heatmap(binaryRegulonActivity,
        col = c("white", "black"),
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

# Split complexHeatmap: zt.scarType
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
sc.kNC_kJID_hNC.FB_0.16.forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
Idents(sc.kNC_kJID_hNC.FB_0.16.forHeatmap) <- "zt.scarType"
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16.forHeatmap$zt.scarType)
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e97d72", "#b49f33", "#53b64c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(binaryRegulonActivity,
        col = c("white", "black"),
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

# Non-split complexHeatmap: dFeatureState
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
sc.kNC_kJID_hNC.FB_0.16.forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
Idents(sc.kNC_kJID_hNC.FB_0.16.forHeatmap) <- "dFeatureState"
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16.forHeatmap$dFeatureState)
cluster_annotation = HeatmapAnnotation(dFeatureState = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info,
                                       col = list(dFeatureState = c("preBranch" = "#979797", "cellFate1" = "#df6166", "cellFate2" = "#7d90c4")))
Heatmap(binaryRegulonActivity,
        col = c("white", "black"),
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

# Split complexHeatmap: dFeatureState
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
sc.kNC_kJID_hNC.FB_0.16.forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
Idents(sc.kNC_kJID_hNC.FB_0.16.forHeatmap) <- "dFeatureState"
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16.forHeatmap$dFeatureState)
cluster_annotation = HeatmapAnnotation(dFeatureState = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info,
                                       col = list(dFeatureState = c("preBranch" = "#979797", "cellFate1" = "#df6166", "cellFate2" = "#7d90c4")))
Heatmap(binaryRegulonActivity,
        col = c("white", "black"),
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


## Split complexHeatmap: zt.scarType & dFeatureState
# HTS
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'HTS')
binAct_subset <- binaryRegulonActivity[, sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS@assays$RNA@data@Dimnames[[2]]]
dim(binAct_subset)
Idents(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS) <- "dFeatureState"
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS$dFeatureState)
# define the annotation labels of heatmap
cluster_annotation <-  HeatmapAnnotation(dFeatureState = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS.cluster_info,
                                         col = list(dFeatureState = c("preBranch" = "#979797", "cellFate1" = "#df6166", "cellFate2" = "#7d90c4")))
Heatmap(binAct_subset, name = 'HTS',
        col = c("white", "black"),
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

# Keloid
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'Keloid')
binAct_subset <- binaryRegulonActivity[, sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Keloid@assays$RNA@data@Dimnames[[2]]]
dim(binAct_subset)
Idents(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Keloid) <- "dFeatureState"
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Keloid.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Keloid$dFeatureState)
# define the annotation labels of heatmap
cluster_annotation <-  HeatmapAnnotation(dFeatureState = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Keloid.cluster_info,
                                         col = list(dFeatureState = c("preBranch" = "#979797", "cellFate1" = "#df6166", "cellFate2" = "#7d90c4")))
Heatmap(binAct_subset, name = 'Keloid',
        col = c("white", "black"),
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Keloid.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))

# Normal
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Normal <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'Normal')
binAct_subset <- binaryRegulonActivity[, sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Normal@assays$RNA@data@Dimnames[[2]]]
dim(binAct_subset)
Idents(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Normal) <- "dFeatureState"
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Normal.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Normal$dFeatureState)
# define the annotation labels of heatmap
cluster_annotation <-  HeatmapAnnotation(dFeatureState = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Normal.cluster_info,
                                         col = list(dFeatureState = c("preBranch" = "#979797", "cellFate1" = "#df6166", "cellFate2" = "#7d90c4")))
Heatmap(binAct_subset, name = 'Normal',
        col = c("white", "black"),
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.Normal.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


#################################################


#################### FB_6clusters.HTS.GeneRegulatoryNetwork #############################
# 提取表达矩阵，并命名为scRNA.csv，注意矩阵一定要转置，不然会报错 serving pySCENIC_linux workflow
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
Idents(sc.kNC_kJID_hNC.FB_0.16) <- 'zt.scarType'
sc.kNC_kJID_hNC.FB_0.16.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'HTS')
write.csv(t(as.matrix(sc.kNC_kJID_hNC.FB_0.16.HTS@assays$RNA@counts)), file = "R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/HTS.pySCENIC_Linux/scRNA.FB6clusters.csv")


library(SCENIC)
library(SCopeLoomR)
library(AUCell)

# 提取sample_SCENIC.loom 信息
loom <- open_loom('R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/pySCENIC_Linux/sample_SCENIC.loom')
# 首先我们要把导入的loom处理成R中的数据
#1.获取regulon（TF与作用的genes）：提取每一个TF与每一个gene作用系数
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")   
regulons_incidMat[1:3,1:3]    #在这里就可以看出 每一个TF与每一个gene的作用数值
# a list of regulons (TF and its targets)
regulons <- regulonsToGeneLists(regulons_incidMat)
# regulons.DF <- rbind.fill(data.frame(t(data.frame(regulons[1]))), data.frame(t(data.frame(regulons[2]))))
write.table(names(regulons), "R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/pySCENIC_Linux/ROutput/regulons.TFname.xlsx", sep = "\t")
#2.获得regulon的AUC  即TF在每一个细胞的激活程度
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAUC[1:3,1:3]  #regulonAUC这个文件含有每一个TF在各个细胞中的表达量  列名为细胞名   行名为TF
#3.找出在这单细胞数据中 高表达的TF，将regulons-thresholds转换为dataframe
regulonAucThresholds <- get_regulon_thresholds(loom)
regulonAucThresholds <- data.frame(threshold = names(regulonAucThresholds), row.names = regulonAucThresholds)

tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])   #这里可以看出哪一些TF是在这个单细胞数据中高表达的
#4.这两步不知道是啥     
embeddings <- get_embeddings(loom)  #好像是什么嵌入   必须要做 否则后面会报错
close_loom(loom)

# 接下来我们就可以挑选自己感兴趣的TF，进行个性化分析
# 导入单细胞数据
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
scRNA <- sc.kNC_kJID_hNC.FB_0.16
Idents(scRNA) <- "zt.scarType"

# Expression pattern of regulons
scRNA@meta.data <- cbind(scRNA@meta.data, "SMAD1(+)" = getAUC(regulonAUC)["SMAD1(+)",])
DotPlot(scRNA, features = "SMAD1(+)")

# ComplexHeatmap
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- scRNA
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- t(scale(t(getAUC(regulonAUC))))
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix>2]=2
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix< -2]= -2
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_forHeatmap$zt.scarType)
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e97d72", "#b49f33", "#53b64c", "#54bcc2", "#6c9cf8"))
col <- matrix(c("#e97d72", "#b49f33", "#53b64c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


#################################################


#################### FB_6clusters.Keloid.GeneRegulatoryNetwork #############################
# 提取表达矩阵，并命名为scRNA.csv，注意矩阵一定要转置，不然会报错 serving pySCENIC_linux workflow
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
Idents(sc.kNC_kJID_hNC.FB_0.16) <- 'zt.scarType'
sc.kNC_kJID_hNC.FB_0.16.Keloid <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'Keloid')
write.csv(t(as.matrix(sc.kNC_kJID_hNC.FB_0.16.Keloid@assays$RNA@counts)), file = "R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/Keloid.pySCENIC_Linux/scRNA.FB6clusters.csv")


library(SCENIC)
library(SCopeLoomR)
library(AUCell)

# 提取sample_SCENIC.loom 信息
loom <- open_loom('R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/Keloid.pySCENIC_Linux/sample_SCENIC.loom')
# 首先我们要把导入的loom处理成R中的数据
#1.获取regulon（TF与作用的genes）：提取每一个TF与每一个gene作用系数
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")   
regulons_incidMat[1:3,1:3]    #在这里就可以看出 每一个TF与每一个gene的作用数值
# a list of regulons (TF and its targets)
regulons <- regulonsToGeneLists(regulons_incidMat)
# regulons.DF <- rbind.fill(data.frame(t(data.frame(regulons[1]))), data.frame(t(data.frame(regulons[2]))))
write.table(names(regulons), "R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/Keloid.pySCENIC_Linux/ROutput/regulons.TFname.xlsx", sep = "\t")
#2.获得regulon的AUC  即TF在每一个细胞的激活程度
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAUC[1:3,1:3]  #regulonAUC这个文件含有每一个TF在各个细胞中的表达量  列名为细胞名   行名为TF
#3.找出在这单细胞数据中 高表达的TF，将regulons-thresholds转换为dataframe
regulonAucThresholds <- get_regulon_thresholds(loom)
regulonAucThresholds <- data.frame(threshold = names(regulonAucThresholds), row.names = regulonAucThresholds)

tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])   #这里可以看出哪一些TF是在这个单细胞数据中高表达的
#4.这两步不知道是啥     
embeddings <- get_embeddings(loom)  #好像是什么嵌入   必须要做 否则后面会报错
close_loom(loom)

# 接下来我们就可以挑选自己感兴趣的TF，进行个性化分析
# 导入单细胞数据
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
scRNA <- sc.kNC_kJID_hNC.FB_0.16
Idents(scRNA) <- "zt.scarType"

# Expression pattern of regulons
scRNA@meta.data <- cbind(scRNA@meta.data, "SMAD1(+)" = getAUC(regulonAUC)["SMAD1(+)",])
DotPlot(scRNA, features = "SMAD1(+)")

# ComplexHeatmap
sc.kNC_kJID_hNC.FB_0.16_forHeatmap <- scRNA
# extract the scaled expression matrix of all genes in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix <- t(scale(t(getAUC(regulonAUC))))
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix>2]=2
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix[sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix< -2]= -2
# extract the cluster information of all cells in FB
sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16_forHeatmap$zt.scarType)
# define the annotation labels of heatmap
cluster_annotation <- HeatmapAnnotation(labels = levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info)))
# define the colors of cluster annotations
col <- matrix(c("#e97d72", "#b49f33", "#53b64c", "#54bcc2", "#6c9cf8"))
col <- matrix(c("#e97d72", "#b49f33", "#53b64c"))
rownames(col) <- levels(factor(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info))
# set the style of cluster annotation 
cluster_annotation <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col, lwd = 0)))
Heatmap(sc.kNC_kJID_hNC.FB_0.16_forHeatmap.exprMatrix, 
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16_forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))


#################################################


#################### FB_6clusters.randomSelected0.1_10times: GeneRegulatoryNetwork #############################
library(SCENIC)
library(SCopeLoomR)
library(AUCell)



## FB
# 提取表达矩阵，并命名为scRNA.csv，注意矩阵一定要转置，不然会报错 serving pySCENIC_linux workflow
for (i in 1:10) {
  load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
  selectedCells <- sample(colnames(sc.kNC_kJID_hNC.FB_0.16), floor(dim(sc.kNC_kJID_hNC.FB_0.16)[2] * 0.1))
  sc.kNC_kJID_hNC.FB_0.16 <- sc.kNC_kJID_hNC.FB_0.16[, selectedCells]
  filename <- paste0('R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/randomSelected0.1_10times/FB/', i, '.pySCENIC_Linux/scRNA.FB6clusters.csv')
  write.csv(t(as.matrix(sc.kNC_kJID_hNC.FB_0.16@assays$RNA@counts)), file = filename)
}
# get information from loom file
for (i in 1:10) {
  ## Read information from loom file
  loomName <- paste0('R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/randomSelected0.1_10times/FB/', i, '.pySCENIC_Linux/sample_SCENIC.loom')
  loom <- open_loom(loomName)
  # a list of regulons (TF and its targets)
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  fileName <- paste0('R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/randomSelected0.1_10times/FB/', i, '.pySCENIC_Linux/ROutput/FB', i, '-regulons.TFname.csv')
  write.csv(names(regulons), file = fileName)
  # # regulon activity (& the corresponding threshold) in each cell
  # regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
  # regulonAucThresholds <- get_regulon_thresholds(loom)
  # regulonAucThresholds <- data.frame(names(regulonAucThresholds), row.names = regulonAucThresholds)
  # embeddings <- get_embeddings(loom)
  close_loom(loom)
}


# HTS
for (i in 1:10) {
  load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
  Idents(sc.kNC_kJID_hNC.FB_0.16) <- 'zt.scarType'
  sc.kNC_kJID_hNC.FB_0.16 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'HTS')
  selectedCells <- sample(colnames(sc.kNC_kJID_hNC.FB_0.16), floor(dim(sc.kNC_kJID_hNC.FB_0.16)[2] * 0.1))
  sc.kNC_kJID_hNC.FB_0.16 <- sc.kNC_kJID_hNC.FB_0.16[, selectedCells]
  filename <- paste0('R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/randomSelected0.1_10times/HTS/', i, '.pySCENIC_Linux/scRNA.FB6clusters.csv')
  write.csv(t(as.matrix(sc.kNC_kJID_hNC.FB_0.16@assays$RNA@counts)), file = filename)
}


# Keloid
for (i in 1:10) {
  load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
  Idents(sc.kNC_kJID_hNC.FB_0.16) <- 'zt.scarType'
  sc.kNC_kJID_hNC.FB_0.16 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'Keloid')
  selectedCells <- sample(colnames(sc.kNC_kJID_hNC.FB_0.16), floor(dim(sc.kNC_kJID_hNC.FB_0.16)[2] * 0.1))
  sc.kNC_kJID_hNC.FB_0.16 <- sc.kNC_kJID_hNC.FB_0.16[, selectedCells]
  filename <- paste0('R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/randomSelected0.1_10times/Keloid/', i, '.pySCENIC_Linux/scRNA.FB6clusters.csv')
  write.csv(t(as.matrix(sc.kNC_kJID_hNC.FB_0.16@assays$RNA@counts)), file = filename)
}










## Read information from loom file
loom <- open_loom('R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/HTS.pySCENIC_Linux/sample_SCENIC.loom')
# a list of regulons (TF and its targets)
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
write.table(names(regulons), "R/sc.KNC_kJID_hNC/FB/GRN_TFanalysis/HTS.pySCENIC_Linux/ROutput/regulons.TFname.xlsx", sep = "\t")
# regulon activity (& the corresponding threshold) in each cell
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
regulonAucThresholds <- data.frame(names(regulonAucThresholds), row.names = regulonAucThresholds)
embeddings <- get_embeddings(loom)
close_loom(loom)




#################################################


#################### FB_6clusters.CellCellComunication: CellPhoneDB #############################
#### 10 cellTypes communications in HTS/Keloid/Normal
load('R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata')

### HTS
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.HTS <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "HTS")
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.HTS@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/HTS/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.HTS@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.HTS@meta.data[,'zt.cellType', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/HTS/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/HTS/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FB", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "HTS") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8))

## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/HTS/out/count_network.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$SOURCE)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
coords <- layout_in_circle(net, order = c('LL', 'EC', 'FB', 'spKC', 'ML', 'SMC', 'LEC', 'baKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 20 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/HTS/out/interaction_count.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)


### Keloid
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid")
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/Keloid/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid@meta.data[,'zt.cellType', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/Keloid/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/Keloid/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FB", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "Keloid") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8))

## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/Keloid/out/count_network.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$SOURCE)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
coords <- layout_in_circle(net, order = c('LL', 'EC', 'FB', 'spKC', 'ML', 'SMC', 'LEC', 'baKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 20 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/Keloid/out/interaction_count.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
V(net)$label.color <- "black"
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)


### selectedKeloid: randomly selected (the numbers of selectedKeloid.cells are the same as those of HTS)
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid")
selectedNumbers <- sample(table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['Keloid'], size = table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['HTS'])
sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid <- sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid[,selectedNumbers]
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data[,'zt.cellType', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FBs", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "Keloid") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8))


## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid/out/count_network.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$SOURCE)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
coords <- layout_in_circle(net, order = c('LL', 'EC', 'FB', 'spKC', 'ML', 'SMC', 'LEC', 'baKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 20 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid/out/interaction_count-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
V(net)$label.color <- "black"
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)


### Normal
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Normal <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Normal")
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.Normal@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/Normal/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.Normal@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.Normal@meta.data[,'zt.cellType', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/Normal/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/Normal/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FB", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "Normal") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8))

## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/Normal/out/count_network.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$cluster)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
# coords <- layout_in_circle(net, order = c('FB_0', 'FB_1', 'FB_2', 'FB_3', 'FB_4', 'FB_5', 'SMC', 'EC', 'LEC', 'baKC', 'LL', 'ML', 'spKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 20 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/Normal/out/interaction_count.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
V(net)$label.color <- "black"
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)








#### 9 cellTypes + FB6clusters communications in HTS/Keloid/Normal
# label cells
load('R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata')
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")

# add a new column/tag(cellType_FB6cluster) to the metadata
z <- merge(data.frame(sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType), data.frame(sc.kNC_kJID_hNC.FB_0.16$integrated_snn_res.0.16), by="row.names", all = T)
sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType_FB6cluster <- paste(z$sc.kNC_kJID_hNC.integrated_2000_0.1.zt.cellType, z$sc.kNC_kJID_hNC.FB_0.16.integrated_snn_res.0.16, sep = "_")
sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType_FB6cluster <- gsub("_NA", "", sc.kNC_kJID_hNC.integrated_2000_0.1$zt.cellType_FB6cluster)
save(sc.kNC_kJID_hNC.integrated_2000_0.1, file = "R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")



### HTS
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.HTS <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "HTS")
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.HTS@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/HTS/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.HTS@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.HTS@meta.data[,'zt.cellType_FB6cluster', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/HTS/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/HTS/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "HTS") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8.1))

## BubblePlot of interaction
## common interacting_pairs (selected) between HTS and Keloid
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/HTS/commonHTS.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/HTS/commonHTS.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'commonHTS') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## specific interacting_pairs (selected) in HTS
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/HTS/onlyHTS.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/HTS/onlyHTS.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyHTS') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## specific interacting_pairs (selected, eliminateProteinComplex) in HTS
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/HTS/onlyHTS.FBc012-means_significant_eliminate protein complex.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/HTS/onlyHTS.FBc012-pvalues_significant_eliminate protein complex.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.EC', 'FB_0.LL', 'FB_0.ML', 'FB_1.EC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.LL', 'FB_1.ML', 'FB_2.EC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.LL', 'FB_2.ML')) # change the order of axis labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyHTS') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.4, limit = c(0, 0.8)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## common interacting_pairs (selected) between HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/HTS/commonHTS.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/HTS/commonHTS.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('EGFR_HBEGF', 'EGFR_MIF', 'EGFR_TGFB1', 'EGFR_COPA', 'EGFR_AREG', 'TNFRSF1A_GRN', 'MIF_TNFRSF14', 'TNFSF12_TNFRSF12A', 'NRP2_VEGFA', 'VEGFA_FLT1', 'JAG1_NOTCH3', 'CD46_JAG1', 'ESAM_ESAM', 'CD44_SELE', 'CD55_ADGRE5', 'CCL2_ACKR1', 'COL4A2_a2b1 complex', 'COL18A1_a2b1 complex', 'CXCL8_ACKR1', 'ACKR3_CXCL12', 'PROS1_AXL')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'commonHTS') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## interacting_pairs (selected) only in HTS (comparing HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations)
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/HTS/onlyHTS.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/HTS/onlyHTS.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('FGFR1_FGFR2', 'FGF7_FGFR2', 'CD44_FGFR2', 'TIMP1_FGFR2', 'TGFB1_TGFbeta receptor1', 'TGFB3_TGFbeta receptor1', 'PDGFA_PDGFRA', 'PDGFRB_PDGFD', 'PDGFR complex_PDGFD', 'NRP1_VEGFA', 'NRP1_VEGFB', 'CD74_APP', 'CD74_MIF', 'CD74_COPA', 'CD44_HBEGF', 'TNFRSF10B_TNFSF10', 'IGF1_IGF1R', 'MDK_LRP1', 'LAMP1_FAM3C', 'NRP1_SEMA3A', 'AXL_GAS6', 'COL4A1_a2b1 complex', 'COL4A1_a11b1 complex', 'COL4A2_a11b1 complex', 'COL5A1_a11b1 complex', 'COL6A2_a11b1 complex', 'COL6A3_a11b1 complex', 'COL12A1_a11b1 complex', 'COL14A1_a11b1 complex', 'COL15A1_a11b1 complex', 'COL16A1_a2b1 complex', 'COL16A1_a11b1 complex', 'COL18A1_a11b1 complex', 'FN1_a11b1 complex', 'JAM2_JAM3', 'JAM3_JAM3', 'LAMC1_a2b1 complex', 'LAMC1_a6b1 complex', 'ICAM1_AREG', 'CADM3_CADM1', 'CXCL2_DPP4', 'CXCL12_CXCR4', 'NRP2_SEMA3C', 'HLA-DRB1_OGN', 'HLA-C_FAM3C', 'HLA-DPB1_NRG1', 'IL6_HRH1', 'SELE_GLG1')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyHTS') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))


## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/HTS/out/count_network-HTS-orderChanged.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$cluster)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
# coords <- layout_in_circle(net, order = c('FB_0', 'FB_1', 'FB_2', 'FB_3', 'FB_4', 'FB_5', 'SMC', 'EC', 'LEC', 'baKC', 'LL', 'ML', 'spKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 10 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/HTS/out/interaction_count-HTS-orderChanged.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
V(net)$label.color <- "black"
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)


## Keloid
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid")
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/Keloid/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid@meta.data[,'zt.cellType_FB6cluster', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/Keloid/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/Keloid/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "Keloid") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8.1))

## BubblePlot of interaction
## common interacting_pairs (selected) between HTS and Keloid
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/commonKeloid.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/commonKeloid.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'commonKeloid') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## specific interacting_pairs (selected) in Keloid
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/onlyKeloid.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/onlyKeloid.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyKeloid') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## specific interacting_pairs (selected, eliminateProteinComplex) in Keloid
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/onlyKeloid.FBc012-means_significant_eliminate protein complex.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/onlyKeloid.FBc012-pvalues_significant_eliminate protein complex.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.EC', 'FB_0.LL', 'FB_0.ML', 'FB_1.EC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.LL', 'FB_1.ML', 'FB_2.EC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.LL', 'FB_2.ML')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('FGFR3_EPHA4', 'FGF18_FGFR3', 'FGFR3_NECTIN1', 'FGF7_FGFR3', 'FGF2_FGFR3', 'FGFR1_FGF7', 'FGF18_FGFR1', 'FGF2_CD44', 'FGF2_FGFR1', 'DLL1_NOTCH1', 'DLL1_NOTCH2', 'DLL1_NOTCH3', 'DLL1_NOTCH4', 'DLL4_NOTCH3', 'NOTCH1_DLL4', 'NOTCH2_DLL4', 'NOTCH1_JAG1', 'NOTCH1_JAG2', 'NOTCH2_JAG2', 'NOTCH3_JAG2', 'NOTCH4_JAG2', 'JAG1_NOTCH4', 'NOTCH1_NOV', 'EFNA1_EPHA4', 'EFNA1_EPHA2', 'EPHA2_EFNA3', 'EPHA4_EFNB1', 'EFNB2_EPHB4', 'EFNB2_EPHB3', 'EFNB2_EPHA4', 'EPHB3_EFNB1', 'EPHB4_EFNB1', 'TGFB1_TGFbeta receptor2', 'TGFB1_TGFBR3', 'TGFB3_TGFBR3', 'ACVR_1B2A receptor_GDF11', 'ACVR_1A2A receptor_INHBB', 'FLT1 complex_PGF', 'NRP1_PGF', 'NRP2_PGF', 'FLT1_PGF', 'PDGFB_PDGFRB', 'ADRB2_VEGFB', 'VEGFA_KDR', 'FLT1 complex_VEGFA', 'LGALS9_PTPRK', 'LGALS9_CD44', 'LGALS9_CD47', 'SIRPA_CD47', 'TNFRSF1B_GRN', 'CXADR_FAM3C', 'SELP_CD34', 'SEMA4A_PLXND1', 'PLXNB2_SEMA4C', 'ANGPT2_TEK', 'LRP6_CKLF', 'CSF1_SLC7A1', 'MIF_TNFRSF10D', 'NECTIN1_NECTIN3', 'CDH1_a2b1 complex', 'DSG1_DSC2', 'DSG1_DSC3', 'IL6 receptor_IL6', 'IGF1_a6b4 complex', 'LAMC1_a7b1 complex', 'FN1_a2b1 complex')) # change the order of axis labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyKeloid') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.4, limit = c(0, 0.8)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## common interacting_pairs (selected) between HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/Keloid/commonKeloid.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/Keloid/commonKeloid.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('EGFR_HBEGF', 'EGFR_MIF', 'EGFR_TGFB1', 'EGFR_COPA', 'EGFR_AREG', 'TNFRSF1A_GRN', 'MIF_TNFRSF14', 'TNFSF12_TNFRSF12A', 'NRP2_VEGFA', 'VEGFA_FLT1', 'JAG1_NOTCH3', 'CD46_JAG1', 'ESAM_ESAM', 'CD44_SELE', 'CD55_ADGRE5', 'CCL2_ACKR1', 'COL4A2_a2b1 complex', 'COL18A1_a2b1 complex', 'CXCL8_ACKR1', 'ACKR3_CXCL12', 'PROS1_AXL')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'commonKeloid') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## interacting_pairs (selected) only in Keloid (comparing HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations)
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/Keloid/onlyKeloid.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/Keloid/onlyKeloid.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('FGFR3_EPHA4', 'FGFR3_NECTIN1', 'FGF7_FGFR3', 'FGF18_FGFR1', 'FGF18_FGFR3', 'EFNA1_EPHA2', 'EFNB2_EPHB3', 'EFNB2_EPHB4', 'EFNB2_EPHA4', 'EPHA4_EFNB1', 'EPHB3_EFNB1', 'EPHB4_EFNB1', 'NOTCH1_JAG1', 'NOTCH2_JAG2', 'NOTCH1_DLL4', 'NOTCH2_DLL4', 'DLL1_NOTCH1', 'DLL1_NOTCH2', 'DLL1_NOTCH3', 'DLL1_NOTCH4', 'DLL4_NOTCH3', 'LGALS9_CD47', 'LGALS9_LRP1', 'LGALS9_SLC1A5', 'LGALS9_PTPRK', 'LGALS9_CD44', 'LGALS9_MRC2', 'TNFRSF1B_GRN', 'TNFRSF10A_TNFSF10', 'TNFSF12_TNFRSF25', 'MIF_TNFRSF10D', 'ADRB2_VEGFB', 'VEGFA_KDR', 'FLT1 complex_VEGFA', 'FLT1 complex_VEGFB', 'NRP2_PGF', 'FLT1_PGF', 'FLT1 complex_PGF', 'PDGFB_PDGFRB', 'PDGFB_LRP1', 'ANGPT2_TEK', 'TGFB3_TGFBR3', 'COL5A3_a1b1 complex', 'COL5A3_a2b1 complex', 'COL7A1_a1b1 complex', 'COL7A1_a2b1 complex', 'COL7A1_a10b1 complex', 'COL17A1_a1b1 complex', 'COL17A1_a2b1 complex', 'COL17A1_a10b1 complex', 'COL21A1_a1b1 complex', 'COL21A1_a2b1 complex', 'COL27A1_a1b1 complex', 'COL27A1_a2b1 complex', 'DSG1_DSC2', 'DSG1_DSC3', 'CXADR_FAM3C', 'NECTIN1_NECTIN3', 'CDH1_a2b1 complex', 'SEMA4A_PLXND1', 'PLXNB2_SEMA4C', 'SELP_CD34', 'IL6 receptor_IL6', 'CSF1_SLC7A1', 'LRP6_CKLF', 'SIRPA_CD47')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyKeloid') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))


## selectedKeloid: randomly selected (the numbers of selectedKeloid.cells are the same as those of HTS)
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid")
selectedNumbers <- sample(table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['Keloid'], size = table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['HTS'])
sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid <- sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid[,selectedNumbers]
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data[,'zt.cellType_FB6cluster', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "selectedKeloid") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8.1))

## BubblePlot of interaction
## common interacting_pairs (selected) between HTS and Keloid
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/commonKeloid.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/commonKeloid.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'commonKeloid') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## specific interacting_pairs (selected) in Keloid
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/onlyKeloid.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/onlyKeloid.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyKeloid') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## specific interacting_pairs (selected, eliminateProteinComplex) in Keloid
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/onlyKeloid.FBc012-means_significant_eliminate protein complex.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/Keloid/onlyKeloid.FBc012-pvalues_significant_eliminate protein complex.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.EC', 'FB_0.LL', 'FB_0.ML', 'FB_1.EC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.LL', 'FB_1.ML', 'FB_2.EC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.LL', 'FB_2.ML')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('FGFR3_EPHA4', 'FGF18_FGFR3', 'FGFR3_NECTIN1', 'FGF7_FGFR3', 'FGF2_FGFR3', 'FGFR1_FGF7', 'FGF18_FGFR1', 'FGF2_CD44', 'FGF2_FGFR1', 'DLL1_NOTCH1', 'DLL1_NOTCH2', 'DLL1_NOTCH3', 'DLL1_NOTCH4', 'DLL4_NOTCH3', 'NOTCH1_DLL4', 'NOTCH2_DLL4', 'NOTCH1_JAG1', 'NOTCH1_JAG2', 'NOTCH2_JAG2', 'NOTCH3_JAG2', 'NOTCH4_JAG2', 'JAG1_NOTCH4', 'NOTCH1_NOV', 'EFNA1_EPHA4', 'EFNA1_EPHA2', 'EPHA2_EFNA3', 'EPHA4_EFNB1', 'EFNB2_EPHB4', 'EFNB2_EPHB3', 'EFNB2_EPHA4', 'EPHB3_EFNB1', 'EPHB4_EFNB1', 'TGFB1_TGFbeta receptor2', 'TGFB1_TGFBR3', 'TGFB3_TGFBR3', 'ACVR_1B2A receptor_GDF11', 'ACVR_1A2A receptor_INHBB', 'FLT1 complex_PGF', 'NRP1_PGF', 'NRP2_PGF', 'FLT1_PGF', 'PDGFB_PDGFRB', 'ADRB2_VEGFB', 'VEGFA_KDR', 'FLT1 complex_VEGFA', 'LGALS9_PTPRK', 'LGALS9_CD44', 'LGALS9_CD47', 'SIRPA_CD47', 'TNFRSF1B_GRN', 'CXADR_FAM3C', 'SELP_CD34', 'SEMA4A_PLXND1', 'PLXNB2_SEMA4C', 'ANGPT2_TEK', 'LRP6_CKLF', 'CSF1_SLC7A1', 'MIF_TNFRSF10D', 'NECTIN1_NECTIN3', 'CDH1_a2b1 complex', 'DSG1_DSC2', 'DSG1_DSC3', 'IL6 receptor_IL6', 'IGF1_a6b4 complex', 'LAMC1_a7b1 complex', 'FN1_a2b1 complex')) # change the order of axis labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyKeloid') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.4, limit = c(0, 0.8)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## common interacting_pairs (selected) between HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid/commonselectedKeloid.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid/commonselectedKeloid.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('EGFR_HBEGF', 'EGFR_MIF', 'EGFR_TGFB1', 'EGFR_COPA', 'EGFR_AREG', 'TNFRSF1A_GRN', 'MIF_TNFRSF14', 'TNFSF12_TNFRSF12A', 'NRP2_VEGFA', 'VEGFA_FLT1', 'JAG1_NOTCH3', 'CD46_JAG1', 'ESAM_ESAM', 'CD44_SELE', 'CD55_ADGRE5', 'CCL2_ACKR1', 'COL4A2_a2b1 complex', 'COL18A1_a2b1 complex', 'CXCL8_ACKR1', 'ACKR3_CXCL12', 'PROS1_AXL')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'commonselectedKeloid') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## interacting_pairs (selected) only in Keloid (comparing HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations)
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid/onlyselectedKeloid.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid/onlyselectedKeloid.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('FGFR3_EPHA4', 'FGFR3_NECTIN1', 'FGF7_FGFR3', 'FGF18_FGFR1', 'FGF18_FGFR3', 'EFNA1_EPHA2', 'EFNB2_EPHB3', 'EFNB2_EPHB4', 'EFNB2_EPHA4', 'EPHA4_EFNB1', 'EPHB3_EFNB1', 'EPHB4_EFNB1', 'NOTCH1_JAG1', 'NOTCH2_JAG2', 'NOTCH1_DLL4', 'NOTCH2_DLL4', 'DLL1_NOTCH1', 'DLL1_NOTCH2', 'DLL1_NOTCH3', 'DLL1_NOTCH4', 'DLL4_NOTCH3', 'LGALS9_CD47', 'LGALS9_LRP1', 'LGALS9_SLC1A5', 'LGALS9_PTPRK', 'LGALS9_CD44', 'LGALS9_MRC2', 'TNFRSF1B_GRN', 'TNFRSF10A_TNFSF10', 'TNFSF12_TNFRSF25', 'MIF_TNFRSF10D', 'ADRB2_VEGFB', 'VEGFA_KDR', 'FLT1 complex_VEGFA', 'FLT1 complex_VEGFB', 'NRP2_PGF', 'FLT1_PGF', 'FLT1 complex_PGF', 'PDGFB_PDGFRB', 'PDGFB_LRP1', 'ANGPT2_TEK', 'TGFB3_TGFBR3', 'COL5A3_a1b1 complex', 'COL5A3_a2b1 complex', 'COL7A1_a1b1 complex', 'COL7A1_a2b1 complex', 'COL7A1_a10b1 complex', 'COL17A1_a1b1 complex', 'COL17A1_a2b1 complex', 'COL17A1_a10b1 complex', 'COL21A1_a1b1 complex', 'COL21A1_a2b1 complex', 'COL27A1_a1b1 complex', 'COL27A1_a2b1 complex', 'DSG1_DSC2', 'DSG1_DSC3', 'CXADR_FAM3C', 'NECTIN1_NECTIN3', 'CDH1_a2b1 complex', 'SEMA4A_PLXND1', 'PLXNB2_SEMA4C', 'SELP_CD34', 'IL6 receptor_IL6', 'CSF1_SLC7A1', 'LRP6_CKLF', 'SIRPA_CD47')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyselectedKeloid') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))

## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid/out/count_network-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$cluster)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
# coords <- layout_in_circle(net, order = c('FB_0', 'FB_1', 'FB_2', 'FB_3', 'FB_4', 'FB_5', 'SMC', 'EC', 'LEC', 'baKC', 'LL', 'ML', 'spKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 10 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid/out/interaction_count-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
V(net)$label.color <- "black"
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)



## Normal
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Normal <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Normal")
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.Normal@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/Normal/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.Normal@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.Normal@meta.data[,'zt.cellType_FB6cluster', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/Normal/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/Normal/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "Normal") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8.1))


#################################################


#################### FB_6clusters.CellCellComunication: CellPhoneDB--selectedKeloid_3iterations #############################
#### 10 cellTypes communications in selectedKeloid, whose numbers are the same as the numbers of HTS

### Iteration 1
load('R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata')
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid")
selectedNumbers <- sample(table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['Keloid'], size = table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['HTS'])
sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid <- sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid[,selectedNumbers]
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration1/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data[,'zt.cellType', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration1/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration1/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FBs", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "Keloid") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8))

## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration1/out/count_network-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$SOURCE)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
coords <- layout_in_circle(net, order = c('LL', 'EC', 'FB', 'spKC', 'ML', 'SMC', 'LEC', 'baKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 20 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration1/out/interaction_count-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
V(net)$label.color <- "black"
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)


### Iteration 2
load('R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata')
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid")
selectedNumbers <- sample(table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['Keloid'], size = table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['HTS'])
sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid <- sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid[,selectedNumbers]
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration2/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data[,'zt.cellType', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration2/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration2/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FBs", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "Keloid") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8))

## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration2/out/count_network-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$SOURCE)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
coords <- layout_in_circle(net, order = c('LL', 'EC', 'FB', 'spKC', 'ML', 'SMC', 'LEC', 'baKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 20 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration2/out/interaction_count-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
V(net)$label.color <- "black"
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)


### Iteration 3
load('R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata')
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid")
selectedNumbers <- sample(table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['Keloid'], size = table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['HTS'])
sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid <- sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid[,selectedNumbers]
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration3/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data[,'zt.cellType', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration3/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration3/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FBs", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "Keloid") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8))

## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration3/out/count_network-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$SOURCE)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
coords <- layout_in_circle(net, order = c('LL', 'EC', 'FB', 'spKC', 'ML', 'SMC', 'LEC', 'baKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 20 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/selectedKeloid_3iterations/Iteration3/out/interaction_count-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
V(net)$label.color <- "black"
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)




#### 9 cellTypes + FB6clusters communications in HTS/Keloid/Normal

### Iteration 1
load('R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata')
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid")
selectedNumbers <- sample(table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['Keloid'], size = table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['HTS'])
sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid <- sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid[,selectedNumbers]
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration1/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data[,'zt.cellType_FB6cluster', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration1/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration1/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "selectedKeloid") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8.1))

## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration1/out/count_network-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$cluster)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
# coords <- layout_in_circle(net, order = c('FB_0', 'FB_1', 'FB_2', 'FB_3', 'FB_4', 'FB_5', 'SMC', 'EC', 'LEC', 'baKC', 'LL', 'ML', 'spKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 10 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration1/out/interaction_count-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
V(net)$label.color <- "black"
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)

## BubblePlot of interaction
## common interacting_pairs (selected) between HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration1/commonselectedKeloid_iteration1.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration1/commonselectedKeloid_iteration1.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('EGFR_HBEGF', 'EGFR_MIF', 'EGFR_TGFB1', 'EGFR_COPA', 'EGFR_AREG', 'TNFRSF1A_GRN', 'MIF_TNFRSF14', 'TNFSF12_TNFRSF12A', 'NRP2_VEGFA', 'VEGFA_FLT1', 'JAG1_NOTCH3', 'CD46_JAG1', 'ESAM_ESAM', 'CD44_SELE', 'CD55_ADGRE5', 'CCL2_ACKR1', 'COL4A2_a2b1 complex', 'COL18A1_a2b1 complex', 'CXCL8_ACKR1', 'ACKR3_CXCL12', 'PROS1_AXL')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'commonselectedKeloid_iteration1') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## interacting_pairs (selected) only in Keloid (comparing HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations)
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration1/onlyselectedKeloid_iteration1.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration1/onlyselectedKeloid_iteration1.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('FGFR3_EPHA4', 'FGFR3_NECTIN1', 'FGF7_FGFR3', 'FGF18_FGFR1', 'FGF18_FGFR3', 'EFNA1_EPHA2', 'EFNB2_EPHB3', 'EFNB2_EPHB4', 'EFNB2_EPHA4', 'EPHA4_EFNB1', 'EPHB3_EFNB1', 'EPHB4_EFNB1', 'NOTCH1_JAG1', 'NOTCH2_JAG2', 'NOTCH1_DLL4', 'NOTCH2_DLL4', 'DLL1_NOTCH1', 'DLL1_NOTCH2', 'DLL1_NOTCH3', 'DLL1_NOTCH4', 'DLL4_NOTCH3', 'LGALS9_CD47', 'LGALS9_LRP1', 'LGALS9_SLC1A5', 'LGALS9_PTPRK', 'LGALS9_CD44', 'LGALS9_MRC2', 'TNFRSF1B_GRN', 'TNFRSF10A_TNFSF10', 'TNFSF12_TNFRSF25', 'MIF_TNFRSF10D', 'ADRB2_VEGFB', 'VEGFA_KDR', 'FLT1 complex_VEGFA', 'FLT1 complex_VEGFB', 'NRP2_PGF', 'FLT1_PGF', 'FLT1 complex_PGF', 'PDGFB_PDGFRB', 'PDGFB_LRP1', 'ANGPT2_TEK', 'TGFB3_TGFBR3', 'COL5A3_a1b1 complex', 'COL5A3_a2b1 complex', 'COL7A1_a1b1 complex', 'COL7A1_a2b1 complex', 'COL7A1_a10b1 complex', 'COL17A1_a1b1 complex', 'COL17A1_a2b1 complex', 'COL17A1_a10b1 complex', 'COL21A1_a1b1 complex', 'COL21A1_a2b1 complex', 'COL27A1_a1b1 complex', 'COL27A1_a2b1 complex', 'DSG1_DSC2', 'DSG1_DSC3', 'CXADR_FAM3C', 'NECTIN1_NECTIN3', 'CDH1_a2b1 complex', 'SEMA4A_PLXND1', 'PLXNB2_SEMA4C', 'SELP_CD34', 'IL6 receptor_IL6', 'CSF1_SLC7A1', 'LRP6_CKLF', 'SIRPA_CD47')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyselectedKeloid_iteration1') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))


### Iteration 2
load('R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata')
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid")
selectedNumbers <- sample(table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['Keloid'], size = table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['HTS'])
sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid <- sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid[,selectedNumbers]
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration2/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data[,'zt.cellType_FB6cluster', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration2/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration2/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "selectedKeloid") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8.1))

## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration2/out/count_network-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$cluster)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
# coords <- layout_in_circle(net, order = c('FB_0', 'FB_1', 'FB_2', 'FB_3', 'FB_4', 'FB_5', 'SMC', 'EC', 'LEC', 'baKC', 'LL', 'ML', 'spKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 10 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration2/out/interaction_count-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
V(net)$label.color <- "black"
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)

## BubblePlot of interaction
## common interacting_pairs (selected) between HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration2/commonselectedKeloid_iteration2.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration2/commonselectedKeloid_iteration2.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('EGFR_HBEGF', 'EGFR_MIF', 'EGFR_TGFB1', 'EGFR_COPA', 'EGFR_AREG', 'TNFRSF1A_GRN', 'MIF_TNFRSF14', 'TNFSF12_TNFRSF12A', 'NRP2_VEGFA', 'VEGFA_FLT1', 'JAG1_NOTCH3', 'CD46_JAG1', 'ESAM_ESAM', 'CD44_SELE', 'CD55_ADGRE5', 'CCL2_ACKR1', 'COL4A2_a2b1 complex', 'COL18A1_a2b1 complex', 'CXCL8_ACKR1', 'ACKR3_CXCL12', 'PROS1_AXL')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'commonselectedKeloid_iteration2') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## interacting_pairs (selected) only in Keloid (comparing HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations)
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration2/onlyselectedKeloid_iteration2.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration2/onlyselectedKeloid_iteration2.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('FGFR3_EPHA4', 'FGFR3_NECTIN1', 'FGF7_FGFR3', 'FGF18_FGFR1', 'FGF18_FGFR3', 'EFNA1_EPHA2', 'EFNB2_EPHB3', 'EFNB2_EPHB4', 'EFNB2_EPHA4', 'EPHA4_EFNB1', 'EPHB3_EFNB1', 'EPHB4_EFNB1', 'NOTCH1_JAG1', 'NOTCH2_JAG2', 'NOTCH1_DLL4', 'NOTCH2_DLL4', 'DLL1_NOTCH1', 'DLL1_NOTCH2', 'DLL1_NOTCH3', 'DLL1_NOTCH4', 'DLL4_NOTCH3', 'LGALS9_CD47', 'LGALS9_LRP1', 'LGALS9_SLC1A5', 'LGALS9_PTPRK', 'LGALS9_CD44', 'LGALS9_MRC2', 'TNFRSF1B_GRN', 'TNFRSF10A_TNFSF10', 'TNFSF12_TNFRSF25', 'MIF_TNFRSF10D', 'ADRB2_VEGFB', 'VEGFA_KDR', 'FLT1 complex_VEGFA', 'FLT1 complex_VEGFB', 'NRP2_PGF', 'FLT1_PGF', 'FLT1 complex_PGF', 'PDGFB_PDGFRB', 'PDGFB_LRP1', 'ANGPT2_TEK', 'TGFB3_TGFBR3', 'COL5A3_a1b1 complex', 'COL5A3_a2b1 complex', 'COL7A1_a1b1 complex', 'COL7A1_a2b1 complex', 'COL7A1_a10b1 complex', 'COL17A1_a1b1 complex', 'COL17A1_a2b1 complex', 'COL17A1_a10b1 complex', 'COL21A1_a1b1 complex', 'COL21A1_a2b1 complex', 'COL27A1_a1b1 complex', 'COL27A1_a2b1 complex', 'DSG1_DSC2', 'DSG1_DSC3', 'CXADR_FAM3C', 'NECTIN1_NECTIN3', 'CDH1_a2b1 complex', 'SEMA4A_PLXND1', 'PLXNB2_SEMA4C', 'SELP_CD34', 'IL6 receptor_IL6', 'CSF1_SLC7A1', 'LRP6_CKLF', 'SIRPA_CD47')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyselectedKeloid_iteration2') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))


### Iteration 3
load('R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata')
## Heatmap
# create data for linux
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid")
selectedNumbers <- sample(table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['Keloid'], size = table(sc.kNC_kJID_hNC.integrated_2000_0.1@meta.data$zt.scarType)['HTS'])
sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid <- sc.kNC_kJID_hNC.integrated_2000_0.1.Keloid[,selectedNumbers]
write.table(as.matrix(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@assays$RNA@data), 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration3/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data), sc.kNC_kJID_hNC.integrated_2000_0.1.selectedKeloid@meta.data[,'zt.cellType_FB6cluster', drop=F])  
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration3/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# symetrical heatmap of interaction
heatmap.data <- read.table('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration3/out/count_network.txt', header = T)
heatmap.data$SOURCE <- factor(heatmap.data$SOURCE, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU")) # change the order of axis labels
heatmap.data$TARGET <- factor(heatmap.data$TARGET, levels = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "EC", "baKC", "spKC", "SMC", "ML", "LL", "LEC", "SGC", "MELA/NEU"))
ggplot(data = heatmap.data, aes(x=SOURCE, y=TARGET, fill=log2(count))) + labs(caption = "selectedKeloid") +
  geom_tile() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6",
                       midpoint = 4, 
                       limit = c(0, 8.1))

## Network
#输入文件，并确定数据中所有行的count值均大于0
network_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration3/out/count_network-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE)
network_df <- network_df %>% filter(count > 0)
#建立网络
net <- graph_from_data_frame(network_df)
#提取顶点位置节点的坐标
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
#  选取ggsci保重的d3-category20c调色板
col <- pal_d3("category20c")(length(unique(network_df$SOURCE)))
names(col) <- unique(network_df$cluster)
#计算社区结构，并按照社区结构对位点进行排列
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = membership(group))
# coords <- layout_in_circle(net, order = c('FB_0', 'FB_1', 'FB_2', 'FB_3', 'FB_4', 'FB_5', 'SMC', 'EC', 'LEC', 'baKC', 'LL', 'ML', 'spKC', 'SGC', 'MELA/NEU'))
#计算edge的权重
E(net)$width <- E(net)$count / 10 #根据count值设置边的宽
# define the node size, in proportion to its total interaction counts
nodeInteractionCount_df <- read.table("R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/selectedKeloid_3iterations/Iteration3/out/interaction_count-selectedKeloid-orderChanged.txt", sep = "\t", quote = "", header = TRUE, col.names = c('cellType', 'totalCount'))
vextex.size <- nodeInteractionCount_df[, 'totalCount'] / 20
#调节线条和label的配色
E(net)$color <- col #用前面设置好的颜色赋给连线
V(net)$label.color <- "black"
#调整节点位置的线条角度，这里一开始也挺挠头的，因为默认定点的线是一致向右的，感谢cellchat，在这个软件的源代码中找到了如何调整，很巧妙的一个方法
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
plot(net,
     edge.arrow.size = 0, #连线箭头
     edge.curved = 0.2, #连线弯曲
     vertex.color = col,
     vertex.frame.color = F, #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = vextex.size)

## BubblePlot of interaction
## common interacting_pairs (selected) between HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration3/commonselectedKeloid_iteration3.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration3/commonselectedKeloid_iteration3.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('EGFR_HBEGF', 'EGFR_MIF', 'EGFR_TGFB1', 'EGFR_COPA', 'EGFR_AREG', 'TNFRSF1A_GRN', 'MIF_TNFRSF14', 'TNFSF12_TNFRSF12A', 'NRP2_VEGFA', 'VEGFA_FLT1', 'JAG1_NOTCH3', 'CD46_JAG1', 'ESAM_ESAM', 'CD44_SELE', 'CD55_ADGRE5', 'CCL2_ACKR1', 'COL4A2_a2b1 complex', 'COL18A1_a2b1 complex', 'CXCL8_ACKR1', 'ACKR3_CXCL12', 'PROS1_AXL')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'commonselectedKeloid_iteration3') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))
## interacting_pairs (selected) only in Keloid (comparing HTS, Keloid, selectedKeloid, and selectedKeloid_3iterations)
# data cleansing
means.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration3/onlyselectedKeloid_iteration3.FBc012-means.xlsx', sheetIndex = 1)
row.names(means.data) <- means.data[,1]
means.data <- means.data[,-1]
means <- vector()
for (i in 1:length(colnames(means.data))) {
  means <- append(means, means.data[,i])
}
interacting_pair <- rep(rownames(means.data), length(colnames(means.data)))
cell_pair <- vector()
for (i in 1:length(colnames(means.data))) {
  cell_pair <- append(cell_pair, rep(colnames(means.data)[i], length(rownames(means.data))))
}
pvalues.data <- read.xlsx('R/sc.KNC_kJID_hNC/CellPhoneDB/9cellType_FB6cluster/FBc012-CellPhoneDB_bubblePlot.data/interactions_after3iterationsOfselectedKeloid/selectedKeloid_iteration3/onlyselectedKeloid_iteration3.FBc012-pvalues.xlsx', sheetIndex = 1)
row.names(pvalues.data) <- pvalues.data[,1]
pvalues.data <- pvalues.data[,-1]
pvalues <- vector()
for (i in 1:length(colnames(pvalues.data))) {
  pvalues <- append(pvalues, pvalues.data[,i])
}
pvalues[pvalues>0.05] <- 1
bubblePlot.data <- data.frame(interacting_pair=interacting_pair, cell_pair=cell_pair, means=means, pvalues=pvalues)
bubblePlot.data$cell_pair <- factor(bubblePlot.data$cell_pair, levels = c('FB_0.FB_0', 'FB_0.FB_1', 'FB_0.FB_2', 'FB_0.FB_3', 'FB_0.FB_4', 'FB_0.FB_5', 'FB_0.EC', 'FB_0.LEC', 'FB_0.LL', 'FB_0.MELA.NEU', 'FB_0.ML', 'FB_0.SGC', 'FB_0.SMC', 'FB_0.baKC', 'FB_0.spKC', 'FB_1.FB_0', 'FB_1.FB_1', 'FB_1.FB_2', 'FB_1.FB_3', 'FB_1.FB_4', 'FB_1.FB_5', 'FB_1.EC', 'FB_1.LEC', 'FB_1.LL', 'FB_1.MELA.NEU', 'FB_1.ML', 'FB_1.SGC', 'FB_1.SMC', 'FB_1.baKC', 'FB_1.spKC', 'FB_2.FB_0', 'FB_2.FB_1', 'FB_2.FB_2', 'FB_2.FB_3', 'FB_2.FB_4', 'FB_2.FB_5', 'FB_2.EC', 'FB_2.LEC', 'FB_2.LL', 'FB_2.MELA.NEU', 'FB_2.ML', 'FB_2.SGC', 'FB_2.SMC', 'FB_2.baKC', 'FB_2.spKC')) # change the order of axis labels
bubblePlot.data$interacting_pair <- factor(bubblePlot.data$interacting_pair, levels = c('FGFR3_EPHA4', 'FGFR3_NECTIN1', 'FGF7_FGFR3', 'FGF18_FGFR1', 'FGF18_FGFR3', 'EFNA1_EPHA2', 'EFNB2_EPHB3', 'EFNB2_EPHB4', 'EFNB2_EPHA4', 'EPHA4_EFNB1', 'EPHB3_EFNB1', 'EPHB4_EFNB1', 'NOTCH1_JAG1', 'NOTCH2_JAG2', 'NOTCH1_DLL4', 'NOTCH2_DLL4', 'DLL1_NOTCH1', 'DLL1_NOTCH2', 'DLL1_NOTCH3', 'DLL1_NOTCH4', 'DLL4_NOTCH3', 'LGALS9_CD47', 'LGALS9_LRP1', 'LGALS9_SLC1A5', 'LGALS9_PTPRK', 'LGALS9_CD44', 'LGALS9_MRC2', 'TNFRSF1B_GRN', 'TNFRSF10A_TNFSF10', 'TNFSF12_TNFRSF25', 'MIF_TNFRSF10D', 'ADRB2_VEGFB', 'VEGFA_KDR', 'FLT1 complex_VEGFA', 'FLT1 complex_VEGFB', 'NRP2_PGF', 'FLT1_PGF', 'FLT1 complex_PGF', 'PDGFB_PDGFRB', 'PDGFB_LRP1', 'ANGPT2_TEK', 'TGFB3_TGFBR3', 'COL5A3_a1b1 complex', 'COL5A3_a2b1 complex', 'COL7A1_a1b1 complex', 'COL7A1_a2b1 complex', 'COL7A1_a10b1 complex', 'COL17A1_a1b1 complex', 'COL17A1_a2b1 complex', 'COL17A1_a10b1 complex', 'COL21A1_a1b1 complex', 'COL21A1_a2b1 complex', 'COL27A1_a1b1 complex', 'COL27A1_a2b1 complex', 'DSG1_DSC2', 'DSG1_DSC3', 'CXADR_FAM3C', 'NECTIN1_NECTIN3', 'CDH1_a2b1 complex', 'SEMA4A_PLXND1', 'PLXNB2_SEMA4C', 'SELP_CD34', 'IL6 receptor_IL6', 'CSF1_SLC7A1', 'LRP6_CKLF', 'SIRPA_CD47')) # change the order of column labels
ggplot(bubblePlot.data, aes(x=cell_pair, y=interacting_pair)) + labs(caption = 'onlyselectedKeloid_iteration3') +
  geom_point(aes(color=means, size=-log10(pvalues+0.001))) +
  scale_color_gradient2(low = "#014aa8", high = "#98054f",  mid = "#fed9b6", midpoint = 0.7, limit = c(0, 1.4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_line(colour = 'lightgrey'))


#################################################


###################### FB_8 clusters ###########################
# ______Explore 8 clusters using 0.168 resolution
DefaultAssay(sc.kNC_kJID_hNC.FB) <- "integrated"
sc.kNC_kJID_hNC.FB_0.168 <- FindClusters(sc.kNC_kJID_hNC.FB, resolution = 0.168)
DimPlot(sc.kNC_kJID_hNC.FB_0.168, reduction = "umap", label = T, repel = T)
DimPlot(sc.kNC_kJID_hNC.FB_0.168, reduction = "umap", label = F, group.by = "zt.scarType", shuffle = T)
DimPlot(sc.kNC_kJID_hNC.FB_0.168, reduction = "umap", label = F, group.by = "zt.dataSource", shuffle = T)
DimPlot(sc.kNC_kJID_hNC.FB_0.168, reduction = "umap", label = F, group.by = "orig.ident", shuffle = T)
DimPlot(sc.kNC_kJID_hNC.FB_0.168, reduction = "umap", label = T, split.by = "zt.scarType")
save(sc.kNC_kJID_hNC.FB_0.168, file = "R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.168resolution_8clusters.Rdata")

# identify the DEGs of each cluster in FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.168resolution_8clusters.Rdata")
options(future.globals.maxSize= 3210612736)
DefaultAssay(sc.kNC_kJID_hNC.FB_0.168) <- "RNA"
# cluster 0 in 8 clusters
c0.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.168, ident.1 = 0, verbose = FALSE, logfc.threshold = 0.1)
write.table(c0.markers, file = "R/sc.KNC_kJID_hNC/FB/DEGs_8clusters/c0.markers", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 1 in 8 clusters
c1.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.168, ident.1 = 1, verbose = FALSE, logfc.threshold = 0.1)
write.table(c1.markers, file = "R/sc.KNC_kJID_hNC/FB/DEGs_8clusters/c1.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 2 in 8 clusters
c2.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.168, ident.1 = 2, verbose = FALSE, logfc.threshold = 0.1)
write.table(c2.markers, file = "R/sc.KNC_kJID_hNC/FB/DEGs_8clusters/c2.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 3 in 8 clusters
c3.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.168, ident.1 = 3, verbose = FALSE, logfc.threshold = 0.1)
write.table(c3.markers, file = "R/sc.KNC_kJID_hNC/FB/DEGs_8clusters/c3.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 4 in 8 clusters
c4.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.168, ident.1 = 4, verbose = FALSE, logfc.threshold = 0.1)
write.table(c4.markers, file = "R/sc.KNC_kJID_hNC/FB/DEGs_8clusters/c4.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 5 in 8 clusters
c5.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.168, ident.1 = 5, verbose = FALSE, logfc.threshold = 0.1)
write.table(c5.markers, file = "R/sc.KNC_kJID_hNC/FB/DEGs_8clusters/c5.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 6 in 8 clusters
c6.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.168, ident.1 = 6, verbose = FALSE, logfc.threshold = 0.1)
write.table(c6.markers, file = "R/sc.KNC_kJID_hNC/FB/DEGs_8clusters/c6.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
# cluster 7 in 8 clusters
c7.markers <- FindMarkers(sc.kNC_kJID_hNC.FB_0.168, ident.1 = 7, verbose = FALSE, logfc.threshold = 0.1)
write.table(c7.markers, file = "R/sc.KNC_kJID_hNC/FB/DEGs_8clusters/c7.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

## review top10 DEGs of each FBcluster in UMAP
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.168resolution_8clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.168) <- "RNA"
# top10 DEGs of FB cluster0
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("COL1A1", "CTHRC1", "ASPN", "POSTN", "SPARC", "COL1A2", "COL3A1", "ELN", "WISP2", "FN1"), ncol = 5)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("COL1A1", "CTHRC1", "ASPN", "POSTN", "SPARC"), split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("COL1A2", "COL3A1", "ELN", "WISP2", "FN1"), split.by = "zt.scarType")
# top10 DEGs of FB cluster1
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("CCL19", "APOE", "IGFBP7", "CXCL12", "C3", "PTGDS", "GGT5", "ABCA8", "TXNIP", "C7"), ncol = 5)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("CCL19", "APOE", "IGFBP7", "CXCL12", "C3"), split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("PTGDS", "GGT5", "ABCA8", "TXNIP", "C7"), split.by = "zt.scarType")
# top10 DEGs of FB cluster2
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("CXCL3", "CXCL1", "C11orf96", "CXCL2", "IL6", "GEM", "SOD2", "CYP1B1", "TNFAIP6", "HMOX1"), ncol = 5)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("CXCL3", "CXCL1", "C11orf96", "CXCL2", "IL6"), split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("GEM", "SOD2", "CYP1B1", "TNFAIP6", "HMOX1"), split.by = "zt.scarType")
# top10 DEGs of FB cluster3
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("APCDD1", "ID1", "F13A1", "COL18A1", "CD9", "COL6A5", "WIF1", "TCF4", "PTGDS", "NKD2"), ncol = 5)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("APCDD1", "ID1", "F13A1", "COL18A1", "CD9"), split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("COL6A5", "WIF1", "TCF4", "PTGDS", "NKD2"), split.by = "zt.scarType")
# top10 DEGs of FB cluster4
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("HBB", "HBA2", "HLA-DRA", "TM4SF1", "HBA1", "KRT14", "SELE", "CD74", "HLA-DRB1", "HLA-DPA1"), ncol = 5)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("HBB", "HBA2", "HLA-DRA", "TM4SF1", "HBA1"), split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("KRT14", "SELE", "CD74", "HLA-DRB1", "HLA-DPA1"), split.by = "zt.scarType")
# top10 DEGs of FB cluster5
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("LRRC38", "SCN7A", "IGFBP2", "RAMP1", "SEMA3G", "FGFBP2", "COL26A1", "NECAB1", "OLFML2A", "ALX4"), ncol = 5)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("LRRC38", "SCN7A", "IGFBP2", "RAMP1", "SEMA3G"), split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("FGFBP2", "COL26A1", "NECAB1", "OLFML2A", "ALX4"), split.by = "zt.scarType")
# top10 DEGs of FB cluster6
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("APOD", "ANGPTL7", "C2orf40", "TM4SF1", "DCD", "NR2F2", "CLDN1", "PTGDS", "ITGA6", "GPC3"), ncol = 5)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("APOD", "ANGPTL7", "C2orf40", "TM4SF1", "DCD"), split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("NR2F2", "CLDN1", "PTGDS", "ITGA6", "GPC3"), split.by = "zt.scarType")
# top10 DEGs of FB cluster7
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("IL24", "CSF2", "MMP3", "MMP1", "FCMR", "IL13RA2", "SERPINB2", "CXCL8", "CXCL5", "IL7R"), ncol = 5)
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("IL24", "CSF2", "MMP3", "MMP1", "FCMR"), split.by = "zt.scarType")
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("IL13RA2", "SERPINB2", "CXCL8", "CXCL5", "IL7R"), split.by = "zt.scarType")
# IL11
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("IL11"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.168, features = c("IL11"), split.by = "zt.scarType")

# APCDD1
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("APCDD1"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("APCDD1"), split.by = "zt.scarType")
# COL6A5
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("COL6A5"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("COL6A5"), split.by = "zt.scarType")
# PTGDS
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("PTGDS"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("PTGDS"), split.by = "zt.scarType")
# ENTPD1
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("ENTPD1"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("ENTPD1"), split.by = "zt.scarType")

#################################################


#################### explore a particular gene: selectedGenes (importantTF) #############################

## deep dive into subclustered FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"

sc.kNC_kJID_hNC.FB_0.16.c0 <- subset(sc.kNC_kJID_hNC.FB_0.16, idents='0')
selectedGenes <- c('ALX3', 'USF2', 'HOXC5', 'ATF2', 'ZNF76', 'ZFX', 'NRF1', 'FOXM1', 'HOXA3', 'CRTC2', 'CARF', 'ZNF319', 'MAX', 'NPAS2', 'FOXA1', 'WT1', 'FOXL1', 'ZNF562', 'ATF1', 'ZBTB33', 'TFF3', 'FOXD3', 'DLX5', 'ZNF823', 'USF1', 'UBE2V1', 'SOX12', 'BARX2', 'POU2F2', 'ZNF12', 'ZNF689', 'RFX3', 'ATF6B', 'ZNF639', 'TBP', 'VDR', 'HOXA4', 'ARID3A')
# expression in different scarTypeselectedGenes
VlnPlot(sc.kNC_kJID_hNC.FB_0.16.c0, features = selectedGenes, combine = T, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16.c0, features = selectedGenes, combine = T, group.by = "zt.scarType", pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.FB_0.16.c0, features = selectedGenes, dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
DotPlot(sc.kNC_kJID_hNC.FB_0.16.c0, features = selectedGenes, dot.scale = 8, group.by = "zt.FB6cluster_scarType") + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
# expression of selected TF, which is the intersection of 10 iterations of GeneRegulatoryNetwork, in FB6clusters_scarType
TF.commonHTSkeloid <- c('CREB3L1', 'MXD4', 'BHLHE40', 'JUND', 'TRIM28', 'CEBPB', 'FOSB', 'XBP1', 'CHD2', 'FOS', 'JUNB', 'NR2F2', 'TBX15', 'ATF3', 'MEF2C', 'REL')
TF.onlyHTS <- c('MAFG', 'NFIL3', 'NFKB1', 'ETS1', 'IRF1', 'MAFF', 'NFE2L2', 'PRDM1', 'CREM', 'ELF1', 'ETS2', 'EZH2', 'FOSL2', 'YY1')
TF.onlyKeloid <- c('CEBPD', 'EGR1', 'EGR2', 'EGR3', 'FOXD1', 'JUN', 'STAT3', 'TFAP2A', 'CREB5', 'IRF8')
TF <- c(TF.commonHTSkeloid, TF.onlyHTS, TF.onlyKeloid)
DotPlot(sc.kNC_kJID_hNC.FB_0.16, features = TF, dot.scale = 8, group.by = "zt.FB6cluster_scarType") + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
# expression in different subclusters
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = selectedGenes, combine = T)
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = selectedGenes, combine = T, pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.FB_0.16, features = selectedGenes, dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = selectedGenes, combine = T, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = selectedGenes, combine = T, split.by = "zt.scarType", pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.FB_0.16, features = selectedGenes, dot.scale = 8, group.by = "zt.FB6cluster_scarType") + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')



## Split complexHeatmap: FB6clusters
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
sc.kNC_kJID_hNC.FB_0.16.forHeatmap <- sc.kNC_kJID_hNC.FB_0.16
selectedGenes <- c('ALX3', 'USF2', 'HOXC5', 'ATF2', 'ZNF76', 'ZFX', 'NRF1', 'FOXM1', 'HOXA3', 'CRTC2', 'CARF', 'ZNF319', 'MAX', 'NPAS2', 'FOXA1', 'WT1', 'FOXL1', 'ZNF562', 'ATF1', 'ZBTB33', 'TFF3', 'FOXD3', 'DLX5', 'ZNF823', 'USF1', 'UBE2V1', 'SOX12', 'BARX2', 'POU2F2', 'ZNF12', 'ZNF689', 'RFX3', 'ATF6B', 'ZNF639', 'TBP', 'VDR', 'HOXA4', 'ARID3A')
sc.kNC_kJID_hNC.FB_0.16.forHeatmap <- sc.kNC_kJID_hNC.FB_0.16.forHeatmap[selectedGenes, ]
sc.kNC_kJID_hNC.FB_0.16.forHeatmap <- ScaleData(sc.kNC_kJID_hNC.FB_0.16.forHeatmap, features = rownames(sc.kNC_kJID_hNC.FB_0.16.forHeatmap))
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16.forHeatmap, slot = "scale.data")
dim(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.exprMatrix)
Idents(sc.kNC_kJID_hNC.FB_0.16.forHeatmap) <- "integrated_snn_res.0.16"
# extract the cluster information
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16.forHeatmap$integrated_snn_res.0.16)
# define the annotation labels of heatmap
cluster_annotation <-  HeatmapAnnotation(FB6cluster = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info,
                                         col = list(FB6cluster = c('0'="#e77d71", '1'="#b49f33", '2'="#64b74c", '3'="#53bbc2", '4'="#6b9cf8", '5'="#e46fdc")))
Heatmap(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.exprMatrix, name = 'FB',
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))



## Split complexHeatmap: FB6clusters
# HTS
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS <- subset(sc.kNC_kJID_hNC.FB_0.16, idents = 'HTS')
selectedGenes <- c('ALX3', 'USF2', 'HOXC5', 'ATF2', 'ZNF76', 'ZFX', 'NRF1', 'FOXM1', 'HOXA3', 'CRTC2', 'CARF', 'ZNF319', 'MAX', 'NPAS2', 'FOXA1', 'WT1', 'FOXL1', 'ZNF562', 'ATF1', 'ZBTB33', 'TFF3', 'FOXD3', 'DLX5', 'ZNF823', 'USF1', 'UBE2V1', 'SOX12', 'BARX2', 'POU2F2', 'ZNF12', 'ZNF689', 'RFX3', 'ATF6B', 'ZNF639', 'TBP', 'VDR', 'HOXA4', 'ARID3A')
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS <- sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS[selectedGenes, ]
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS <- ScaleData(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS, features = rownames(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS))
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS.exprMatrix <- GetAssayData(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS, slot = "scale.data")
dim(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS.exprMatrix)
Idents(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS) <- "integrated_snn_res.0.16"
# extract the cluster information
sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS.cluster_info <- sort(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS$integrated_snn_res.0.16)
# define the annotation labels of heatmap
cluster_annotation <-  HeatmapAnnotation(FB6cluster = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS.cluster_info,
                                         col = list(FB6cluster = c('0'="#e77d71", '1'="#b49f33", '2'="#64b74c", '3'="#53bbc2", '4'="#6b9cf8", '5'="#e46fdc")))


Heatmap(sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS.exprMatrix, name = 'HTS',
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = T, 
        column_split = sc.kNC_kJID_hNC.FB_0.16.forHeatmap.HTS.cluster_info, 
        top_annotation = cluster_annotation, 
        row_names_gp = gpar(fontsize = 5))



### deep dive into SDC2
## review SDC2 in all cell types
load("R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "RNA"
# expression distribution
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), split.by = "zt.scarType")
# expression distribution by kernel density estimation
plot_density(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), pal = "viridis")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
plot_density(subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "HTS"), features = c("SDC2"), pal = "viridis")
plot_density(subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid"), features = c("SDC2"), pal = "viridis")
plot_density(subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Normal"), features = c("SDC2"), pal = "viridis")
# expression in different scarType
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), combine = T, group.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), combine = T, group.by = "zt.scarType", pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), cols = c("blue", "red"), dot.scale = 8, group.by = "zt.scarType") + RotatedAxis()
# expression in different cellType
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), combine = T)
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), combine = T, pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), combine = T, split.by = "zt.scarType")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType_scarType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')

## deep dive into subseted FB
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
sc.kNC_kJID_hNC.integrated_2000_0.1.FB <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "FB")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1.FB) <- "zt.cellType_scarType"
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("SDC2"), split.by = "zt.cellType_scarType")
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("SDC2"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("SDC2"), combine = T)
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("SDC2"), combine = T, pt.size = 0)


#################################################


###################### 骆申英焦亡课题: GSDMD, GSDMB ###########################
## review GSDMD, GSDMB in all cell types
load("R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "RNA"
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("GSDMD", "GSDMB"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("GSDMD", "GSDMB"), split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("GSDMD", "GSDMB"), combine = T, group.by = "zt.scarType")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("GSDMD", "GSDMB"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("GSDMD", "GSDMB"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()+ scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType_scarType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("GSDMD"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("GSDMD"), cols = c("blue", "grey", "red"), dot.scale = 8, split.by = "zt.scarType") + RotatedAxis()

## deep dive into FB
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
sc.kNC_kJID_hNC.integrated_2000_0.1.FB <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "FB")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1.FB) <- "zt.cellType_scarType"
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("GSDMD"), split.by = "zt.cellType_scarType")
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("GSDMD"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("GSDMD", "GSDMB"), combine = T)
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("GSDMD", "GSDMB"), combine = T, pt.size = 0)
# DEGs between FB_Keloid and FB_Normal
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1.FB) <- "zt.cellType_scarType"
DEG.FB.KeloidvsNormal <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, ident.1 = "FB_Keloid", ident.2 = "FB_Normal", verbose = FALSE, logfc.threshold = 0.01)
write.table(DEG.FB.KeloidvsNormal, file = "R/sc.KNC_kJID_hNC/ShenyingLuo/DEG.FB.KeloidvsNormal.markers_FC0.01.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)

## deep dive into ML
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
sc.kNC_kJID_hNC.integrated_2000_0.1.ML <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "ML")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1.ML) <- "zt.cellType_scarType"
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1.ML, features = c("GSDMD"), split.by = "zt.cellType_scarType")
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.ML, features = c("GSDMD"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.ML, features = c("GSDMD", "GSDMB"), combine = T)
# DEGs between ML_Keloid and ML_Normal
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1.ML) <- "zt.cellType_scarType"
DEG.ML.KeloidvsNormal <- FindMarkers(sc.kNC_kJID_hNC.integrated_2000_0.1.ML, ident.1 = "ML_Keloid", ident.2 = "ML_Normal", verbose = FALSE, logfc.threshold = 0.1)
write.table(DEG.ML.KeloidvsNormal, file = "R/sc.KNC_kJID_hNC/ShenyingLuo/DEG.ML.KeloidvsNormal.markers.xlsx", quote = F, sep = "\t", row.names = T, col.names = T)
#################################################

#################### explore a particular gene: SDC2 #############################
### deep dive into SDC2
## review SDC2 in all cell types
load("R/sc.kNC_kJID_hNC.integrated2000_0.1resolution_10clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "RNA"
# expression distribution
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"))
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), split.by = "zt.scarType")
# expression distribution by kernel density estimation
plot_density(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), pal = "viridis")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.scarType"
plot_density(subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "HTS"), features = c("SDC2"), pal = "viridis")
plot_density(subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Keloid"), features = c("SDC2"), pal = "viridis")
plot_density(subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "Normal"), features = c("SDC2"), pal = "viridis")
# expression in different scarType
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), combine = T, group.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), combine = T, group.by = "zt.scarType", pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), cols = c("blue", "red"), dot.scale = 8, group.by = "zt.scarType") + RotatedAxis()
# expression in different cellType
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), combine = T)
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), combine = T, pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), combine = T, split.by = "zt.scarType")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType_scarType"
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1, features = c("SDC2"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')

## deep dive into subseted FB
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1) <- "zt.cellType"
sc.kNC_kJID_hNC.integrated_2000_0.1.FB <- subset(sc.kNC_kJID_hNC.integrated_2000_0.1, idents = "FB")
Idents(sc.kNC_kJID_hNC.integrated_2000_0.1.FB) <- "zt.cellType_scarType"
FeaturePlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("SDC2"), split.by = "zt.cellType_scarType")
DotPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("SDC2"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("SDC2"), combine = T)
VlnPlot(sc.kNC_kJID_hNC.integrated_2000_0.1.FB, features = c("SDC2"), combine = T, pt.size = 0)

## deep dive into subclustered FB
load("R/sc.KNC_kJID_hNC/FB/sc.kNC_kJID_hNC.FB_0.16resolution_6clusters.Rdata")
DefaultAssay(sc.kNC_kJID_hNC.FB_0.16) <- "RNA"
# expression distribution
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"))
FeaturePlot(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"), split.by = "zt.scarType")
# expression distribution by kernel density estimation
plot_density(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"), pal = "viridis")
Idents(sc.kNC_kJID_hNC.FB_0.16) <- "zt.scarType"
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16, idents = "HTS"), features = c("SDC2"), pal = "viridis")
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16, idents = "Keloid"), features = c("SDC2"), pal = "viridis")
plot_density(subset(sc.kNC_kJID_hNC.FB_0.16, idents = "Normal"), features = c("SDC2"), pal = "viridis")
# expression in different scarType
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"), combine = T, group.by = "zt.cellType_scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"), combine = T, group.by = "zt.cellType_scarType", pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"), dot.scale = 8, group.by = "zt.cellType_scarType") + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
# expression in different subclusters
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"), combine = T)
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"), combine = T, pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"), dot.scale = 8) + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"), combine = T, split.by = "zt.scarType")
VlnPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"), combine = T, split.by = "zt.scarType", pt.size = 0)
DotPlot(sc.kNC_kJID_hNC.FB_0.16, features = c("SDC2"), dot.scale = 8, group.by = "zt.FB6cluster_scarType") + RotatedAxis() + scale_color_gradient2(low = "blue", mid = 'grey', high = 'red')
#################################################

