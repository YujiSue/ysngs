library(Seurat)
library(ggplot2)
library(dplyr)

# Load data
loadH5Matrix <- function(path, # path-to-matrix
                    title = "Seurat Analysis",
                    min_cells = 3,
                    min_features = 200) 
{
    counts <- Read10X_h5(path)
    seurat_obj <- CreateSeuratObject(counts = counts, 
                                      project = title, 
                                      min.cells = min_cells, 
                                      min.features = min_features)
}

# Filter cells
additionalFilter <- function(seurat_obj,
                    mt_prefix = "MT-",
                    min_features = 500, 
                    max_features = 6000, 
                    max_mt = 10) 
{
    ## Calculate percentage of mitochondrial genes 
    mt_pattern <- paste0("^", mt_prefix, separator = "")
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
    ## Filter cells
    seurat_obj_sub <- subset(seurat_obj, 
                             subset = nFeature_RNA > min_features & 
                                      nFeature_RNA < max_features & 
                                      percent.mt < max_mt)
    return(seurat_obj_sub)
}
# Export QC metrics
qcMetrics <- function(seurat_obj, out_dir=".", show_plot=TRUE, export_image=TRUE) 
{
    qc_data <- data.frame(
        nFeature_RNA = seurat_obj$nFeature_RNA,
        nCount_RNA = seurat_obj$nCount_RNA,
        percent.mt = seurat_obj$percent.mt
    )
    ## Export QC metrics to CSV
    write.csv(qc_data, file = file.path(out_dir, "qc_metrics.csv"), row.names = FALSE)
    ## Visualize QC metrics
    if (show_plot) {
        p <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        print(p)
        if (export_image) {
            ggsave(file.path(out_dir, "qc_violin_plot.png"), plot = p, width = 10, height = 4, dpi = 150)
        }
    }
}
# Normalization & select HVGs
hvgSelection <- function(seurat_obj, 
                          normalization_method = "LogNormalize", # "LogNormalize" or "CLR" or "RC" or "SCT"
                          scale = 10000,
                          selection_method = "vst", # "vst" or "mean.var.plot" or "dispersion"
                          count = 3000) 
{
    seurat_obj_norm <- NormalizeData(seurat_obj, normalization.method = normalization_method, scale.factor = scale)
    seurat_obj_norm <- FindVariableFeatures(seurat_obj_norm, selection.method = selection_method, nfeatures = count)
    return(seurat_obj_norm)
}

# Perform PCA
seuratPCA <- function(seurat_obj,
                       dims = 50) 
{
    seurat_obj_scaled <- ScaleData(seurat_obj)
    seurat_obj_scaled <- RunPCA(seurat_obj_scaled, npcs = dims, verbose = FALSE)
    return(seurat_obj_scaled)
}

# Clustering
scCluster <- function(seurat_obj,
                      dims = 1:20,
                      cluster_res = 0.5,
                      harmonized = FALSE,
                      sample_ids = "",
                      umap_reduction = TRUE,
                      umap_name = "umap",
                      show_plot = TRUE,
                      out_dir = ".",
                      export_image = TRUE) {
    #
    if (harmonized) {
        seurat_obj_classified <- FindNeighbors(seurat_obj, reduction = "harmony", dims = dims)
        seurat_obj_classified <- FindClusters(seurat_obj_classified, resolution = cluster_res, cluster.name = "clusters")
        seurat_obj_classified <- RunUMAP(seurat_obj_classified, reduction = "harmony", dims = dims, reduction.name = umap_name)
    }
    else {
        seurat_obj_classified <- FindNeighbors(seurat_obj, dims = dims)
        seurat_obj_classified <- FindClusters(seurat_obj_classified, resolution = cluster_res)
    }
    if (umap_reduction) seurat_obj_classified <- RunUMAP(seurat_obj_classified, dims = dims, reduction.name = umap_name)
    if (show_plot) {
        if (harmonized) p <- DimPlot(seurat_obj, reduction = "umap.harmony", group.by = c("clusters", sample_ids), label = TRUE)
        else p <- DimPlot(seurat_obj, reduction = umap_name, label = TRUE)
        print(p)
        if (export_image) {
            ggsave(file.path(out_dir, "umap_plot.png"), plot = p, width = 8, height = 6, dpi = 150)
        }
    }
    return(seurat_obj_classified)
}

# Select markers for each cluster and visualize top markers
selectMarkers <- function(seurat_obj, 
                      only_pos = TRUE, 
                      min_pct = 0.25, 
                      threshold = 0.25,
                      show_plot = TRUE,
                      show_count = 10,
                      out_dir = ".", 
                      export_image = TRUE) 
{
    markers <- FindAllMarkers(seurat_obj, only.pos = only_pos, min.pct = min_pct, logfc.threshold = threshold)
    if (show_plot) {
        top_genes <- markers %>% group_by(cluster) %>% top_n(n = show_count, wt = avg_log2FC)
        p <- DoHeatmap(seurat_obj, features = top_genes$gene) + NoLegend()
        print(p)
        if (export_image) {
            ggsave(file.path(out_dir, paste("top", show_count, "_markers_heatmap.png", sep = "")), plot = p, width = 8, height = 6, dpi = 150)
        }
    }
    return(markers)
}

# Dot plot original markers
plotMarkers <- function(seurat_obj, 
                        markers,
                        marker_title = "custom_marker_set",
                        out_dir = ".", 
                        export_image = TRUE) 
{
    p <- DotPlot(seurat_obj, features = markers) + RotatedAxis()
    print(p)
    ggsave(file.path(out_dir, paste(marker_title, "_plot.png", sep = "")), plot = p, width = 6, height = 5, dpi = 150)
}

# Correcting batch effects
correctBatch <- function(seurat_obj, sample_ids = "orig.ident") {
    library(harmony)
    seurat_obj_corrected <- RunHarmony(seurat_obj, group.by.vars = sample_ids)
    return(seurat_obj_corrected)
} 

# Rename clusters
setClusterName <- function(seurat_obj, names) {
    mycluster.ids <- names
    names(mycluster.ids) <- levels(seurat_obj)
    seurat_obj_renamed <- RenameIdents(seurat_obj, mycluster.ids)
    return(seurat_obj_renamed)
}

# Cell annotation
cellAnnotation <- function(seurat_obj, 
                            method = "azimuth", # "azimuth" or "singler"
                            ref_db = "pbmcref", # Azimuth reference database
                            # ex. ref_db <- celldex::HumanPrimaryCellAtlasData() for SingleR
                            out_dir=".", 
                            show_plot=TRUE, 
                            export_image=TRUE) 
{
    if (method == "azimuth") {
        ## Load Azimuth reference database
        library(Azimuth)
        ## Run Azimuth for cell annotation
        seurat_obj_annot <- RunAzimuth(query = seurat_obj, reference = ref_db)
        ## Visualize the results
        p <- DimPlot(seurat_obj_annot, 
                    reduction = "ref.umap", 
                    group.by = "predicted.celltype.l2", 
                    label = TRUE) +
            NoLegend() + 
            ggtitle("Cell Annotation (Azimuth Reference Space)")
        if (show_plot) print(p)
        if (export_image) {
            ggsave(file.path(out_dir, "cell_annotation_azimuth.png"), plot = p, width = 8, height = 6, dpi = 150)
        }
        return(seurat_obj_annot)
    } else if (method == "singler") {
        ## Load SingleR and celldex libraries
        library(SingleR)
        library(celldex)
        ## Extract counts matrix from Seurat object
        counts <- GetAssayData(seurat_obj, slot = "data")
        prediced <- SingleR(test = counts, ref = ref_db, labels = ref_db$label.main)
        seurat_obj_annot <- seurat_obj
        seurat_obj_annot$predicted.singler <- prediced$labels
        ## Visualize the results
        p <- DimPlot(seurat_obj_annot, 
                    reduction = "umap", 
                    group.by = "predicted.singler", 
                    label = TRUE) +
            NoLegend() + 
            ggtitle("Cell Annotation (SingleR Reference Space)")
        if (show_plot) print(p)
        if (export_image) {
            ggsave(file.path(out_dir, "cell_annotation_singler.png"), plot = p, width = 8, height = 6, dpi = 150)
        }
        return(seurat_obj_annot)
    } else {
        stop("Invalid method. Choose 'azimuth' or 'singler'.")
    }
}
