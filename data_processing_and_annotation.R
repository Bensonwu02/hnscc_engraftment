### Data processing and annotation of Choi et al. (https://doi.org/10.1038/s41467-023-36691-x) and Kurten et al. (https://doi.org/10.1038/s41467-021-27619-4) datasets

library(dplyr)
library(ggplot2)
library(Seurat)
library(data.table)
library(stringr)
library(leiden)
library(glmGamPoi)
library(reticulate)
library(infercnv)

out_dir <- "path_to_output_directory"

plot_dir <- paste0(out_dir, "/processing_and_annotation/")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---------------------------------------------------------------------------- #
#                                 Choi dataset                                 #
# ---------------------------------------------------------------------------- #

# This dataset is already annotated. Simply load in raw counts and metadata to create a Seurat object

umi_counts <- as.data.frame(fread("/GSE181919_UMI_counts.txt"))
rownames(umi_counts) <- umi_counts$V1
umi_counts <- umi_counts[,-1]

metadata <- read.table("/GSE181919_Barcode_metadata.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rownames(metadata) <- gsub("-", ".", rownames(metadata))
metadata$Malignant <- ifelse(metadata$cell.type == "Malignant.cells", TRUE, FALSE)
metadata$cell.type <- case_when(
    metadata$cell.type == "B_Plasma.cells" ~ "B_cell",
    metadata$cell.type == "Dendritic.cells" ~ "Dendritic_cell",
    metadata$cell.type == "Endothelial.cells" ~ "Endothelial",
    metadata$cell.type == "Epithelial.cells" ~ "Epithelial",
    metadata$cell.type == "Fibroblasts" ~ "Fibroblast",
    metadata$cell.type == "Macrophages" ~ "Macrophage",
    metadata$cell.type == "Malignant.cells" ~ "Epithelial",
    metadata$cell.type == "Mast.cells" ~ "Mast",
    metadata$cell.type == "Myocytes" ~ "Myocyte",
    metadata$cell.type == "T.cells" ~ "T_cell",
)

choi <- CreateSeuratObject(umi_counts, project = "Choi")
choi <- AddMetaData(choi, metadata = metadata)
choi$dataset <- "Choi"
choi <- NormalizeData(choi)
print(choi)
print(table(choi$hpv))

# ---------------------------------------------------------------------------- #
#                                Kurten dataset                                #
# ---------------------------------------------------------------------------- #

# This dataset only provides cellranger output.
# We will need to annotate these samples. I will follow their processing and annotation workflow as close as possible
hpv_neg_patients <- c("HN01", "HN02", "HN03", "HN04", "HN05", "HN06", "HN07", "HN08", "HN09", "HN10", "HN11", "HN15")

path_to_dirs <- "/Kurten-Data/"

dirs <- list.files(path_to_dirs, recursive = FALSE, full.names = TRUE)
list_of_so <- list()

for (dir in dirs) { # Create object for each sample
    expr_mtx <- Read10X(dir)
    so <- CreateSeuratObject(expr_mtx, project = "Kurten")
    print(so)

    sample.id <- basename(dir)
    print(paste0("sample id: ", sample.id))
    so$sample.id <- sample.id

    patient.id <- str_split(basename(dir), "_")[[1]][1]
    print(paste0("patient id: ", patient.id))
    so$patient.id <- patient.id

    tissue.type <- str_split(basename(dir), "_")[[1]][2]
    print(paste0("tissue type: ", tissue.type))
    so$tissue.type <- tissue.type

    hpv.status <- ifelse(patient.id %in% hpv_neg_patients, "HPV-", "HPV+")
    print(paste0("HPV status: ", hpv.status))
    so$hpv <- hpv.status

    print(so)

    ### QC, following steps reported in 'Quality control (QC) and filtering the data set' section of methods 
    so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

    # Filtering
    so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

    # Filter genes
    so <- subset(so, features = rowSums(so@assays$RNA@counts > 0) >= 3)

    list_of_so[[sample.id]] <- so

}

kurten <- merge(list_of_so[[1]], list_of_so[-1]) # Merge all samples into one object
kurten$dataset <- "Kurten"
print(kurten)

# We are interested in epithelial and fibroblast population. I keep only CD45neg samples because they contain these populations. I also keep PBL for inferCNV ref.
kurten <- subset(kurten, subset = tissue.type %in% c("CD45n", "PBL"))

# ------------------------- Normalization, dim reduc ------------------------- #
# follow processing steps reported in 'Normalization, dimensionality reduction, and data visualization' section of methods

kurten <- NormalizeData(kurten, normalization.method = "LogNormalize", scale.factor = 10000)
kurten <- FindVariableFeatures(kurten, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
kurten <- ScaleData(kurten, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE)
kurten <- RunPCA(kurten, verbose = FALSE)

# Using the first 10 PCs, as per methods
kurten <- RunUMAP(kurten, dims = 1:10)
kurten <- FindNeighbors(kurten, dims = 1:10)

# ----------------------- Clustering, cell type assign ----------------------- #

# Clustering
kurten <- FindClusters(kurten, resolution = 0.3, method = "igraph", algorithm = 4)

# UMAP visualization
p <- DimPlot(kurten, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(filename = "kurten_umap.pdf", path = plot_dir)

saveRDS(kurten, "/kurten.rds")

# Cell type assignment using marker genes from Supp. Fig1A
markers <- c("CD3D", "CD8A", "CD4", "IL7R", "FOXP3", "CD19", "CD79A", "CD79B", "IL3RA", "NCAM1", "NCR1", "KLRD1",
           "IL5RA", "CSF1R", "MS4A7", "CD14", "MS4A7", "CD86", "CD83", "LAMP3", "FCGRA1", "FPR1", "THBD", "CD80",
           "CD86", "CD1C", "CD209", "LYZ", "TPSB2", "TPSAB1",
           "CLDN5", "FLT1", "CDH5", "RAMP2", "EPCAM", "KRT14", "KRT17", "COL1A1", "COL1A2", "DCN", "COL1A2", "RGS5")
markers <- unique(markers[markers %in% rownames(kurten)])

kurten$seurat_clusters <- factor(kurten$seurat_clusters, levels = unique(kurten$seurat_clusters))
p <- DotPlot(kurten, features = markers) +
    theme(axis.text.x = element_text(angle = 90))
ggsave(filename = "kurten_marker_genes_by_cluster.pdf", path = plot_dir)

# Assigning cell types based on marker genes 
kurten$cell.type <- case_when(
    kurten$seurat_clusters %in% c("1", "3", "12", "14") ~ "T_cell",
    kurten$seurat_clusters == "7" ~ "B_cell",
    kurten$seurat_clusters == "2" ~ "Monocyte",
    kurten$seurat_clusters == "11" ~ "Macrophage",
    kurten$seurat_clusters %in% c("4", "5", "10") ~ "Epithelial",
    kurten$seurat_clusters %in% c("6", "9") ~ "Endothelial",
    kurten$seurat_clusters == "13" ~ "Pericyte",
    kurten$seurat_clusters == "8" ~ "Fibroblast"
)
p <- DimPlot(kurten, reduction = "umap", group.by = "cell.type", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(filename = "kurten_umap_annot.pdf", path = plot_dir)

p <- DotPlot(kurten, features = markers, group.by = "cell.type") +
    theme(axis.text.x = element_text(angle = 90))
ggsave(filename = "kurten_marker_genes_by_cell_type.pdf", path = plot_dir)

# --------------------------------- infercnv --------------------------------- #

# Load gene ordering file
gene_order <- "/gencode_v21_gen_pos.complete.txt"

kurten_infercnv_epi <- subset(kurten, subset = cell.type == "Epithelial")
set.seed(123)
kurten_infercnv_pbl <- subset(kurten, subset = tissue.type == "PBL")
kurten_infercnv_pbl <- subset(kurten_infercnv_pbl, cells = sample(Cells(kurten_infercnv_pbl), 1000)) # Sample 1000 PBLs to use as reference to reduce computational load
kurten_infercnv <- merge(kurten_infercnv_epi, kurten_infercnv_pbl)

# Get raw counts matrix
counts_matrix <- GetAssayData(kurten_infercnv, slot = "counts", assay = "RNA") # Get raw counts matrix from Seurat object

# Set up sample annotation file
metadata <- kurten_infercnv[[]] %>% select(c("cell.type", "tissue.type"))
metadata$annotations <- ifelse(metadata$cell.type == "Epithelial", "Epithelial", "PBL")
cell_annotations <- metadata %>% select("annotations")

# Set reference cells
ref_names <- c("PBL")

cell_annotations_file_path <- "/kurten_epithelial_cells_bg_PBL_downsample.txt"

write.table(cell_annotations, cell_annotations_file_path, col.names = FALSE, row.names = TRUE, sep = "\t")

# ------------------------------- Run inferCNV ------------------------------- #
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                    annotations_file = cell_annotations_file_path,
                                    delim = "\t",
                                    gene_order_file = gene_order,
                                    ref_group_names = ref_names)

infercnv_output_dir <- "/kurten_epithelial_cells_bg_PBL_downsample"

options(scipen = 100)

leiden_res <- 0.0005
tumour_subclust_pval <- 0.05
bayesmaxp <- 0.25

infercnv_obj <- infercnv::run(infercnv_obj,
    cutoff = 0.1,
    out_dir = infercnv_output_dir,  # output directory
    num_threads = 10,

    cluster_by_groups = TRUE,
    denoise = TRUE,
    sd_amplifier = 1.5,
    HMM = TRUE, # TRUE to get CNV predictions. For visualization purposes, FALSE is fine.

    analysis_mode = "subclusters", # Usually 'subclusters'. For visualization purposes, samples can be better. 
    tumor_subcluster_partition_method = "leiden",
    leiden_resolution = leiden_res, # Default: 0.05. 
    tumor_subcluster_pval = tumour_subclust_pval, # Default: 0.1.
    BayesMaxPNormal = bayesmaxp,

    plot_probabilities = FALSE,
    save_rds = TRUE,
    diagnostics = FALSE,
    output_format = "png",
    plot_chr_scale = TRUE,
    write_phylo = FALSE
)

kurten_infercnv <- add_to_seurat(seurat_obj = kurten_infercnv, infercnv_output_path = infercnv_output_dir)

saveRDS(kurten_infercnv, "/kurten_infercnv_downsample.rds")

# ------------------ Assign malignant cells on basis of CNVs ----------------- #

kurten_infercnv <- readRDS("/kurten_infercnv_downsample.rds")
kurten_infercnv$Malignant <- ifelse(kurten_infercnv$has_cnv_chr9 == FALSE & kurten_infercnv$has_cnv_chr11 == FALSE,
    FALSE,
    TRUE)

malignant_cells <- WhichCells(kurten_infercnv, expression = Malignant == TRUE)

kurten$Malignant <- FALSE
kurten$Malignant[colnames(kurten) %in% malignant_cells] <- TRUE

p <- DimPlot(kurten, reduction = "umap", group.by = "Malignant")
ggsave(filename = "kurten_malignant_cells.pdf", path = plot_dir)

### Combine datasets
print(kurten)
print(choi)
hnscc <- merge(kurten, choi)
print(hnscc)
saveRDS(hnscc, "/hnscc.rds")

print("Script successfully completed")

sessionInfo()