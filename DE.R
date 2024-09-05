### Differential gene expression analysis 

library(dplyr)
library(ggplot2)
library(Seurat)
library(data.table)
library(stringr)
library(DESeq2)
library(SingleCellExperiment)
library(Matrix.utils)
library(patchwork)
library(ggrepel)
library(ggthemes)

out_dir <- "path_to_output_directory"

plot_dir <- paste0(out_dir, "/DE/")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# load data
so <- readRDS("hnscc.rds")

print(table(so$orig.ident, so$cell.type))


# Exclude HPV+ patients, keep only HPV-
so <- subset(so, subset = hpv == "HPV-")

print(table(so$orig.ident, so$cell.type))

# Only HNSCC samples 
so <- subset(so, subset = tissue.type %in% c("CD45n", "CA"))

print(table(so$orig.ident, so$cell.type))

# create cell.type.l2 annot column that separates out malignant epithelial with non-malignant epithelial
so$cell.type.l2 <- ifelse(so$Malignant == TRUE, "Malignant", so$cell.type)
print(table(so$orig.ident, so$cell.type.l2))

# ---------------------------------------------------------------------------- #
#                          function for pseudobulk DE                          #
# ---------------------------------------------------------------------------- #

#' Add annotated subset of non-malignant cells back to original Seurat object
#' 
#' @param seurat_object Seurat object containing cells that you want to perform differential expression testing on
#' @param sample_column Metadata column that contains sample id
#' @param annotation_column Metadata column that contains cell type annotations
#' @param cell_types Cell types you want to include in the test
#' @param contrast Contrast to define groups to test i.e., contrast = c("Condition of interest", "Group A", "Group B"). Group A will have positive LFC in results 
#' @param latent_vars Latent variables to include in design e.g., c("Patient", "Study")
#' @param lcfshrink Bool. TRUE to perform LFC shrinkage
#' @param out_dir Path to output directory
#' 
#' @return Data frame with DE results
perform_pseudobulk_DE <- function(seurat_object, sample_column, annotation_column, cell_types, contrast, latent_vars = NULL, lfcshrink, out_dir) {
    require(DESeq2)
    require(dplyr)
    require(SingleCellExperiment)

    # Create output dir if does not exist
    if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    # Subset Seurat object for cell types of interest
    Idents(seurat_object) <- annotation_column 
    seurat_object <- subset(seurat_object, idents = cell_types, invert = FALSE)
    print(table(seurat_object[[annotation_column]]))

    # Make pseudobulking column
    seurat_object$pseudobulk_group <- paste0(seurat_object@meta.data[[sample_column]], "_", seurat_object@meta.data[[contrast[1]]])
    print(table(seurat_object$pseudobulk_group))
        
    # Convert to SingleCellExperiment
    print("Converting to SCE")
    counts <- seurat_object@assays$RNA@counts
    metadata <- seurat_object[[]]
    sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
    
    # Aggregate counts by sample
    print("Pseudobulking")
    groups <- colData(sce)[, "pseudobulk_group"]
    aggr_counts <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum")
    aggr_counts <- t(aggr_counts)
    
    # Prepare metadata
    print("Preparing metadata")
    metadata <- colData(sce) %>%
        as.data.frame() %>%
        select(c("pseudobulk_group", contrast[1], latent_vars)) %>%
        distinct()
    rownames(metadata) <- metadata[["pseudobulk_group"]]
    t <- table(colData(sce)[["pseudobulk_group"]])
    df <- data.frame(pseudobulk_group = names(t), cell_count = as.integer(t))
    df <- plyr::join(df, metadata, by = "pseudobulk_group")
    df <- df[match(colnames(aggr_counts), df[["pseudobulk_group"]]), ]
    
    write.csv(df, file.path(out_dir, "samples_tested.csv"))
    
    # Convert col metadata to factors
    df[[contrast[1]]] <- factor(df[[contrast[1]]], levels = unique(df[[contrast[1]]]))
    for (var in latent_vars) {
        df[[var]] <- factor(df[[var]], levels = unique(df[[var]]))
    }
    
    # Create DESeq2 object
    if (!is.null(latent_vars)) {
        deseq <- DESeqDataSetFromMatrix(aggr_counts, colData = df, design = as.formula(paste0("~ ", paste(latent_vars, collapse = "+"), "+", contrast[1])))
    } else {
        deseq <- DESeqDataSetFromMatrix(aggr_counts, colData = df, design = as.formula(paste0("~ ", contrast[1])))
    }
    
    # Filter low count genes
    print("Filtering genes")
    num_in_groups <- df %>% group_by(get(contrast[1])) %>% summarise(n = n())
    print(num_in_groups)
    keep <- rowSums(counts(deseq) >= 3) >= min(num_in_groups$n)
    deseq <- deseq[keep, ]
    
    genes_used <- rownames(deseq)
    write.csv(data.frame(gene = genes_used), file.path(out_dir, "genes_used.csv"))
    
    # PCA
    print("Plotting PCA")
    rld <- rlog(deseq, blind = TRUE)
    pca_plots <- lapply(c(contrast[1], latent_vars), function(x) {
        p <- DESeq2::plotPCA(rld, ntop = 1000, intgroup = x)
        return(p)
    })
    plots <- wrap_plots(pca_plots)
    ggsave(plots, filename = "PCA_plot.pdf", path = out_dir)
    
    # Run DE analysis
    print("Performing DE")
    deseq <- DESeq(deseq)
    pdf(file = paste0(out_dir, "/deseq_disp_ests.pdf"), width = 10, height = 10)
    plotDispEsts(deseq)
    dev.off()
    
    res <- results(deseq, contrast = contrast, alpha = 0.05)

    if (lfcshrink) {
        res <- lfcShrink(deseq, contrast = contrast, type = "ashr", res = res)
    }
    
    res_tbl <- as_tibble(data.frame(gene = rownames(res), res))
    
    # Write all results to file
    print("Writing results")
    write.csv(res_tbl, file.path(out_dir, "de_results.csv"), quote = FALSE, row.names = FALSE)
    
    return(res_tbl)
}

volcano <- function(de_results, fc_threshold = 1, padj_threshold = 0.05, out_dir) {

    de_results <- as.data.frame(de_results)
    
    # Thresholds
    fc <- fc_threshold
    p_valadj <- padj_threshold

    de_results$diffexpressed <- "NO"
    de_results$diffexpressed[de_results$log2FoldChange > fc & de_results$padj < p_valadj] <- "YES"
    de_results$diffexpressed[de_results$log2FoldChange < -fc & de_results$padj < p_valadj] <- "YES"
    de_results$rank <- rank(de_results$padj)

    rownames(de_results) <- de_results$gene

    de_results$diffexpressed <- factor(de_results$diffexpressed, levels = c("NO", "YES"))
    de_results$label <- ifelse(de_results$diffexpressed == "YES", rownames(de_results), NA)

    p <- ggplot(de_results, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(colour = diffexpressed), alpha = 0.05, size = 3.5, stroke = NA) +
    scale_colour_manual(values = c("black", "red4")) +
    geom_label_repel(aes(label = label), size = 3, colour = "darkorange", na.rm = TRUE) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_vline(xintercept = -1, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    (theme_foundation(base_size = 14, base_family = "Helvetica") 
        + theme(
                plot.title = element_text(
                    face = "bold",
                    size = rel(1.2), hjust = 0.5
                ),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold", size = rel(1.1)),
                axis.title.y = element_text(angle = 90, vjust = 2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(1)),
                axis.line = element_line(colour = "black"),
                axis.ticks = element_line(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.key.size = unit(0.2, "cm"),
                legend.spacing = unit(0, "cm"),
                legend.title = element_text(face = "plain"),
                plot.margin = unit(c(10, 5, 5, 5), "mm"),
                strip.background = element_rect(colour = "#ffffff", fill = "#ffffff"),
                strip.text = element_text(face = "bold", size = rel(1.1))
            )) +
    xlab(bquote(~Log[2] ~ FoldChange)) +
    ylab(bquote(~-Log[10] ~ italic(FDR)))
    ggsave(filename = "volcano.pdf", path = out_dir, height = 10, width = 7)
}

# ---------------------------------------------------------------------------- #
#                     1: Fibroblast vs malignant epithelial                    #
# ---------------------------------------------------------------------------- #

res <- perform_pseudobulk_DE(so, 
    sample_column = "sample.id", 
    annotation_column = "cell.type.l2",
    cell_types = c("Fibroblast", "Malignant"),
    contrast = c("cell.type.l2", "Fibroblast", "Malignant"),
    latent_vars = c("orig.ident"),
    lfcshrink = TRUE,
    out_dir = file.path(plot_dir, "fibro_vs_malignant"))

print(head(res))
write.csv(res, file.path(plot_dir, "fibro_vs_malignant.csv"))

volcano(res, out_dir = file.path(plot_dir, "fibro_vs_malignant"))

# ---------------------------------------------------------------------------- #
#                   2: Fibroblast vs non-malignant epithelial                  #
# ---------------------------------------------------------------------------- #

res <- perform_pseudobulk_DE(so, 
    sample_column = "sample.id", 
    annotation_column = "cell.type.l2",
    cell_types = c("Fibroblast", "Epithelial"),
    contrast = c("cell.type.l2", "Fibroblast", "Epithelial"),
    latent_vars = c("orig.ident"),
    lfcshrink = TRUE,
    out_dir = file.path(plot_dir, "fibro_vs_epi"))

print(head(res))
write.csv(res, file.path(plot_dir, "fibro_vs_epi.csv"))

volcano(res, out_dir = file.path(plot_dir, "fibro_vs_epi"))

# ---------------------------------------------------------------------------- #
#              3: Malignant epithelial vs non-malignant epithelial             #
# ---------------------------------------------------------------------------- #

res <- perform_pseudobulk_DE(so, 
    sample_column = "sample.id", 
    annotation_column = "cell.type.l2",
    cell_types = c("Malignant", "Epithelial"),
    contrast = c("cell.type.l2", "Malignant", "Epithelial"),
    latent_vars = c("orig.ident"),
    lfcshrink = TRUE,
    out_dir = file.path(plot_dir, "malignant_vs_epi"))

print(head(res))
write.csv(res, file.path(plot_dir, "malignant_vs_epi.csv"))

volcano(res, out_dir = file.path(plot_dir, "malignant_vs_epi"))


# ---------------------------------------------------------------------------- #
#                4: Malignant epithelial LAM2C high vs LAM2C low               #
# ---------------------------------------------------------------------------- #

# Only looking at malignant cells now
malignant <- subset(so, subset = cell.type.l2 == "Malignant")

### Annotate cells as LAMC2 positive or negative
lamc2_expression <- FetchData(malignant, vars = "LAMC2")

threshold <- 0

# Annotate malignant cells as positive or negative
malignant$LAMC2_status_pos_neg <- case_when(
  lamc2_expression$LAMC2 > threshold ~ "Positive",
  lamc2_expression$LAMC2 <= threshold ~ "Negative",
  TRUE ~ "NA"
)

# Verify the annotation
table(malignant$LAMC2_status_pos_neg)

# DE
res <- perform_pseudobulk_DE(malignant, 
    sample_column = "sample.id", 
    annotation_column = "LAMC2_status_pos_neg",
    cell_types = c("Positive", "Negative"),
    contrast = c("LAMC2_status_pos_neg", "Positive", "Negative"),
    latent_vars = c("orig.ident"),
    lfcshrink = TRUE,
    out_dir = file.path(plot_dir, "lamc2_pos_vs_neg"))

print(head(res))
write.csv(res, file.path(plot_dir, "lamc2_pos_vs_neg.csv"))

volcano(res, out_dir = file.path(plot_dir, "lamc2_pos_vs_neg"))

### Using normalized data > 1 as threshold

threshold <- 1

# Annotate malignant cells as positive or negative
malignant$LAMC2_status <- case_when(
  lamc2_expression$LAMC2 > threshold ~ "Positive",
  lamc2_expression$LAMC2 <= threshold ~ "Negative",
  TRUE ~ "NA"
)

# Verify the annotation
table(malignant$LAMC2_status)

# DE
res <- perform_pseudobulk_DE(malignant, 
    sample_column = "sample.id", 
    annotation_column = "LAMC2_status",
    cell_types = c("Positive", "Negative"),
    contrast = c("LAMC2_status", "Positive", "Negative"),
    latent_vars = c("orig.ident"),
    lfcshrink = TRUE,
    out_dir = file.path(plot_dir, "lamc2_pos_vs_neg_threshold1"))

print(head(res))
write.csv(res, file.path(plot_dir, "lamc2_pos_vs_neg_threshold1.csv"))

volcano(res, out_dir = file.path(plot_dir, "lamc2_pos_vs_neg_threshold1"))

### Trying Q1 vs Q4

quartiles <- quantile(lamc2_expression$LAMC2, probs = c(0.25, 0.75))
print(quartiles)

# Annotate malignant cells as positive or negative
malignant$LAMC2_status <- case_when(
  lamc2_expression$LAMC2 >= quartiles[2] ~ "Positive",
  lamc2_expression$LAMC2 <= quartiles[1] ~ "Negative",
  TRUE ~ "NA"
)

# Verify the annotation
table(malignant$LAMC2_status)

# DE
res <- perform_pseudobulk_DE(malignant, 
    sample_column = "sample.id", 
    annotation_column = "LAMC2_status",
    cell_types = c("Positive", "Negative"),
    contrast = c("LAMC2_status", "Positive", "Negative"),
    latent_vars = c("orig.ident"),
    lfcshrink = TRUE,
    out_dir = file.path(plot_dir, "lamc2_pos_vs_neg_threshold_q1_vs_q4"))

print(head(res))
write.csv(res, file.path(plot_dir, "lamc2_pos_vs_neg_threshold_q1_vs_q4.csv"))

volcano(res, out_dir = file.path(plot_dir, "lamc2_pos_vs_neg_threshold_q1_vs_q4"))

### Trying Q1 vs Q10

quartiles <- quantile(lamc2_expression$LAMC2, probs = c(0.1, 0.9))

# Annotate malignant cells as positive or negative
malignant$LAMC2_status_q1_q10 <- case_when(
  lamc2_expression$LAMC2 >= quartiles[2] & malignant$Malignant == TRUE ~ "Positive",
  lamc2_expression$LAMC2 <= quartiles[1] & malignant$Malignant == TRUE ~ "Negative",
  TRUE ~ "NA"
)

# Verify the annotation
table(malignant$LAMC2_status_q1_q10)

# DE
res <- perform_pseudobulk_DE(malignant, 
    sample_column = "sample.id", 
    annotation_column = "LAMC2_status_q1_q10",
    cell_types = c("Positive", "Negative"),
    contrast = c("LAMC2_status_q1_q10", "Positive", "Negative"),
    latent_vars = c("orig.ident"),
    lfcshrink = TRUE,
    out_dir = file.path(plot_dir, "lamc2_pos_vs_neg_threshold_q1_vs_q10"))

print(head(res))
write.csv(res, file.path(plot_dir, "lamc2_pos_vs_neg_threshold_q1_vs_q10.csv"))

volcano(res, out_dir = file.path(plot_dir, "lamc2_pos_vs_neg_threshold_q1_vs_q10"))

### save object with lamc2 annotation

# pos-neg

pos <- WhichCells(malignant, expression = LAMC2_status_pos_neg == "Positive")
neg <- WhichCells(malignant, expression = LAMC2_status_pos_neg == "Negative")

so$celltype_lamc2_pos_neg <- so$cell.type.l2
so$celltype_lamc2_pos_neg[colnames(so) %in% pos] <- "malignant_LAMC2_positive"
so$celltype_lamc2_pos_neg[colnames(so) %in% neg] <- "malignant_LAMC2_negative"

print(table(so$celltype_lamc2_pos_neg))

# q1 q10

pos <- WhichCells(malignant, expression = LAMC2_status_q1_q10 == "Positive")
neg <- WhichCells(malignant, expression = LAMC2_status_q1_q10 == "Negative")

so$celltype_lamc2_q1_q10 <- so$cell.type.l2
so$celltype_lamc2_q1_q10[colnames(so) %in% pos] <- "malignant_LAMC2_positive"
so$celltype_lamc2_q1_q10[colnames(so) %in% neg] <- "malignant_LAMC2_negative"

print(table(so$celltype_lamc2_q1_q10))

saveRDS(so, "/cci_object.rds")