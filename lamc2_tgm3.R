### Check distribution of LAMC2 and TGM3, divide cells as high/low

library(dplyr)
library(ggplot2)
library(Seurat)
library(stringr)
library(ggpubr)

out_dir <- "path_to_output_directory"

plot_dir <- paste0(out_dir, "/lamc2_tgm3/")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

### Load object
hnscc <- readRDS("/hnscc.rds")

# Remove HPV+ samples
hnscc <- subset(hnscc, subset = hpv == "HPV-")

### What is global expression pattern of LAMC2 and TGM3?

p <- VlnPlot(hnscc, features = c("LAMC2", "TGM3"), group.by = "cell.type")
ggsave(filename = "global_LAMC2_TGM3_expr.pdf", path = plot_dir)

p <- DotPlot(hnscc, features = c("LAMC2", "TGM3"), group.by = "cell.type")
ggsave(filename = "dotplot_global_LAMC2_TGM3_expr.pdf", path = plot_dir)

### LAMC2 and TGM3 expression in malignant vs non-malignant epithelial cells
epi <- subset(hnscc, subset = cell.type == "Epithelial")

p <- VlnPlot(epi, features = c("LAMC2", "TGM3"), group.by = "Malignant")
ggsave(filename = "epi_LAMC2_TGM3_expr.pdf", path = plot_dir)

### What is the distribution of LAMC2 expression within malignant cells?

malignant_epi <- subset(epi, subset = Malignant == TRUE)

LAMC2_expr <- FetchData(object = malignant_epi, vars = "LAMC2")

p <- ggplot(LAMC2_expr, aes(x = LAMC2)) +
  geom_histogram() +
  xlab("Expression Level") +
  ylab("Frequency") +
  ggtitle("LAMC2 Expression Histogram")
ggsave(filename = "LAMC2_expr_histo.pdf", path = plot_dir)

### What is the distribution of TGM3 expression within malignant cells?

TGM3_expr <- FetchData(object = malignant_epi, vars = "TGM3")

ggplot(TGM3_expr, aes(x = TGM3)) +
  geom_histogram() +
  xlab("Expression Level") +
  ylab("Frequency") +
  ggtitle("TGM3 Expression Histogram")
ggsave(filename = "TGM3_expr_histo.pdf", path = plot_dir)  

# ---------------------------------------------------------------------------- #
#                   Compare LAMC2 expression malignant vs epi                  #
# ---------------------------------------------------------------------------- #

epi <- subset(hnscc, subset = cell.type == "Epithelial")
epi$sample_malignant <- paste0(epi$sample.id, "_", epi$Malignant)
avg_expr <- AverageExpression(epi, group.by = "sample_malignant", features = c("LAMC2"))
write.csv(avg_expr, file.path(plot_dir, "mal_vs_epi_lamc2_expr.csv"))
avg_expr <- read.csv("/mal_vs_epi_lamc2_expr.csv", row.names = 1)

df <- as.data.frame(t(avg_expr))
colnames(df) <- c("avg_lamc2_expression")
group <- sapply(str_split(rownames(df), "_"), function(x) {
  if ("CD45n" %in% x) {
    return(x[3])
  } else {
    return(x[2])
  }
  
})
df$group <- ifelse(group == "TRUE", "Malignant", "Non-malignant")
print(head(df))

p <- ggplot(df, aes(x=group, y=avg_lamc2_expression)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  stat_compare_means(comparisons = list(c("Malignant", "Non-malignant")))
ggsave(filename = "mal_vs_epi_lamc2_expr.pdf", path = plot_dir)
