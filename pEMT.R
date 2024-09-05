### Analyzing Arora et al. ST data to investigate expansion of pEMT signature at LE compared to core

library(dplyr)
library(ggplot2)
library(Seurat)
library(data.table)
library(stringr)

## functions for lmm test

#' perform a nested t-test using LMMs.
#' 
#' @param data_df data frame containing data
#' @param response name of column in data_df that corresponds to the response variable
#' @param condition name of column in data_df that corresponds to condition of interest
#' @param latent_vars vector containing columns in data_df that correspond to variables that describe the nested structure of the data e.g. Patient
#' 
#' @return p-value derived from likelihood ratio test of model incorporating condition vs null model
#' @export
LMM_test <- function(data_df, response, condition, latent_vars) {
    
    require(lme4)
    
    full_formula <- as.formula(
        paste(response, "~", condition, "+", paste("(1 |", latent_vars, ")", collapse = "+"))
    )
    
    null_formula <- as.formula(
        paste(response, "~ 1 +", paste("(1 |", latent_vars, ")", collapse = "+"))
    )
    
    # Fit the models
    model <- lmer(full_formula, data = data_df)
    model_null <- lmer(null_formula, data = data_df)
    
    # Perform the likelihood ratio test
    res <- anova(model_null, model)
    
    # Return the p-value from the Chi-square test
    return(res$`Pr(>Chisq)`[2])
}


#' geom signif with lmm test. See geom_signif from ggsignif for more information
#' 
#' @param data_df data frame containing data
#' @param response name of column in data_df that corresponds to the response variable
#' @param condition name of column in data_df that corresponds to condition of interest
#' @param latent_vars vector containing columns in data_df that correspond to variables that describe the nested structure of the data e.g. Patient
#' @param comparisons list of vectors of length two that define the comparisons to make
#' @param y_position position of brackets 
#' @param tip_length length of bracket tips
#' @param size size of brackets
#' @param step_increase distance between brackets
#' 
#' @return geom signif
#' @export
geom_signif_lmm <- function(data_df, response, condition, latent_vars,
    comparisons, y_position = NULL, tip_length = 0.03, size = 0.5, step_increase = 0) {
    
    require(lme4)
    require(ggsignif)
    require(dplyr)

    p_vals <- sapply(seq_along(comparisons), function(comp) {
        
        curr_comparison <- comparisons[[comp]]
        curr_df <- data_df %>% 
            filter(!!sym(condition) %in% curr_comparison)
        
        return(round(LMM_test(curr_df, response, condition, latent_vars), digits = 5))
    })

    return(
        geom_signif(
            comparisons = comparisons,
            map_signif_level = FALSE,
            annotations = p_vals,
            y_position = y_position,
            tip_length = tip_length,
            size = size,
            step_increase = step_increase
        )
    )
}

out_dir <- "path_to_output_directory"

plot_dir <- paste0(out_dir, "/pEMT/")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# load data
files <- list.files("/Arora-Data", full.names = TRUE)

# signatures
sigs <- get(load("/puram.sigs.rDa"))

list_so <- list()
for (file in files) {
    so <- get(load(file))
    
    # subset for SCC spots i.e. tumor
    so <- subset(so, subset = pathologist_anno.x == "SCC")

    so <- AddModuleScore(so, features = list(sigs$p.EMT), name = "p_EMT", nbin = 15)

    list_so[[file]] <- so
}

so <- merge(list_so[[1]], list_so[-1])
so <- subset(so, subset = cluster_annotations != "nc")

so$cluster_annotations <- factor(so$cluster_annotations, levels = c("edge", "transitory", "core"))
p <- ggplot(so[[]], aes(x = cluster_annotations, y = p_EMT1)) +
    geom_boxplot(aes(fill = cluster_annotations)) +
    theme_classic() +
    scale_fill_manual(values = c("edge" = "#4c5760", "transitory" = "#93a8ac", "core" = "#d7ceb2")) +
    geom_signif_lmm(
        data_df = so[[]],
        response = "p_EMT1", 
        condition = "cluster_annotations",
        latent_vars = c("sample_id.x"),
        comparisons = list(c("edge", "transitory"), c("edge", "core")),
        step_increase = c(0, 0.1)
    )
ggsave(filename = "pEMT_edge_vs_other.pdf", path = plot_dir)

so$sample_clust <- paste0(so$sample_id.x, "_", so$cluster_annotations)
df <- AverageExpression(so, features = sigs$p.EMT, group.by = "sample_clust")
write.csv(df, file.path(plot_dir, "pEMT_avg_expression.csv"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

pal <- brewer.pal(12, "Set3")

df <- read.csv("/pEMT_avg_expression.csv", row.names = 1)
df <- log2(df + 1)
df <- df[, grepl("edge|core", colnames(df))]

fc_df <- sapply(seq(1,12), function(sample_num) {
    log_fc <- df[[paste0("SCT.sample_", sample_num, "_edge")]] - df[[paste0("SCT.sample_", sample_num, "_core")]]
    names(log_fc) <- rownames(df)
    return (log_fc)
}) 

fc_df <- as.data.frame(fc_df)
colnames(fc_df) <- paste0("Sample", seq(1,12))
fc_df$gene <- rownames(fc_df)
plot_df <- fc_df %>% 
    pivot_longer(-gene, names_to = "sample", values_to = "log_fc") %>% 
    group_by(gene) %>% 
    mutate(total = sum(log_fc)) %>% 
    arrange(total)
plot_df$gene <- factor(plot_df$gene, levels = unique(plot_df$gene))
plot_df$sample <- factor(plot_df$sample, levels = paste0("Sample", seq(1,12)))

p <- ggplot(plot_df, aes(x = gene, y = log_fc, fill = sample)) +
    geom_bar(position = "stack", stat = "identity", color = "black") +
    scale_fill_manual(values = pal) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 90)
    )
ggsave(filename = "pEMT_by_gene.pdf", path = plot_dir, width = 11)


