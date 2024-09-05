### pathway enrichment analysis

library(dplyr)
library(ggplot2)
library(patchwork)
library(fgsea)
library(data.table)

out_dir <- "path_to_output_directory"

plot_dir <- paste0(out_dir, "/pathway/")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---------------------------------------------------------------------------- #
#                             Load MsigDB pathways                             #
# ---------------------------------------------------------------------------- #

c5_go <- gmtPathways("/c5.go.v2023.2.Hs.symbols.gmt")
c2_cp <- gmtPathways("/c2.cp.v2023.2.Hs.symbols.gmt")
hallmark <- gmtPathways("/h.all.v2023.2.Hs.symbols.gmt")

msigdb <- list("c5_go" = c5_go, "c2_cp" = c2_cp, "hallmark" = hallmark)

# ---------------------------------------------------------------------------- #
#                                      ORA                                     #
# ---------------------------------------------------------------------------- #

de_res <- list.files("/DE", recursive = FALSE, full.names = TRUE)

de_res <- de_res[grepl(".csv", de_res)]

for (res in de_res){

    for (direction in c("up", "down")) {
        df <- read.csv(res)
        if(direction == "up") {
            curr_df <- df %>% filter(padj < 0.05 & log2FoldChange > 1)
        } else {
            curr_df <- df %>% filter(padj < 0.05 & log2FoldChange < -1)
        }

        universe_genes <- df$gene

        rownames(df) <- df$gene

        print(head(df))
        ranked_genes_response <- df %>%
            mutate(rank = log2FoldChange) %>%
            filter(rank != Inf) %>% filter(rank != -Inf) %>% arrange(-rank)

        ranked_genes_response$ranking <- rank(-ranked_genes_response$rank, ties.method = "first")

        response_level_stats <- ranked_genes_response$rank
        names(response_level_stats) <- rownames(ranked_genes_response)

        all_sig_results <- c()
        for(j in seq(1, length(msigdb))){
            curr_cat <- msigdb[[j]]

            ### GSEA

            # Run gsea and filter the pathways with a p-val less than 0.5
            fgseaRes_ctrl_multilevel <- fgseaMultilevel(pathways = curr_cat,
                                                        stats = response_level_stats,
                                                        minSize=10,
                                                        maxSize=2000,
                                                        nPermSimple = 10000) %>% 
                filter(padj < 0.5) %>% arrange(-abs(NES))


            gsea_res <- as.data.frame(fgseaRes_ctrl_multilevel)
            gsea_res <- gsea_res[!grepl("leadingEdge", colnames(fgseaRes_ctrl_multilevel))]
            write.csv(gsea_res, file.path(plot_dir, paste0("fgseaRes", "_", gsub(".csv", "", basename(res)), "_", names(msigdb)[j], ".csv")))

            ### ORA
            
            curr_ora_results <- fora(curr_cat, curr_df$gene, universe_genes)
            sig_ora_results <- curr_ora_results[padj < 0.05]
            
            if(nrow(sig_ora_results) > 10){
            setorder(sig_ora_results, padj)
            sig_ora_results <- sig_ora_results[1:10]
            }
            
            sig_ora_results$msigdb_cat <- paste0(names(msigdb)[j])
            all_sig_results <- rbindlist(list(all_sig_results, sig_ora_results))
        }
        
        if(nrow(all_sig_results) == 0){
            sig_result_plots <- NULL
        }else {
            if(length(unique(all_sig_results$msigdb_cat)) == 1){
            sig_result_plots <- ggplot(all_sig_results, aes(y = reorder(pathway, -padj), x = -log10(padj), color = -log10(padj))) + 
                geom_point(aes(size = size)) + 
                scale_color_gradient(low = "orange", high = "red") +
                labs(title = "Pathway Enrichment", y = "Pathway") +
                theme_classic()

            ggsave(filename = paste0(gsub(".csv", "", basename(res)), "_", direction, ".pdf"), plot =  sig_result_plots, path = plot_dir, width = 8, height = 15)
            }else {
            plots <- list()
            for (cat in unique(all_sig_results$msigdb_cat)) {

                ora_df <- all_sig_results %>% 
                filter(msigdb_cat == cat)

                plots[[cat]] <- ggplot(ora_df, aes(y = reorder(pathway, -padj), x = -log10(padj), color = -log10(padj))) +
                geom_point(aes(size = size)) +
                scale_color_gradient(low = "orange", high = "red") +
                labs(title = "Pathway Enrichment", y = "Pathway") +
                ggtitle(cat) +
                theme_classic()

            }
            plots <- wrap_plots(plots, ncol = 3, nrow = 3)
            ggsave(plots, filename = paste0("fora_indiv_", gsub(".csv", "", basename(res)), "_", direction, ".pdf"), path = plot_dir, height = 15, width = 25)

            all_sig_results <- all_sig_results %>% arrange(padj)
            }

            fwrite(all_sig_results, file = paste0(plot_dir, "/", gsub(".csv", "", basename(res)), "_", direction, ".csv"))
        }
    }

}