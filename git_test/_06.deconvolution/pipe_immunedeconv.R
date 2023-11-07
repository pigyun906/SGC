source("/home/jyhong906/Project/SGC/Script/RNA/git/_05.TCGA/pipe_pan-cancer.R")

gene_cov <- get(load("/home/jyhong906/Public/Gene_cov/gene_cov.rda"))
normalized_TPM <- countToTpm(Keep_SGC_matrix, # raw counts #
                             keyType = "SYMBOL",
                             gene_cov = gene_cov)

deconvolution_methods
immunedeconv::timer_available_cancers

tool_df_list <- c()
for (idx in seq(length(deconvolution_methods))){
          tool <- deconvolution_methods[idx]
          
          if (tool %in% c('cibersort', "cibersort_abs")) {
                    set_cibersort_binary('/home/jyhong906/Project/CRC/script/RNA/05.Immune cell composition/Immunedeconv/CIBERSORT/CIBERSORT.R')
                    set_cibersort_mat('/home/jyhong906/Project/CRC/script/RNA/05.Immune cell composition/Immunedeconv/CIBERSORT/LM22.txt')
                    print(tool)
                    assign(tool, as.data.frame(deconvolute(normalized_TPM, tool)))
                    assign(paste0(tool, '_df'), data.frame(get(tool),
                                                           row.names = get(tool)$cell_type)[,-1])
                    
          }
          
          else if (tool %in% 'timer') {
                    print(tool)
                    assign(tool, deconvolute(normalized_TPM, tool,
                                             indications=c(rep("hnsc", ncol(normalized_TPM)))))
                    assign(paste0(tool, '_df'), data.frame(get(tool),
                                                           row.names = get(tool)$cell_type)[,-1])
                    
          }
          
          else {
                    print(tool)
                    assign(tool, as.data.frame(deconvolute(normalized_TPM, tool)))
                    assign(paste0(tool, '_df'), data.frame(get(tool),
                                                           row.names = get(tool)$cell_type)[,-1])
          }
          
          tool_df_list <- append(tool_df_list, paste0(tool, '_df'))
}

total_immune_df <- rbind(mcp_counter_df, epic_df, quantiseq_df, xcell_df, cibersort_abs_df, timer_df)
immune_scaled_matrix <- t(scale(t(total_immune_df)))

immune_scaled_matrix[!is.na(immune_scaled_matrix[,1]),]
method_df <- as.data.frame(matrix(nrow = nrow(total_immune_df), ncol = 2))
colnames(method_df) <- c('type', 'method')
method_df$type <- rownames(total_immune_df)
method_df$method <- c(rep('MCP_counter', nrow(mcp_counter_df)),
                      rep('Epic', nrow(epic_df)),
                      rep('quantiseq', nrow(quantiseq_df)),
                      rep('xCell', nrow(xcell_df)),
                      rep('CIBERSORT', nrow(cibersort_abs_df)),
                      rep('TIMER', nrow(timer_df)))
method_df <- method_df[!is.na(immune_scaled_matrix[,1]),]

new_immune_scaled_mat <- immune_scaled_matrix[method_df$type,]
cols <- colorRamp2(c(min(new_immune_scaled_mat[!is.na(new_immune_scaled_mat[,1]),]), mean(new_immune_scaled_mat[!is.na(new_immune_scaled_mat[,1]),]), max(new_immune_scaled_mat[!is.na(new_immune_scaled_mat[,1]),])),
                   c('#392BC0','white','#C0392B'))
top_annotation <- HeatmapAnnotation('Condition' = SGC_groups$Condition,
                                    col          = list("Condition" = c("ACC" = "#FEC260", "MECA" = "#3FA796", "SDC" = "#2A0944", "MEC" = "#A10035"),
                                                        gp = gpar(col = "grey")))

left_annotation <- rowAnnotation(Method = method_df$method,
                                 # P_value = method_df[order(method_df$method, method_df$p_value, decreasing = T),]$p_value,
                                 # Immune_cell_type = method_df$immune_cell_type,
                                 col = list("Method" = c("xCell" = "#DB9D85", "TIMER" = "#BAAC65", "quantiseq" = "#86B875", "MCP_counter" = "#47BEA2", "Epic" = "#4CB9CC", "CIBERSORT" = "#96AAE1")),
                                 width = unit(1, "cm"))

# save in Fig. 1K. Immune cell composition.pdf (15 x 12)
Immunedeconv_heatmap <- Heatmap(new_immune_scaled_mat,
                                # column_title           = "Immune ssGSEA", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                                # column_split           = SGC_groups$Condition,
                                row_split = as.factor(method_df$method),
                                col                    = cols,
                                clustering_distance_columns = "pearson",
                                clustering_method_rows = "complete",
                                show_column_names      = F,
                                show_row_names         = T,
                                show_row_dend          = T,
                                show_column_dend       = T, 
                                name                   = "Z-score",
                                row_names_gp           = gpar(fontsize = 8),
                                cluster_columns        = T,
                                cluster_rows           = F,
                                # column_km = 2,
                                # row_km = 2,
                                left_annotation        = left_annotation,
                                top_annotation         = top_annotation); Immunedeconv_heatmap
