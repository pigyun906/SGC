source("/home/jyhong906/Project/SGC/Script/RNA/git/_02.exp/pipe_deseq2.R")

Hallmark_order <- read.table("/home/jyhong906/Public/GSEA_gmt/MSigDB hallmark order.txt",
                             sep = '\t',
                             row.names = 1); Hallmark_order$V2 <- str_to_title(Hallmark_order$V2); Hallmark_order$V2 <- gsub("Dna Damage", "DNA damage", Hallmark_order$V2)
Hallmark_order$V2 <- gsub("Cellular Component", "Cellular component", Hallmark_order$V2)
Hallmark_list <- gmtPathways("/home/jyhong906/Public/GSEA_gmt/h.all.v7.2.symbols.gmt")
Hallmark_score <- gsva(Cor_obj$x[,c(1:5)],
                       Hallmark_list,
                       min.sz=1,
                       method="ssgsea",
                       ssgsea.norm=F, 
                       max.sz=999999,
                       abs.ranking=F,
                       verbose=T); rownames(Hallmark_score) <- sapply(str_split(rownames(Hallmark_score), "HALLMARK_"), function(x) {x[2]})

Hallmark_scaled_matrix <- t(scale(t(abs(Hallmark_score))))
cols <- colorRamp2(c(0, max(Hallmark_scaled_matrix)),
                   c("#FFFFFF", "#FFB72B"))

left_annotation <- rowAnnotation(Category = Hallmark_order$V2,
                                 col = list("Category" = c("Cellular component" = "#d31e25",
                                                           "Development" = "#d7a32e",
                                                           "DNA damage" = "#d1c02b",
                                                           "Immune" = "#369e4b",
                                                           "Metabolic" = "#5db5b7",
                                                           "Pathway" = "#31407b",
                                                           "Proliferation" = "#8a3f64",
                                                           "Signaling" = "#4f2e39")),
                                 width = unit(1, "cm"))

# save in .pdf (10 x 5.1)
Correct_mat <- Hallmark_scaled_matrix[rownames(Hallmark_order),]
rownames(Correct_mat) <- gsub("_", " ", rownames(Correct_mat))
pca_heatmap <- Heatmap(Correct_mat,
                       col = cols,
                       clustering_method_rows = "complete",
                       show_column_names      = T,
                       show_row_names         = T,
                       show_row_dend          = T,
                       show_column_dend       = T, 
                       name                   = "GSVA score",
                       row_names_gp           = gpar(fontsize = 8),
                       # row_km = 5,
                       cluster_columns        = F,
                       cluster_rows           = F,
                       left_annotation        = left_annotation); pca_heatmap