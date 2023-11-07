source("/home/jyhong906/Project/SGC/Script/RNA/git/_10.csc/pipe_pancanstem.R")

merge_immune_df <- rbind(as.data.frame(Immune_ssGSEA[c(1:24),]), new_immune_scaled_mat[c("T cell NK","T cell CD8+ naive"),], TIDE_ssGSEA[c(7,8),])
CSC_cor_df <- data.frame("cell" = seq(nrow(merge_immune_df)),
                         "cor" = seq(nrow(merge_immune_df)),
                         "pval" = seq(nrow(merge_immune_df))); rownames(CSC_cor_df) <- rownames(merge_immune_df); CSC_cor_df$cell <- rownames(merge_immune_df)

for (idx in seq(nrow(merge_immune_df))) {
          cor_test <- cor.test(CS_df$Stemness_Index, as.numeric(merge_immune_df[idx,]))
          CSC_cor_df$cor[idx] <- cor_test$estimate; CSC_cor_df$pval[idx] <- cor_test$p.value             
}; CSC_cor_df$padj <- p.adjust(CSC_cor_df$pval, method = "fdr")

CSC_cor_df$cor <- as.numeric(CSC_cor_df$cor); CSC_cor_df$color <- "gray"
CSC_cor_df[CSC_cor_df$padj <= 0.05 & CSC_cor_df$cor > 0,]$color <- "red"
CSC_cor_df[CSC_cor_df$padj <= 0.05 & CSC_cor_df$cor < 0,]$color <- "blue"

CSC_cor_df$cell <- factor(CSC_cor_df$cell, levels = CSC_cor_df[order(CSC_cor_df$cor, decreasing = T),]$cell)

# Save in /home/jyhong906/Project/SGC/Data/RNA/12.Cancer stem cell/Fig. CSC-immune.cells.correlation.pdf (5 x 8)
CSC_cor_df %>% ggplot(aes(cell, cor, fill = color)) +
          
          geom_col() +
          scale_fill_manual(values = c("gray" = "gray",
                                       "red" = "#C0392B",
                                       "blue" = "#392BC0")) +
          
          # coord_flip() +

          # white background
          theme_classic() +
          
          # axis, main title
          labs(x     = "Immune cells",
               y     = "Correlation with CSCs") +
          
          # eliminates background, gridlines, and chart border
          theme(
                    # Plot background
                    plot.background  = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    
                    # Legend 
                    legend.position = "none",
                    # legend.title = element_blank(),
                    # legend.text = element_text(size = 15),
                    
                    
                    axis.line = element_line(color = 'black'),
                    axis.text.x = element_text(size = 12,
                                               colour = "black",
                                               angle = 90),
                    axis.text.y = element_text(size = 10,
                                               colour = "black"),
                    axis.title.x = element_text(size = 18,
                                                face = "bold",
                                                vjust = -1,
                                                margin = margin(0,0,10,0)),
                    axis.title.y = element_text(size = 18,
                                                face = "bold",
                                                vjust= 3,
                                                margin = margin(0,5,0,15)))

                    