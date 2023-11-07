source("/home/jyhong906/Project/SGC/Script/RNA/git/_01.count/pipe_htseq.R")
SGC_cols <- c("ACC"="#FFB200","MECA"="#3D8361","SDC"="#25316D","MEC"="#D2001A")

# Make DESeqDataSet & running DESeq
SGC_dds <- DESeqDataSetFromMatrix(countData = Keep_SGC_matrix, # 19055 genes
                              colData = SGC_groups,
                              design = ~ Condition)

          # Un-normalized expression boxplot #
          par(mar=c(8,5,2,2))
          Un_normalized_boxplot <- boxplot(log10(Keep_SGC_matrix + 1),
                                           range=0,
                                           las=2,
                                           main='Unnormalized boxplot')

SGC_dds <- DESeq(SGC_dds)

SGC_dat <- counts(SGC_dds, normalized=T)
mod <- model.matrix(~ Condition, colData(SGC_dds))
mod0 <- model.matrix(~ 1, colData(SGC_dds))
SGC_svseq <- svaseq(SGC_dat, mod, mod0, n.sv=2)

SGC_svadds <- SGC_dds
SGC_svadds$SV1 <- SGC_svseq$sv[,1]
SGC_svadds$SV2 <- SGC_svseq$sv[,2]
design(SGC_svadds) <- ~ SV1 + SV2 + Condition
SGC_svadds <- DESeq(SGC_svadds)

SGC_vsd <- vst(SGC_svadds)
SGC_mat <- assay(SGC_vsd); print("SGC_mat")

          par(mar=c(8,5,2,2))
          boxplot(vst(assay(SGC_dds)),
                  range=0,
                  las=2,
                  main='Normalized boxplot')
          
# write.table(SGC_mat, "/home/jyhong906/Project/SGC/Data/RNA/02.DESeq/matrix.txt",
#             sep="\t",
#             row.names = T,
#             col.names = T,
#             quote=F,
#             append=F,
#             na="NA")
SGC_norm_mat <- scale(SGC_mat)

# PCA plot after ComBat-seq #
var_genes <- apply(SGC_norm_mat, 1, var)
select_var <- names(sort(var_genes, decreasing=T))[1:3000]
SGC_top_cor_mat <- SGC_norm_mat[select_var,]

# Fig. 1A. PCA plot #
# Cor_obj <- prcomp(SGC_mat, rank = 10)
Cor_obj <- prcomp(SGC_top_cor_mat, rank = 10) # 논문용
Cor <- as.data.frame(Cor_obj[2]$rotation)
Cor[,"Condition"] = SGC_groups$Condition
Cor[,"Batch"] = SGC_groups$Batch
Cor_percent <- round(Cor_obj$sdev / sum(Cor_obj$sdev) * 100, 2)
Cor_percent <- paste0(colnames(Cor)[which(substr(colnames(Cor),1,2) == "PC")],"(",paste0(as.character(Cor_percent),"%",")", sep=""))
Cor$label <- rownames(Cor)

# save in 7 x 6.8
Cor_plot <- ggplot(data = Cor,
                   aes(x = PC1,
                       y = PC2,
                       color = Condition,
                       label = rownames(Cor))) + # , label = Batch
          geom_point(size=3) +
          # geom_text(size=3) + 
          stat_ellipse() +
          labs(title="PCA plot - Total genes") +
          scale_colour_manual(values = SGC_cols) +
          theme_bw() +
          
          #eliminates background, gridlines, and chart border
          theme(axis.title       = element_text(size = 15,
                                                face = 'plain'),
                title            = element_text(size = 0,
                                                face = 'bold'),
                plot.title       = element_text(hjust = 0.5),
                plot.background  = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border     = element_rect(colour = "black",
                                                fill   = NA,
                                                size   = 1,),
                legend.position = "bottom",
                legend.title = element_blank(),
                legend.text = element_text(size = 10)) +
          
          #draws x and y axis line
          theme(axis.line = element_line(color = 'black'),
                axis.text.x = element_text(size = 10),
                axis.text.y = element_text(size = 10)) +
          
          xlab(Cor_percent[1]) +
          ylab(Cor_percent[2]); Cor_plot

# Cluster dendrogram #
Fig.1B_df <- t(SGC_top_cor_mat)
Fig.1B_dist <- dist(x      = Fig.1B_df,
                   method = "euclidean")
Fig.1B_cls <- hclust(d      = Fig.1B_dist,
                     method = "complete")

# Customized colors
Fig.1B_dend <- as.dendrogram(Fig.1B_cls, hang = -1) %>% set("branches_k_color", 
                                   value = c("#674747", "#829460"), k = 2) %>% set("branches_lwd", 1.2) %>%
          set("labels", "") %>% set("branches_lwd", 1) %>% set("leaves_pch", 19) %>% set("leaves_cex", 2) %>% set("leaves_col", SGC_cols[SGC_groups[Fig.1B_cls$order,]$Condition])

ggplot(as.ggdend(Fig.1B_dend)) +
          theme_bw()

top_annotation <- HeatmapAnnotation('Condition' = SGC_groups[Fig.1B_cls$order,]$Condition,
                                    col          = list("Condition" = SGC_cols[SGC_groups$Condition],
                                                        gp = gpar(col = "grey")))