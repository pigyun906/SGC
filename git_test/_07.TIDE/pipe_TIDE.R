source("/home/jyhong906/Project/SGC/Script/RNA/git/_06.deconvolution/pipe_immunedeconv.R")

SGC_TIDE <- read.csv('/home/jyhong906/Project/SGC/Data/RNA/07.TIDE/TIDE_deseq_scale.csv',
                     header=T,
                     sep=",",
                     as.is=TRUE); rownames(SGC_TIDE) <- SGC_TIDE[,1]

SGC_TIDE <- as.data.frame(SGC_TIDE)[rownames(SGC_duct_groups),]
SGC_TIDE <- SGC_TIDE[,-c(1,2,3,4,6,10,14,15)]
correct_SGC_TIDE <- t(SGC_TIDE)

TIDE_ssGSEA <- t(scale(t(correct_SGC_TIDE)))
cols <- colorRamp2(c(min(TIDE_ssGSEA), mean(TIDE_ssGSEA), max(TIDE_ssGSEA)),
                   c('#392BC0','white','#C0392B'))
top_annotation <- HeatmapAnnotation('Condition' = SGC_groups$Condition,
                                    col          = list("Condition" = c("ACC" = "#FEC260", "MECA" = "#3FA796", "SDC" = "#2A0944", "MEC" = "#A10035"),
                                                        gp = gpar(col = "grey")))

TIDE_ssGSEA <- TIDE_ssGSEA[,rownames(SGC_groups)]

SGC_groups$Condition <- factor(SGC_groups$Condition, levels = c("ACC","MECA","SDC","MEC"))

# save in /home/jyhong906/Project/SGC/Data/RNA/07.TIDE/Fig 3A. TIDE (4 x 7)
Immune_heatmap <- Heatmap(TIDE_ssGSEA,
                          # column_title           = "Immune ssGSEA", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                          column_split           = SGC_groups$Condition,
                          col                    = cols,
                          clustering_method_rows = "complete",
                          show_column_names      = T,
                          show_row_names         = T,
                          show_row_dend          = T,
                          show_column_dend       = T, 
                          name                   = "Z-score",
                          row_names_gp           = gpar(fontsize = 8),
                          cluster_columns        = T,
                          cluster_rows           = T,
                          top_annotation         = top_annotation); Immune_heatmap