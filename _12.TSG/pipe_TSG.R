source("/home/jyhong906/Project/SGC/Script/RNA/git/_11.TCR/pipe_MIXCR.R")

TSG_geneset <- as.data.frame(read.csv('/home/jyhong906/Public/TSG/TSGene set.csv',
                                         header=F,
                                         sep=",",
                                         as.is=T))
rownames(TSG_geneset) <- TSG_geneset[,1]
TSG_geneset <- TSG_geneset[,-1]
TSG_geneset <- data.frame(t(TSG_geneset))
TSG_geneset <- lapply(TSG_geneset, function(x) x[x != "" ])
for (i in 1:length(TSG_geneset)){
          TSG_geneset[[i]] <- as.character(TSG_geneset[[i]])
}

TSG_ssGSEA <- gsva(SGC_mat,
                      TSG_geneset,
                      min.sz=1,
                      method="ssgsea",
                      ssgsea.norm=F, 
                      max.sz=999999,
                      abs.ranking=F,
                      verbose=T)

TSG_score_df <- as.data.frame(matrix(nrow = 94 * 1, ncol = 3))
colnames(TSG_score_df) <- c("group", "weight", "TSG")
TSG_score_df$weight <- scale(TSG_ssGSEA[1,])
TSG_score_df$group <- c(rep("ID", 58), rep("ED", 36))
TSG_score_df$facet <- rep("TSG", 94)
TSG_score_df$group <- factor(TSG_score_df$group, levels = c("ID", "ED"))

# Save in /home/jyhong906/Project/SGC/Data/RNA/14.TSG/Fig. TSG score.pdf (5 x 4.5)
TSG_p <- ggplot(TSG_score_df,
               aes(x = group,
                   y = weight,
                   fill = group)) +
          
          # dot
          geom_boxplot(width = 0.5,
                       alpha = 0.5) +
          
          # Add an color pallate
          scale_fill_manual(values = c("ID" = "#800080", # red
                                       "ED" = "#FF884B")) +
          
          geom_dotplot(method   = "dotdensity",
                       alpha    = 1,
                       binwidth = (max(TSG_score_df$weight) - min(TSG_score_df$weight))/50,
                       colour   = "black",
                       binaxis  = 'y',
                       stackdir = 'center',
                       position = position_dodge(0.5)) +
          
          # significance
          stat_compare_means(comparisons = list(c("ID","ED")),
                             label = "p.signif",
                             size = 6,
                             position = position_nudge(y = 1),
                             method = "wilcox.test") + # Add significance levels
          # stat_compare_means(position = position_nudge(y = 0.2),
          #                    size = 4,
          #                    method = "wilcox.test") +
          
          # white backgroud
          theme_bw() +
          
          #eliminates background, gridlines, and chart border
          theme(axis.title       = element_text(size = 20,
                                                face = 'plain'),
                title            = element_text(size = 25,
                                                face = 'bold'),
                plot.title       = element_text(hjust = 0.5),
                plot.background  = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border     = element_rect(colour = "black",
                                                fill   = NA,
                                                size   = 1,),
                legend.position = "bottom",
                legend.title = element_blank()) +
          
          # Draws x and y axis line
          theme(axis.line = element_line(color = 'black'),
                axis.text.x = element_text(size = 10,
                                           color = "black"),
                axis.text.y = element_text(size = 20,
                                           color = "black")); TSG_p

# CSC-TSG correlation #
# Save in /home/jyhong906/Project/SGC/Data/RNA/14.TSG/Fig. Correlation between TSG and stemness.pdf (6.5 x 7)
CSC_TSG_df <- data.frame("CSC" = CS_df$Stemness_Index,
                         "TSG" = TSG_score_df$weight,
                         "origin" = SGC_duct_groups$Duct)
CSC_TSG_p <- ggscatterhist(CSC_TSG_df,
                           x = "CSC",
                           y = "TSG",
                           color = "origin",
                           palette = c("ID" = "#800080",
                                       "ED" = "#FFA500"),
                           size = 3, alpha = 1,
                           xlab = "Cancer stemness",
                           ylab = "Tumor suppressor gene score",
                           legend = "none",
                           add = "reg.line",
                           margin.plot.size = 0,
                           margin.params = list(fill = "origin", color = "black", size = 0.1, alpha = 0.7)); CSC_TSG_p
CSC_TSG_p$sp +

          stat_cor(label.x = 0.5,
                   label.y = max(CSC_TSG_df$TSG),
                   size = 7) +

          # white background
          theme_classic() +

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
                    axis.text.x = element_text(size = 20,
                                               colour = "black",
                                               margin = margin(5,0,0,0)),
                    axis.text.y = element_text(size = 20,
                                               colour = "black",
                                               hjust = -0.1,
                                               margin = margin(0,5,0,0)),
                    axis.title.x = element_text(size = 20,
                                                face = "bold",
                                                vjust = -1,
                                                margin = margin(10,0,10,0)),
                    axis.title.y = element_text(size = 20,
                                                face = "bold",
                                                vjust= 4,
                                                margin = margin(0,10,0,10)))
