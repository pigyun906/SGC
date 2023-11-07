# install.packages("estimate", repos="http://R-Forge.R-project.org")

source("/home/jyhong906/Project/SGC/Script/RNA/git/_03.cancer hallmark/pipe_pca_hallmark.R")
Immune_geneset <- as.data.frame(read.csv('/home/jyhong906/Public/TME/TIS, IIS score gene set.csv',
                                         header=F,
                                         sep=",",
                                         as.is=T))
rownames(Immune_geneset) <- Immune_geneset[,1]
Immune_geneset <- Immune_geneset[,-1]
Immune_geneset <- data.frame(t(Immune_geneset))
Immune_geneset <- lapply(Immune_geneset, function(x) x[x != "" ])
for (i in 1:length(Immune_geneset)){
          Immune_geneset[[i]] <- as.character(Immune_geneset[[i]])
}

Immune_ssGSEA <- gsva(SGC_mat,
                      Immune_geneset,
                      min.sz=1,
                      method="ssgsea",
                      ssgsea.norm=F, 
                      max.sz=999999,
                      abs.ranking=F,
                      verbose=T)

IIS_score <- colMeans(Immune_ssGSEA[c(1:18, 21:24),]) # IIS_score, no scaled
TIS_score <- colMeans(Immune_ssGSEA[c(3,15,16,17,18,21,22,23,24),]) # TIS_score, no scaled
CYT_score <- sqrt(SGC_mat[rownames(SGC_mat)=='GZMA',] * SGC_mat[rownames(SGC_mat)=='PRF1',]) # CYT_Score, no scaled
ANG_score <- Immune_ssGSEA[26,]
APM_score <- Immune_ssGSEA[27,] 

Custom_immune_score_df <- as.data.frame(matrix(nrow = 94 * 5, ncol = 3))
colnames(Custom_immune_score_df) <- c("group", "weight", "Immune")
Custom_immune_score_df$weight <- c(scale(IIS_score), scale(TIS_score), scale(CYT_score), scale(ANG_score), scale(APM_score))
Custom_immune_score_df$group <- rep(c(rep("ACC", 18), rep("MECA", 40), rep("SDC", 16), rep("MEC", 20)), 5)
Custom_immune_score_df$Immune <- c(rep("IIS", 94), rep("TIS", 94), rep("CYT", 94), rep("ANG", 94), rep("APM", 94))
Custom_immune_score_df$group <- factor(Custom_immune_score_df$group, levels = c("ACC", "MECA", "SDC", "MEC"))
Custom_immune_score_df$Immune <- factor(Custom_immune_score_df$Immune, levels = c("IIS", "TIS", "CYT", "APM","ANG"))
ggplot(Custom_immune_score_df,
       aes(x     = group,
           y     = weight,
           fill  = group)) +
          
          # dot
          geom_boxplot(width = 0.5,
                       alpha = 0.5) +
          
          geom_dotplot(method   = "dotdensity",
                       alpha    = 1,
                       binwidth = (max(Custom_immune_score_df$weight) - min(Custom_immune_score_df$weight))/50,
                       colour   = "black",
                       binaxis  = 'y',
                       stackdir = 'center') +
          
          # color
          scale_color_manual(values=SGC_cols) + # line
          scale_fill_manual(values=SGC_cols) + # fill
          
          # significance
          stat_compare_means(comparisons = list(c("ACC","MECA"),c("MECA","SDC"),c("SDC","MEC")),
                             label = "p.signif",
                             size = 4) + # Add significance levels
          stat_compare_means(label.x = 1,
                             # label.y = max(Custom_immune_score_df$weight) + 2,
                             position = position_nudge(y = 2),
                             size = 4) +
          
          # axis, main title
          labs(x     = "Subtypes",
               y     = "Normalized expression",
               title = "TIME score") +
          
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
                                           color = "black")) +
          
          # Make a multiple boxes
          facet_wrap( ~ Immune, scales="free")