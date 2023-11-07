source("/home/jyhong906/Project/SGC/Script/RNA/git/_12.TSG/pipe_TSG.R")

# TCR / TCR #
# Fig. 3B

# MHC_class_I #
MHC_class_I_df <- as.data.frame(scale(Immune_ssGSEA[rownames(Immune_ssGSEA) == "MCH.class.I",]))
MHC_I <- scale(MHC_class_I_df$V1)

# MHC_II / MHC_II #
MHC_class_II_df <- as.data.frame(scale(Immune_ssGSEA[rownames(Immune_ssGSEA) == "MCH.class.II",]))
MHC_II <- scale(MHC_class_II_df$V1)

Signal_list <- c("MHC_I", "MHC_II", "CD28","CD40","CD40L","CD80","CD86","ICOS","ICOSLG","OX40","OX40L","CD80","CTLA4","GAL9","PD-1","PDL1","PDL2","TIM3")
col_list <- c('variable','subtype','value',"facet")
signal_df <- as.data.frame(matrix(nrow = ncol(Immune_ssGSEA) * length(Signal_list), ncol = length(col_list)))
colnames(signal_df) <- col_list
signal_df$subtype <- rep(SGC_duct_groups$Duct, length(Signal_list))
signal_df$facet <- c(rep("Signal1", ncol(Immune_ssGSEA) * 2),
                     rep("Signal2c", ncol(Immune_ssGSEA) * 9),
                     rep("Signal2i", ncol(Immune_ssGSEA) * 7))

# CD80 / CD80 #
CD80 <- scale(SGC_mat[grep('CD80',rownames(SGC_mat)),])

# CD86 / CD86 #
CD86 <- scale(SGC_mat[grep('CD86',rownames(SGC_mat)),])

# CD40 / CD40 #
CD40 <- scale(SGC_mat[grep('CD40',rownames(SGC_mat)),][1,])

# OX40L / TNFSF4 #
OX40L <- scale(SGC_mat[grep('TNFSF4',rownames(SGC_mat)),])

# CD28 / CD28 #
CD28 <- scale(SGC_mat[grep('CD28',rownames(SGC_mat)),])

# CD40L / CD40LG #
CD40L <- scale(SGC_mat[grep('CD40LG',rownames(SGC_mat)),])

# OX40 / TNFRSF4 #
OX40 <- scale(SGC_mat[grep('TNFRSF4',rownames(SGC_mat)),])

# ICOS / ICOS #
ICOS <- scale(SGC_mat[grep('ICOS',rownames(SGC_mat)),][1,])

# ICOSLG / ICOSG #
ICOSLG <- scale(SGC_mat[grep('ICOSLG',rownames(SGC_mat)),])

# PD-1 / PDCD1 #
PD1 <- scale(SGC_mat[grep('PDCD1',rownames(SGC_mat)),][1,])

# PD-L1 / CD274 #
PDL1 <- scale(SGC_mat[grep('CD274',rownames(SGC_mat)),])

# PD-L2 / PDCD1LG2 #
PDL2 <- scale(SGC_mat[grep('PDCD1LG2',rownames(SGC_mat)),])

# CD80 / CD80 #
CD80 <- scale(SGC_mat[grep('CD80',rownames(SGC_mat)),])

# GAL9 / LGALS9 #
GAL9 <- scale(SGC_mat[grep('LGALS9',rownames(SGC_mat)),])

# TIM3 / HAVCR2 #
TIM3 <- scale(SGC_mat[grep('HAVCR2',rownames(SGC_mat)),])

# CTLA4 / CTLA4 #
CTLA4 <- scale(SGC_mat[grep('CTLA4',rownames(SGC_mat)),])

signal_df$value <-  c(MHC_I, MHC_II, CD28, CD40, CD40L, CD80, CD86, ICOS, ICOSLG, OX40, OX40L, CD80, CTLA4, GAL9, PD1, PDL1, PDL2, TIM3)
signal_df$variable <- c(rep('MHC I',ncol(Immune_ssGSEA)),rep('MHC II',ncol(Immune_ssGSEA)),
                        rep('CD28',ncol(Immune_ssGSEA)),rep('CD40',ncol(Immune_ssGSEA)),rep('CD40L',ncol(Immune_ssGSEA)),
                        rep('CD80',ncol(Immune_ssGSEA)),rep('CD86',ncol(Immune_ssGSEA)),rep('ICOS',ncol(Immune_ssGSEA)),
                        rep('ICOSLG',ncol(Immune_ssGSEA)),rep('OX40',ncol(Immune_ssGSEA)),rep('OX40L',ncol(Immune_ssGSEA)),
                        rep('CD80',ncol(Immune_ssGSEA)),rep('CTLA4',ncol(Immune_ssGSEA)),rep('GAL9',ncol(Immune_ssGSEA)),
                         rep('PD1',ncol(Immune_ssGSEA)),rep('PDL1',ncol(Immune_ssGSEA)),rep('PDL2',ncol(Immune_ssGSEA)),
                         rep('TIM3',ncol(Immune_ssGSEA)))

signal_df$subtype <- factor(signal_df$subtype, levels = c("ID", "ED"))

# save in /home/jyhong906/Project/SGC/Data/RNA/15.Signaling molecules/Immune signaling molecules.pdf (6 x 25)
signal_p <- ggplot(signal_df,
                  aes(x=variable, y=value, fill=subtype)) +
          geom_boxplot(width = .2, alpha = .5, fatten = NULL, show.legend = F, outlier.alpha = 0) +
          introdataviz::geom_split_violin(alpha = .4, trim = F) +
          stat_summary(fun.y = median,
                       geom = "crossbar",
                       position = position_dodge(width = 0.2),
                       width = 0.2) +
          scale_fill_manual(values = Duct_cols) + # fill
          
          stat_compare_means(
                    label = "p.signif",
                    # symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                    #                    symbols = c("****", "***", "**", "*", "ns")),
                    method = "wilcox.test",
                    position = position_nudge(y = 1.5),
                    size = 5) +
          
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
          facet_wrap( ~ facet, scales="free"); signal_p
