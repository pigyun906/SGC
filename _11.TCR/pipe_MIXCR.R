source("/home/jyhong906/Project/SGC/Script/RNA/git/_10.csc/pipe_cscwithimmune.R")

file_path = "/home/jyhong906/Project/SGC/Data/RNA/13.TCR"
immdata_mixcr <- repLoad(file_path) # ACC 4 samples, MECA 2 samples TCR none.
a <- immdata_mixcr$meta$Sample

# ID <- a[c(1:15,36:45)]
# ED <- a[c(46:61,16:35)]
# immdata_mixcr$meta[1:25,] <- ID
# immdata_mixcr$meta[26:61,] <- ED

# immdata$data
immdata_mixcr$meta$Subtype <- c(rep('ACC',15), rep('MEC',20), rep('MECA',10), rep('SDC',16))
immdata_mixcr$meta$Immune <- c(rep('ID',15), rep('ED',20), rep('ID',10), rep('ED',16))

exp_vol <- repExplore(immdata_mixcr$data, .method = "volume"); exp_vol

p1 <- vis(exp_vol, .by = c("Subtype"),
          .meta = immdata_mixcr$meta,
          .color = c("#FEC260","#3FA796","#2A0944","#A10035")); p1
p2 <- vis(exp_vol, .by = c("Immune"), .meta = immdata_mixcr$meta); p2
p1 + p2


TCR_SGC_mat <- SGC_mat[,exp_vol$Sample]
Immune_ssGSEA <- gsva(TCR_SGC_mat,
                      Immune_geneset,
                      min.sz=1,
                      method="ssgsea",
                      ssgsea.norm=F, 
                      max.sz=999999,
                      abs.ranking=F,
                      verbose=T)
IIS_score <- colMeans(Immune_ssGSEA[c(1:18, 21:24),]) # IIS_score, no scaled
TIS_score <- colMeans(Immune_ssGSEA[c(3,15,16,17,18,21,22,23,24),]) # TIS_score, no scaled
CYT_score <- sqrt(TCR_SGC_mat[rownames(TCR_SGC_mat)=='GZMA',] * TCR_SGC_mat[rownames(TCR_SGC_mat)=='PRF1',]) # CYT_Score, no scaled
ANG_score <- Immune_ssGSEA[26,]
APM_score <- Immune_ssGSEA[27,] 

colnames(exp_vol) <- c("Type","TCR")
# exp_vol$Type <- c(rep('ACC',15), rep('MEC',20), rep('MECA',10), rep('SDC',16))
exp_vol$Type <- c(rep('ID',15), rep('ED',20), rep('ID',10), rep('ED',16))
exp_vol$IIS <- as.numeric(scale(IIS_score))
exp_vol$TIS <- as.numeric(scale(TIS_score))
exp_vol$CYT <- as.numeric(scale(CYT_score))
exp_vol$ANG <- as.numeric(scale(ANG_score))
exp_vol$APM <- as.numeric(scale(APM_score))
new_df <- cbind(exp_vol,scale(t(Immune_ssGSEA)))

duct_col <- c("ED"="#AC7D88", "ID"="#97C4B8")
CYT_TCR_cor <- ggscatter(new_df, x = "TCR", y = "TIS", 
                         add = "reg.line", size = 3,
                         cor.coef = T, cor.method = "pearson", legend = 'right', color = 'Type',
                         palette = duct_col,
                         add.params = list(color =  'black', fill = "lightgray"),
                         xlab = "TCR", ylab = "TIS"); CYT_TCR_cor



ACC_idx <- sapply(str_split(rownames(exp_vol), "_T"), function(x) {x[1]}) == "ACC"; exp_vol$Sample[ACC_idx] <- "ACC"
SDC_idx <- sapply(str_split(rownames(exp_vol), "_T"), function(x) {x[1]}) == "SDC"; exp_vol$Sample[SDC_idx] <- "SDC"
MEC_idx <- sapply(str_split(rownames(exp_vol), "_T"), function(x) {x[1]}) == "MEC"; exp_vol$Sample[MEC_idx] <- "MEC"
MECA_idx <- sapply(str_split(rownames(exp_vol), "_T"), function(x) {x[1]}) == "MECA"; exp_vol$Sample[MECA_idx] <- "MECA"

exp_vol$Group <- "ID"; exp_vol$Group[MEC_idx | SDC_idx] <- "ED"
exp_vol$Group <- factor(exp_vol$Group, level = c("ID", "ED"))

CS_p <- ggplot(exp_vol,
               aes(x = Group,
                   y = TCR,
                   fill = Group)) +
          
          # dot
          geom_boxplot(width = 0.5,
                       alpha = 0.5) +
          
          geom_dotplot(method   = "dotdensity",
                       alpha    = 1,
                       binwidth = (max(exp_vol$TCR) - min(exp_vol$TCR))/60,
                       colour   = "black",
                       binaxis  = 'y',
                       stackdir = 'center',
                       position = position_dodge(0.5)) +
          
          # Add an color pallate
          scale_fill_manual(values = c("ID" = "#800080", # red
                                       "ED" = "#FF884B")) +
          
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
          ylim(c(0,200)) +
          
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
                                           color = "black")); CS_p
