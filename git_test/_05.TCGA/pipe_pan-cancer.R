source("/home/jyhong906/Project/SGC/Script/RNA/git/_04.TIME/pipe_TIME.score.R")

CDG_TCGA_dir <- c("/home/jyhong906/Public/GDC/TCGA") # Set directory for MEC
GDC_TCGA_list <- list() # Make GDC_TCGA_list 
GDC_TCGA_files <- list.files(CDG_TCGA_dir) # Make MEC_file 
Correct_GDC_TCGA_files <- gsub(".txt", "", GDC_TCGA_files) 

for (idx in seq_along(Correct_GDC_TCGA_files)) {
          GDC_TCGA_list[[idx]] <- read.table(file.path(CDG_TCGA_dir,GDC_TCGA_files[idx]),
                                             header = T,
                                             sep = '\t')
          GDC_TCGA_list[[idx]] <- as.data.frame(GDC_TCGA_list[[idx]][,1:ncol(GDC_TCGA_list[[idx]])])
          print(GDC_TCGA_files[idx])
          for (idx1 in seq_along(GDC_TCGA_list[[idx]])){
                    colnames(GDC_TCGA_list[[idx]])[idx1] <- paste0(Correct_GDC_TCGA_files[idx],idx1)
          }
}

TCGA_pan_cancer_count_df <- as.data.frame(GDC_TCGA_list)
int_names <- intersect(rownames(TCGA_pan_cancer_count_df), rownames(Keep_SGC_matrix))
TCGA_pan_cancer_count_df <- TCGA_pan_cancer_count_df[int_names,]
TCGA_pan_cancer_count_df <- TCGA_pan_cancer_count_df[rownames(Keep_SGC_matrix),]; rownames(TCGA_pan_cancer_count_df) <- rownames(Keep_SGC_matrix)
Pan_cancer_with_Keep_SGC_matrix <- cbind(TCGA_pan_cancer_count_df, Keep_SGC_matrix)

TCGA_colData <- as.data.frame(matrix(nrow=length(colnames(Pan_cancer_with_Keep_SGC_matrix)), ncol=1))
rownames(TCGA_colData) <- colnames(Pan_cancer_with_Keep_SGC_matrix)
colnames(TCGA_colData) <- 'TumorType'

for (idx in Correct_GDC_TCGA_files) {
          if (substr(idx,1,4) == 'TCGA') {
                    print(idx)
                    print(grep(idx,rownames(TCGA_colData)))
                    TCGA_colData[grep(idx,rownames(TCGA_colData)),] <- idx            
          }
}

grep('ACC',rownames(TCGA_colData), fixed = T)
TCGA_colData[11094:11111,] <- 'SGC_ACC'
TCGA_colData[grep('MECA',rownames(TCGA_colData)),] <- 'SGC_MECA'
TCGA_colData[grep('SDC',rownames(TCGA_colData)),] <- 'SGC_SDC'
grep('MEC',rownames(TCGA_colData), fixed = T)
rownames(TCGA_colData)[grep('MEC',rownames(TCGA_colData))]
TCGA_colData[11168:11187,] <- 'SGC_MEC'
TCGA_colData$TumorType <- factor(TCGA_colData$TumorType)
Pan_cancer_with_SGC_countmatrix <- as.matrix(Pan_cancer_with_Keep_SGC_matrix)

SGC_files <- c('ACC', 'MECA', 'SDC', 'MEC'); SGC_files
Correct_GDC_TCGA_files[34:37] <- SGC_files
new_files <- Correct_GDC_TCGA_files

Pan_cancer_with_SGC_countmatrix[is.na(Pan_cancer_with_SGC_countmatrix)] <- 0
TCGA_dds <- DESeqDataSetFromMatrix(countData = Pan_cancer_with_SGC_countmatrix,
                              colData = TCGA_colData,
                              design = ~ TumorType)
keep <- rowSums(counts(TCGA_dds)) >= length(rownames(TCGA_colData))
TCGA_dds <- TCGA_dds[keep,]
TCGA_dds <- estimateSizeFactors(TCGA_dds, type='poscounts'); TCGA_dds
TCGA_pan_cancer_mat <- assay(vst(TCGA_dds)); TCGA_pan_cancer_mat

write.table(TCGA_pan_cancer_mat, "/home/jyhong906/Project/SGC/Data/RNA/05.TCGA pan-cancer/TCGA_pan_cancer_mat.txt",
            sep="\t",
            row.names = T,
            col.names = T,
            quote=F,
            append=F,
            na="NA")

filterCommonGenes(input.f = '/home/jyhong906/Project/SGC/Data/RNA/05.TCGA pan-cancer/TCGA_pan_cancer_mat.txt', # Filtering for common genes
                  output.f = '/home/jyhong906/Project/SGC/Data/RNA/05.TCGA pan-cancer/TCGA_pan_cancer_estimate_mat.gct',
                  id = 'GeneSymbol')
estimateScore('/home/jyhong906/Project/SGC/Data/RNA/05.TCGA pan-cancer/TCGA_pan_cancer_estimate_mat.gct',
              '/home/jyhong906/Project/SGC/Data/RNA/05.TCGA pan-cancer/TCGA_pan_cancer_estimate_score.gct',
              platform = 'illumina') # Calculate score

TCGA_pan_cancer_estimate_score <- read.table('/home/jyhong906/Project/SGC/Data/RNA/05.TCGA pan-cancer/TCGA_pan_cancer_estimate_score.gct',
                                             header=T,
                                             sep="\t",
                                             skip = 2,
                                             row.names = 1,
                                             stringsAsFactors=F); TCGA_pan_cancer_estimate_score
TCGA_pan_cancer_estimate_immune_score <- TCGA_pan_cancer_estimate_score[2,-1]


#===================================#
# Boxplot for Estimate immune score #
#####################################
TCGA_pan_cancer_estimate_immune_scaled_score <- as.data.frame(t(t(scale(t(TCGA_pan_cancer_estimate_immune_score), center = T, scale = T))))
# Group은 TCGA, SGC_each subtype으로 설정.
TCGA_pan_cancer_estimate_immune_scaled_score <- TCGA_pan_cancer_estimate_immune_scaled_score %>%
          mutate(Cohort = ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),1,4)=='TCGA','TCGA','SGC'))
TCGA_pan_cancer_estimate_immune_scaled_score <- TCGA_pan_cancer_estimate_immune_scaled_score %>%
          mutate(Type = ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,8)=='ACC','ACC',
                               ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='BLCA','BLCA',
                                      ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='BRCA','BRCA',
                                             ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='CESC','CESC',
                                                    ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='CHOL','CHOL',
                                                           ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='COAD','COAD',
                                                                  ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='DLBC','DLBC',
                                                                         ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='ESCA','ESCA',
                                                                                ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,8)=='GBM','GBM',
                                                                                       ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='HNSC','HNSC',
                                                                                              ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='KICH','KICH',
                                                                                                     ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='KIRC','KIRC',
                                                                                                            ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='KIRP','KIRP',
                                                                                                                   ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='LAML','LAML',
                                                                                                                          ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,8)=='LGG','LGG',
                                                                                                                                 ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='LIHC','LIHC',
                                                                                                                                        ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='LUAD','LUAD',
                                                                                                                                               ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='LUSC','LUSC',
                                                                                                                                                      ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='MESO','MESO',
                                                                                                                                                             ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,7)=='OV','OV',
                                                                                                                                                                    ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='PAAD','PAAD',
                                                                                                                                                                           ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='PCPG','PCPG',
                                                                                                                                                                                  ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='PRAD','PRAD',
                                                                                                                                                                                         ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='READ','READ',
                                                                                                                                                                                                ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='SARC','SARC',
                                                                                                                                                                                                       ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='SKCM','SKCM',
                                                                                                                                                                                                              ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='STAD','STAD',
                                                                                                                                                                                                                     ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='TGCT','TGCT',
                                                                                                                                                                                                                            ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='THCA','THCA',
                                                                                                                                                                                                                                   ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='THYM','THYM',
                                                                                                                                                                                                                                          ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,9)=='UCEC','UCEC',
                                                                                                                                                                                                                                                 ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,8)=='UCS','UCS',
                                                                                                                                                                                                                                                        ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),6,8)=='UVM','UVM',
                                                                                                                                                                                                                                                               ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),1,3)=='ACC','*ACC',
                                                                                                                                                                                                                                                                      ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),1,4)=='MECA','*MECA',
                                                                                                                                                                                                                                                                             ifelse(substr(rownames(TCGA_pan_cancer_estimate_immune_scaled_score),1,3)=='SDC','*SDC','*MEC')))))))))))))))))))))))))))))))))))))
# Identify low, intermediate, high-immune score by tertile #
Pan_type <- unique(TCGA_pan_cancer_estimate_immune_scaled_score$Type)[1:33]

median_list <- c()
for (idx in Pan_type) {
          print(idx)
          print(median(TCGA_pan_cancer_estimate_immune_scaled_score[TCGA_pan_cancer_estimate_immune_scaled_score$Type == idx,]$ImmuneScore))
          median_list <- c(median_list, median(TCGA_pan_cancer_estimate_immune_scaled_score[TCGA_pan_cancer_estimate_immune_scaled_score$Type == idx,]$ImmuneScore))
}
order_type <- Pan_type[mixedorder(median_list)]; test_list <- median_list[mixedorder(median_list)]; names(test_list) <- order_type; quantile(median_list, probs = seq(0, 1, by = .01))

tertile_limits <- quantile(TCGA_pan_cancer_estimate_immune_scaled_score$ImmuneScore[1:11093], seq(0, 1, 1/2), na.rm = TRUE)
New_type <- cut(median_list, tertile_limits, c('TCGA-Low', "TCGA-High"), include.lowest = TRUE)


for (idx in seq(Pan_type)) {
          TCGA_pan_cancer_estimate_immune_scaled_score$Cohort[which(TCGA_pan_cancer_estimate_immune_scaled_score$Type == Pan_type[idx])] <- as.character(New_type[idx])
}
TCGA_pan_cancer_estimate_immune_scaled_score$Cohort[which(TCGA_pan_cancer_estimate_immune_scaled_score$Cohort == "SGC")] <- TCGA_pan_cancer_estimate_immune_scaled_score$Type[which(TCGA_pan_cancer_estimate_immune_scaled_score$Cohort == "SGC")]

Pan_cols <- c("TCGA-Low"="#839AA8","TCGA-Intermediate"="#CDC2AE","TCGA-High"="#D1512D","*ACC"="#FFB200","*MECA"="#3D8361","*SDC"="#25316D","*MEC"="#D2001A")
p <- ggplot(TCGA_pan_cancer_estimate_immune_scaled_score, aes(x=Type, y=ImmuneScore, fill=Cohort)) +
          scale_x_discrete(limits=c(order_type, c('*ACC','*MECA','*SDC','*MEC'))) +
          geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=Cohort), alpha=0.9, size = 0.01) +
          geom_boxplot(alpha = 0.5, show.legend = F, outlier.shape = NA) +
          stat_compare_means(method = "wilcox.test",
                             symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                              symbols = c("****", "***", "**", "*", "ns")), label.y = 6.5, label.x = 1.3) +
          geom_signif(comparisons = list(c("*ACC", "*MECA")),
                      map_signif_level=T, textsize=2, y_position = 2) +
          geom_signif(comparisons = list(c("*MECA", "*SDC")),
                      map_signif_level=T, textsize=2, y_position = 2.4) +
          geom_signif(comparisons = list(c("*SDC", "*MEC")),
                      map_signif_level=T, textsize=2, y_position = 2.8) +
          theme(axis.text.x = element_text(size=10, color="black", angle = 90, hjust=1),
                axis.title.x = element_text(size=15, face="bold"),
                axis.title.y = element_text(size=15, face="bold")) +
          
          # color
          scale_color_manual(values = Pan_cols) + # line
          scale_fill_manual(values = Pan_cols) + # fill
          
          # axis, main title
          labs(x     = "Tumor Type",
               y     = "Normalized score",
               title = "Immune score of Salivary Gland Carcinoma with TCGA pan cancer") +
          
          # white backgroud
          theme_bw() +
          
          #eliminates background, gridlines, and chart border
          theme(axis.title       = element_text(size = 15,
                                                face = 'plain'),
                title            = element_text(size = 20,
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
          
          #draws x and y axis line
          theme(axis.line = element_line(color = 'black'),
                axis.text.x = element_text(size = 5,
                                           color = "black"),
                axis.text.y = element_text(size = 15,
                                           color = "black")); p


kruskal.test(ImmuneScore ~ Type, data = TCGA_pan_cancer_estimate_immune_scaled_score[TCGA_pan_cancer_estimate_immune_scaled_score$Type %in% c("*ACC","*MECA","*SDC","*MEC"),])
test_df <- TCGA_pan_cancer_estimate_immune_scaled_score[TCGA_pan_cancer_estimate_immune_scaled_score$Type %in% c("*ACC","*MECA","*SDC","*MEC"),]
test_df$Type <- c(rep("ID",58), rep("ED",36))
wilcox.test(ImmuneScore ~ Type, data = test_df)

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

Pan_immune_ssGSEA <- gsva(TCGA_pan_cancer_mat,
                          Immune_geneset,
                          min.sz=1,
                          method="ssgsea",
                          ssgsea.norm=F,
                          max.sz=999999,
                          abs.ranking=F,
                          verbose=T)

Pan_immune_df <- as.data.frame(matrix(nrow = length(colnames(Pan_immune_ssGSEA)), ncol = 5))
rownames(Pan_immune_df) <- colnames(Pan_immune_ssGSEA)
names(Pan_immune_df) <- c('CD8', 'CTL', 'CYT', 'Tcell', 'Type')
Pan_immune_df$CD8 <- as.numeric(Pan_immune_ssGSEA[3,])
Pan_immune_df$CTL <- as.numeric(Pan_immune_ssGSEA[4,])
Pan_immune_df$Tcell <- as.numeric(Pan_immune_ssGSEA[15,])

CYT_score <- sqrt(TCGA_pan_cancer_mat[rownames(TCGA_pan_cancer_mat)=='GZMA',] * TCGA_pan_cancer_mat[rownames(TCGA_pan_cancer_mat)=='PRF1',])
Pan_immune_df$CYT <- as.numeric(CYT_score)

for (i in 1:length(rownames(Pan_immune_df))) {
          real_i_name <- gsub('\\d','', rownames(Pan_immune_df))[i]
          if (substr(gsub('\\d','', rownames(Pan_immune_df))[i],0,4) <= 'TCGA') {
                    Pan_immune_df$Type[i] <- gsub('\\d','', rownames(Pan_immune_df))[i]
                    print(i)
                    print(real_i_name)}}
for (i in 1:length(rownames(Pan_immune_df))){
          Pan_immune_df$Type[i] <- gsub('\\d','', Pan_immune_df$Type[i])
}

Pan_immune_df$Type[11094:11111] <- 'ACC'
pan_var_uq <- unique(Pan_immune_df$Type)

Pan_immune_df2 <- as.data.frame(matrix(nrow=37, ncol = 6))
names(Pan_immune_df2) <- c('CD8_median', 'CTL_median', 'CYT_median', 'T_median', 'Group', 'Type')
Pan_immune_df2$Type <- rownames(Pan_immune_df2)
rownames(Pan_immune_df2) <- pan_var_uq
Pan_immune_df2$Group <- 'TCGA'
Pan_immune_df2$Type <- substr(rownames(Pan_immune_df2), 6, 10)
Pan_immune_df2$Type[34:37] <- c('ACC', 'MECA', 'SDC', 'MEC')
Pan_immune_df2$Type[c(34,35)] <- c('ACC', 'MECA')
Pan_immune_df2$Type[c(36,37)] <- c('SDC', 'MEC')

for (idx2 in pan_var_uq){
          if (idx2 != 'ACC' & idx2 != 'MEC') {
                    print(idx2)
                    Pan_immune_df2$CD8_median[grep(idx2, pan_var_uq)] <- median(Pan_immune_df$CD8[grep(idx2, Pan_immune_df$Type)])
                    Pan_immune_df2$CTL_median[grep(idx2, pan_var_uq)] <- median(Pan_immune_df$CTL[grep(idx2, Pan_immune_df$Type)])
                    Pan_immune_df2$CYT_median[grep(idx2, pan_var_uq)] <- median(Pan_immune_df$CYT[grep(idx2, Pan_immune_df$Type)])
                    Pan_immune_df2$T_median[grep(idx2, pan_var_uq)] <- median(Pan_immune_df$Tcell[grep(idx2, Pan_immune_df$Type)])
          }
          else if (idx2 == 'ACC') {
                    print(idx2)
                    Pan_immune_df2$CD8_median[34] <- median(Pan_immune_df$CD8[11094:11111])
                    Pan_immune_df2$CTL_median[34] <- median(Pan_immune_df$CTL[11094:11111])
                    Pan_immune_df2$CYT_median[34] <- median(Pan_immune_df$CYT[11094:11111])
                    Pan_immune_df2$T_median[34] <- median(Pan_immune_df$Tcell[11094:11111])
          }
          else if (idx2 == 'MEC') {
                    print(idx2)
                    Pan_immune_df2$CD8_median[37] <- median(Pan_immune_df$CD8[11168:11187])
                    Pan_immune_df2$CTL_median[37] <- median(Pan_immune_df$CTL[11168:11187])
                    Pan_immune_df2$CYT_median[37] <- median(Pan_immune_df$CYT[11168:11187])
                    Pan_immune_df2$T_median[37] <- median(Pan_immune_df$Tcell[11168:11187])
          }
}
Pan_immune_df2$CD8_median <- scale(Pan_immune_df2$CD8_median)
Pan_immune_df2$CTL_median <- scale(Pan_immune_df2$CTL_median)
Pan_immune_df2$CYT_median <- scale(Pan_immune_df2$CYT_median)
Pan_immune_df2$T_median <- scale(Pan_immune_df2$T_median)

Pan_immune_df2$Group <- c(as.character(New_type), c("*ACC", "*MECA", "*SDC", "*MEC"))
TCGA_T_CYT_cor <- ggscatter(Pan_immune_df2, x = "T_median", y = "CYT_median", size = 5,
                              add = "reg.line",
                              cor.coef = T, cor.method = "pearson", color = 'Group', legend = 'right',
                              add.params = list(color =  'black'),
                              palette = Pan_cols,
                              xlab = "T cell fraction", ylab = "Cytolytic Activity Score") +
          geom_text_repel(aes(label = Type), size = 5, hjust= -0.5) +#, hjust = -0.5, vjust = 0.5) +
          font("xlab", size = 20, vjust = -1)+
          font("ylab", size = 20, vjust = +2); TCGA_T_CYT_cor