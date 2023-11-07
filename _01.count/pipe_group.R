pkg_control(c("sva", "ggrepel"))

# Batch groups
SGC_groups <- as.data.frame(matrix(nrow=length(colnames(SGC)), ncol=2))
rownames(SGC_groups) <- colnames(SGC)
colnames(SGC_groups) <- c("Condition", "Batch")

ACC_idx <- substr(rownames(SGC_groups), 1, 4) == "ACC_"
MECA_idx <- substr(rownames(SGC_groups), 1, 4) == "MECA"
SDC_idx <- substr(rownames(SGC_groups), 1, 4) == "SDC_"
MEC_idx <- substr(rownames(SGC_groups), 1, 4) == "MEC_"

SGC_groups$Condition[ACC_idx] <- "ACC"
SGC_groups$Condition[MECA_idx] <- "MECA"
SGC_groups$Condition[SDC_idx] <- "SDC"
SGC_groups$Condition[MEC_idx] <- "MEC"
SGC_groups$Condition <- factor(SGC_groups$Condition, levels = c("ACC", "MECA", "SDC", "MEC"))

SGC_groups[ACC_idx,]$Batch <- "1"; SGC_groups[MECA_idx,]$Batch <- "1"; SGC_groups[SDC_idx,]$Batch <- "1"
SGC_groups[MEC_idx,]$Batch <- "2"

# Filter genes
Keep_SGC_matrix <- SGC[rowSums(SGC[,substr(colnames(SGC),1,4) == "ACC_"] >= 1) >= 4 &   # 1/2 of Sample >= each 1 counts  # ACC
                                 rowSums(SGC[,substr(colnames(SGC),1,4) == "MECA"] >= 1) >= 20 & # MECA 
                                 rowSums(SGC[,substr(colnames(SGC),1,4) == "SDC_"] >= 1) >= 8 & # SDC
                                 rowSums(SGC[,substr(colnames(SGC),1,4) == "MEC_"] >=1 ) >= 10,];Keep_SGC_matrix # MEC