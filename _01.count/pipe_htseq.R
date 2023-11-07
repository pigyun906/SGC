pkg_control(c("fgsea","GSVA","ggpubr","org.Hs.eg.db","clusterProfiler","viridis","circlize","ComplexHeatmap","stringr",
              "msigdbr","ggridges","enrichplot","tidyverse","ggplot2","grid","gridExtra","DESeq2","ggdendro","dplyr","dendextend",
              "GeoTcgaData", "immunedeconv","gelnet","biomaRt","synapser","immunarch","ggpol","psych","ggcorrplot","stringr",
              "Biobase","Seurat","AnnotationHub","EnsDb.Hsapiens.v86","devtools","patchwork","BayesPrism","reshape","InstaPrism"))


ACC_dir <- c("/home/jyhong906/Project/SGC/Data/RNA/00.HTSeq/ACC")
SDC_dir <- c("/home/jyhong906/Project/SGC/Data/RNA/00.HTSeq/SDC") 
MEC_dir <- c("/home/jyhong906/Project/SGC/Data/RNA/00.HTSeq/MEC")
MECA_dir <- c("/home/jyhong906/Project/SGC/Data/RNA/00.HTSeq/MECA") 

ACC_file <- list.files(ACC_dir)
SDC_file <- list.files(SDC_dir)
MEC_file <- list.files(MEC_dir)
MECA_file <- list.files(MECA_dir)


ACC_files <- gsub(".htseq.count.txt", "", ACC_file)
SDC_files <- gsub(".htseq.count.txt", "", SDC_file)
MEC_files <- gsub(".htseq.count.txt", "", MEC_file)
MECA_files <- gsub(".htseq.count.txt", "", MECA_file)

ACC_list <- list()
SDC_list <- list()
MEC_list <- list()
MECA_list <- list()

ACC_list[[1]] <- read.table(file.path(ACC_dir, ACC_file[1]),
                                header = F,
                                sep = "\t",
                                stringsAsFactors = F)
SDC_list[[1]] <- read.table(file.path(SDC_dir, SDC_file[1]),
                                header = F,
                                sep = "\t",
                                stringsAsFactors = F)
MEC_list[[1]] <- read.table(file.path(MEC_dir, MEC_file[1]),
                            header = F,
                            sep = "\t",
                            stringsAsFactors = F)
MECA_list[[1]] <- read.table(file.path(MECA_dir, MECA_file[1]),
                            header = F,
                            sep = "\t",
                            stringsAsFactors = F)

Gene_list <- ACC_list[[1]][,c(1)]

for( i in 1:length(ACC_file)){
          ACC_list[[i]] <- read.table(file.path(ACC_dir, ACC_file[i]),
                                          header = F,
                                          sep = '\t',
                                          stringsAsFactors = F)
          ACC_list[[i]] <- as.data.frame(ACC_list[[i]][,ncol(ACC_list[[i]])])
          colnames(ACC_list[[i]]) <- ACC_files[i]
}
for( i in 1:length(SDC_file)){
          SDC_list[[i]] <- read.table(file.path(SDC_dir, SDC_file[i]),
                                      header = F,
                                      sep = '\t',
                                      stringsAsFactors = F)
          SDC_list[[i]] <- as.data.frame(SDC_list[[i]][,ncol(SDC_list[[i]])])
          colnames(SDC_list[[i]]) <- SDC_files[i]
}
for( i in 1:length(MEC_file)){
          MEC_list[[i]] <- read.table(file.path(MEC_dir, MEC_file[i]),
                                      header = F,
                                      sep = '\t',
                                      stringsAsFactors = F)
          MEC_list[[i]] <- as.data.frame(MEC_list[[i]][,ncol(MEC_list[[i]])])
          colnames(MEC_list[[i]]) <- MEC_files[i]
}
for( i in 1:length(MECA_file)){
          MECA_list[[i]] <- read.table(file.path(MECA_dir, MECA_file[i]),
                                      header = F,
                                      sep = '\t',
                                      stringsAsFactors = F)
          MECA_list[[i]] <- as.data.frame(MECA_list[[i]][,ncol(MECA_list[[i]])])
          colnames(MECA_list[[i]]) <- MECA_files[i]
}

SGC <- as.data.frame(c(ACC_list, MECA_list, SDC_list, MEC_list))
rownames(SGC) <- Gene_list
Tail_removed_SGC <- SGC[1:c(length(rownames(SGC))-5),]
SGC <- Tail_removed_SGC; SGC_countData <- SGC; sum(is.na(SGC_countData))

SGC_duct_groups <- data.frame(
          
          "Duct" = c(rep("ID", sum(sapply(strsplit(colnames(SGC), "_T_"), function(x) {x[1]}) %in% c("ACC", "MECA"))),
                     rep("ED", sum(sapply(strsplit(colnames(SGC), "_T_"), function(x) {x[1]}) %in% c("SDC", "MEC"))))
          
          ); rownames(SGC_duct_groups) <- colnames(SGC)

print("Successfully make SGC count matrix, which called SGC")

source("/home/jyhong906/Project/SGC/Script/RNA/git/_1.count/pipe_group.R")