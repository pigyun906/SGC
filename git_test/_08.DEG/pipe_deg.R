source("/home/jyhong906/Project/SGC/Script/RNA/git/_07.TIDE/pipe_TIDE.R")

# Make DESeqDataSet & running DESeq
SGC_duct_groups <- as.data.frame(matrix(nrow = nrow(SGC_groups), ncol = 1))
colnames(SGC_duct_groups) <- "Duct"; rownames(SGC_duct_groups) <- rownames(SGC_groups)
SGC_duct_groups[which(SGC_groups$Condition %in% c("ACC", "MECA")),] <- "ID"
SGC_duct_groups[which(SGC_groups$Condition %in% c("SDC", "MEC")),] <- "ED"

SGC_duct_dds <- DESeqDataSetFromMatrix(countData = Keep_SGC_matrix, # 19055 genes
                                  colData = SGC_duct_groups,
                                  design = ~ Duct)
SGC_duct_dds <- DESeq(SGC_duct_dds)

SGC_duct_dat <- counts(SGC_duct_dds, normalized=T)
mod <- model.matrix(~ Duct, colData(SGC_duct_dds))
mod0 <- model.matrix(~ 1, colData(SGC_duct_dds))
SGC_duct_svseq <- svaseq(SGC_duct_dat, mod, mod0, n.sv=2)

SGC_duct_svadds <- SGC_duct_dds
SGC_duct_svadds$SV1 <- SGC_duct_svseq$sv[,1]
SGC_duct_svadds$SV2 <- SGC_duct_svseq$sv[,2]
design(SGC_duct_svadds) <- ~ SV1 + SV2 + Duct
SGC_duct_svadds <- DESeq(SGC_duct_svadds)

SGC_duct_svadds <- SGC_duct_dds
SGC_duct_vsd <- vst(SGC_duct_svadds)
SGC_duct_mat <- assay(SGC_duct_vsd); print("SGC_duct_mat")

# DEG volcano plot -------------------------------------------------------------------------
SGC_duct_svadds_res <- results(SGC_duct_svadds, contrast=c("Duct", "ED", "ID")); SGC_duct_svadds_res
SGC_duct_svadds_sigres <- SGC_duct_svadds_res[which(SGC_duct_svadds_res$padj < 0.01 & abs(SGC_duct_svadds_res$log2FoldChange) > 1),]; SGC_duct_svadds_sigres
SGC_duct_svadds_res <- as.data.frame(dplyr::mutate(as.data.frame(SGC_duct_svadds_res), sig=ifelse(SGC_duct_svadds_res$padj<0.01, "FDR<0.01", "Not Sig")), row.names=rownames(SGC_duct_svadds_res))
SGC_duct_svadds_res[SGC_duct_svadds_res$log2FoldChange > 1,]$sig <- 'ED Originated Carcinoma'
SGC_duct_svadds_res[SGC_duct_svadds_res$log2FoldChange < -1,]$sig <- 'ID Originated Carcinoma'
SGC_duct_svadds_res[SGC_duct_svadds_res$padj > 0.01,]$sig <- 'Not sig'
SGC_duct_svadds_res[SGC_duct_svadds_res$sig!='ID Originated Carcinoma' & SGC_duct_svadds_res$sig!='ED Originated Carcinoma',]$sig <- 'Not Sig'
DEG_volcano_plot <- ggplot(SGC_duct_svadds_res, aes(log2FoldChange, -log10(padj))) +
          geom_point(aes(col = sig), size = 0.7, alpha = 0.5) +
          scale_color_manual(values = c('darkred', 'darkgreen', 'gray')) +
          ggtitle("DEG volcano Plot (padj < 0.01)") +
          ggrepel::geom_text_repel(data=SGC_duct_svadds_res[order(SGC_duct_svadds_res$padj, decreasing = F)[1:15],],
                                   ggplot2::aes(label=rownames(SGC_duct_svadds_res[order(SGC_duct_svadds_res$padj, decreasing = F)[1:15],]))) +
          geom_vline(xintercept=c(-1, 1), col="black", linetype = 'dashed') +
          geom_hline(yintercept=-log10(0.01), col="black", linetype = 'dashed'); DEG_volcano_plot

# DEG #
write.table(rownames(SGC_duct_svadds_res[SGC_duct_svadds_res$log2FoldChange > 1,]), "/home/jyhong906/Project/SGC/Data/RNA/10.DEG/ED UP DEG.txt",
            sep="\t",
            row.names = F,
            col.names = F,
            quote=F,
            append=F,
            na="NA")

write.table(rownames(SGC_duct_svadds_res[SGC_duct_svadds_res$log2FoldChange < -1,]), "/home/jyhong906/Project/SGC/Data/RNA/10.DEG/ID UP DEG.txt",
            sep="\t",
            row.names = F,
            col.names = F,
            quote=F,
            append=F,
            na="NA")


# clusterprofiler #
SGC_duct_svadds_res$Symbol <- rownames(SGC_duct_svadds_res)
SGC_duct_svadds_res <- as.data.frame(SGC_duct_svadds_res)

genelist <- SGC_duct_svadds_res$Symbol
genelist <- genelist %>%  bitr(fromType = "SYMBOL",
                               toType = c("ENSEMBL", "ENTREZID"),
                               OrgDb = org.Hs.eg.db)

colnames(genelist)[1]<-"Symbol"
SGC_duct_svadds_res <- genelist %>%
          inner_join(SGC_duct_svadds_res,by='Symbol') %>% 
          dplyr::select(ENTREZID,
                        log2FoldChange,
                        everything())