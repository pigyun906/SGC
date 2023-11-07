source("/home/jyhong906/Project/SGC/Script/RNA/git/_08.DEG/pipe_deg.R")

SGC_duct_svadds_res$Symbol <- rownames(SGC_duct_svadds_res)
SGC_duct_svadds_res <- as.data.frame(SGC_duct_svadds_res)

genelist <- SGC_duct_svadds_res$Symbol
genelist <- genelist %>%  bitr(fromType = "SYMBOL",
                               toType = c("ENSEMBL", "ENTREZID"),
                               OrgDb = org.Hs.eg.db)

colnames(genelist)[1]<-"Symbol"
SGC_duct_svadds_res <- genelist %>%
          inner_join(SGC_duct_svadds_res,by='Symbol') %>% 
          dplyr::select(ENTREZID,log2FoldChange,everything())

geneList <- SGC_duct_svadds_res$log2FoldChange
names(geneList) <- as.character(SGC_duct_svadds_res$ENTREZID)
geneList <- geneList[!duplicated(names(geneList))]
geneList <- sort(geneList,decreasing=TRUE)

# MSigDB
msigdbr_species()
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
          dplyr::select(gs_name, entrez_gene)
m_c5bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = 'BP') %>% 
          dplyr::select(gs_name, entrez_gene)
m_c2kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = 'KEGG') %>% 
          dplyr::select(gs_name, entrez_gene)

labs <- c("ED up (vs ID down)", "ED down (vs ID up)")
names(labs) <- c("activated", "suppressed")

em <- GSEA(geneList,
           TERM2GENE    = m_t2g,
           pvalueCutoff = 0.01,
           verbose      = FALSE)

for (idx in seq(em@result$ID)) {
          em@result$Description[idx] <- str_split(em@result$Description, pattern = "HALLMARK_")[[idx]][2]
}

dotplot(em,
        split        = ".sign") +
          labs(title = "up HALLMARK") +
          facet_grid(~.sign)

GSEA_C2KEGG <- clusterProfiler::GSEA(geneList,
                                   TERM2GENE    = m_c2kegg,
                                   pvalueCutoff = 0.05,
                                   verbose      = FALSE); save(GSEA_C2KEGG, file = "/home/jyhong906/Project/SGC/Data/RNA/11.GSEA/GSEA_C2KEGG.Rdata")

for (idx in seq(GSEA_C2KEGG@result$ID)) {
          GSEA_C2KEGG@result$Description[idx] <- str_split(GSEA_C2KEGG@result$Description, pattern = "KEGG_")[[idx]][2]
}

write.table(GSEA_C2KEGG@result,
            "/home/jyhong906/Project/SGC/Data/RNA/REVISION/#1-3/GSEA_KEGG.txt",
            sep = "\t",
            quote = F,
            row.names = F)

GSEA_C5BP <- clusterProfiler::GSEA(geneList,
                                   TERM2GENE    = m_c5bp,
                                   pvalueCutoff = 0.05,
                                   verbose      = FALSE); save(GSEA_C5BP, file = "/home/jyhong906/Project/SGC/Data/RNA/11.GSEA/GSEA_C5BP.Rdata")

for (idx in seq(GSEA_C5BP@result$ID)) {
        GSEA_C5BP@result$Description[idx] <- str_split(GSEA_C5BP@result$Description, pattern = "GOBP_")[[idx]][2]
}

# write.table(GSEA_C5BP@result,
#             "/home/jyhong906/Project/SGC/Data/RNA/11.GSEA/GSEA_GOBP.txt",
#             sep = "\t",
#             quote = F,
#             row.names = F)

# enrichplot::dotplot(GSEA_C5BP,
#                     split        = ".sign") +
#         labs(title = "GSEA - Biological process") +
#         facet_grid(~.sign,
#                    labeller = labeller(.sign = labs)) +
#         theme(strip.background = element_blank(),
#               strip.text.y = element_blank(),
#               strip.text.x = element_text(size = 13, color = "black", face = "bold.italic"),
#               title = element_text(size = 15, face = "bold"))

load("/home/jyhong906/Project/SGC/Data/RNA/11.GSEA/GSEA_C5BP.Rdata")
ID_geneSetID <- which(GSEA_C5BP@result$Description %in% c("POSITIVE_REGULATION_OF_STEM_CELL_DIFFERENTIATION",
                                                          "HIPPO_SIGNALING",
                                                          "SMOOTHENED_SIGNALING_PATHWAY",
                                                          "STEM_CELL_DIFFERENTIATION",
                                                          "POSITIVE_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY",
                                                          "REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY",
                                                          "REGULATION_OF_WNT_SIGNALING_PATHWAY",
                                                          "CANONICAL_WNT_SIGNALING_PATHWAY",
                                                          "CELL_CELL_SIGNALING_BY_WNT",
                                                          "REGULATION_OF_STEM_CELL_DIFFERENTIATION")); length(ID_geneSetID)
ED_geneSetID <- which(GSEA_C5BP@result$Description %in% c("REGULATION_OF_DENDRITIC_CELL_DIFFERENTIATION",
                                                          "REGULATION_OF_DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION",
                                                          "DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION",
                                                          "NATURAL_KILLER_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
                                                          "POSITIVE_REGULATION_OF_REGULATORY_T_CELL_DIFFERENTIATION",
                                                          "T_CELL_ACTIVATION_VIA_T_CELL_RECEPTOR_CONTACT_WITH_ANTIGEN_BOUND_TO_MHC_MOLECULE_ON_ANTIGEN_PRESENTING_CELL",
                                                          "REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION",
                                                          "T_CELL_MEDIATED_CYTOTOXICITY",
                                                          "POSITIVE_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                                                          "POSITIVE_REGULATION_OF_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY"))

# Save in /home/jyhong906/Project/SGC/Data/RNA/11.GSEA/Fig. ID GSEA plot.pdf (5 x 7)
gseaplot2(GSEA_C5BP,
          geneSetID = ID_geneSetID,
          pvalue_table = TRUE,
          base_size = 15,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1,
          ES_geom = "line",
          color = rep("#B99B6B", 10))

# Save in /home/jyhong906/Project/SGC/Data/RNA/11.GSEA/Fig. ED GSEA plot.pdf (5 x 7)
gseaplot2(GSEA_C5BP,
          geneSetID = ED_geneSetID,
          pvalue_table = TRUE,
          base_size = 15,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1,
          ES_geom = "line",
          color = rep("#7F669D", 10))

# Combined #
gseaplot2(GSEA_C5BP,
          geneSetID = c(ID_geneSetID, ED_geneSetID),
          pvalue_table = TRUE,
          base_size = 15,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1,
          ES_geom = "line",
          color = c(rep("#800080", 10),
                    rep("#FF884B", 10)))

library(metapod)
# Wnt #
ComFDR <- combineParallelPValues(list(0.001021, 0.001368, 0.002597, 0.002644, 0.04),
                                 method = c("fisher")); ComFDR$p.value

# Hippo #
ComFDR <- combineParallelPValues(list(0.001021, 0.001368, 0.002597, 0.002644, 0.04),
                                 method = c("fisher")); ComFDR$p.value

# Smoothened #
ComFDR <- combineParallelPValues(list(0.001021, 0.001368, 0.002597, 0.002644, 0.04),
                                 method = c("fisher")); ComFDR$p.value

# Stem cell #
ComFDR <- combineParallelPValues(list(0.003708, 0.01121, 0.03532),
                                 method = c("fisher")); ComFDR$p.value

# c2kegg_em <- GSEA(geneList,
#                   nPerm        = 10000,
#                   TERM2GENE    = m_c2kegg,
#                   pvalueCutoff = 0.01,
#                   verbose      = FALSE)

# dotplot(c2kegg_em,
#         showCategory = 100,
#         split        = ".sign") +
#           facet_grid(~.sign)
# 
# ridgeplot(c2kegg_em) +
#         labs(x = "enrichment distribution")
