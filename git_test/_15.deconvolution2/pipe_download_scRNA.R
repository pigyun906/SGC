source("/home/jyhong906/Project/SGC/Script/RNA/git/_14.APCs/pipe_APCs.R")

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


SGC_scftp <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210171/suppl/GSE210171_acc_scrnaseq_counts.txt.gz"
SGC_scCounts <- data.table::fread(SGC_scftp); SGC_scCounts <- SGC_scCounts[substr(SGC_scCounts$Geneid,1,4) == "ENSG",]; SGC_scCounts <- as.data.frame(SGC_scCounts)
GeneID <- SGC_scCounts$Geneid; SGC_scCounts <- SGC_scCounts[,which(substr(colnames(SGC_scCounts),1,3) == "ACC")]; SGC_scCounts <- t(SGC_scCounts); colnames(SGC_scCounts) <- GeneID
rownames(SGC_scCounts) <- sapply(str_split(rownames(SGC_scCounts),".sort.bam.gene_counts.tsv"), function(x) {x[1]})
SGC_scCounts_filter <- cleanup.genes (input=SGC_scCounts,
                                      input.type="count.matrix",
                                      species="hs",
                                      gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                      exp.cells=5)

# select the human dataset from Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# define the input Ensembl IDs
ensembl_ids <- sapply(str_split(colnames(SGC_scCounts_filter), "\\."), function(x) {x[1]})

# retrieve the corresponding HGNC symbols using the biomaRt package
hgnc_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                      filters = "ensembl_gene_id", 
                      values = ensembl_ids, 
                      mart = ensembl);

# create a named vector that maps Ensembl gene IDs to HGNC symbols
ensembl_to_hgnc <- setNames(hgnc_symbols$hgnc_symbol, hgnc_symbols$ensembl_gene_id)

# use the named vector to rename the columns of SGC_scCounts_filter
colnames(SGC_scCounts_filter) <- ensembl_to_hgnc[sapply(str_split(colnames(SGC_scCounts_filter), "\\."), function(x) {x[1]})]
SGC_scCounts_filter <- SGC_scCounts_filter[,!c(is.na(colnames(SGC_scCounts_filter)) | colnames(SGC_scCounts_filter) == "")]
SGC_scCounts_filter <- SGC_scCounts_filter[,!duplicated(colnames(SGC_scCounts_filter))]

scSGC_object <- CreateSeuratObject(counts = t(SGC_scCounts_filter),
                                   project = "scSGC",
                                   min.cells = 5,
                                   min.features = 200)

VlnPlot(
          scSGC_object,
          features = c("nFeature_RNA", "nCount_RNA"),
          ncol = 2
); FeatureScatter(scSGC_object, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")

# Normalize data
scSGC_object <- NormalizeData(scSGC_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable genes
scSGC_object <- FindVariableFeatures(scSGC_object, selection.method = "vst")

# Scale data
scSGC_object <- ScaleData(scSGC_object)

# Perform PCA
scSGC_object <- RunPCA(scSGC_object)

# Score the JackStraw procedure
scSGC_object <- JackStraw(scSGC_object, num.replicate = 100)
scSGC_object <- ScoreJackStraw(scSGC_object, dims = 1:20)

# Generate Elbow plot
elbow_plot <- ElbowPlot(scSGC_object); elbow_plot
dims_use <- 13

# Run FindNeighbors using optimal number of dimensions
scSGC_object <- FindNeighbors(scSGC_object, dims = 1:dims_use, k.param = 15)
scSGC_object <- FindClusters(scSGC_object, resolution = 0.5)

# Run clustering using "TSNE"
scSGC_object<- RunTSNE(object = scSGC_object, dims = 1:dims_use)
DimPlot(object = scSGC_object, reduction = "tsne")

# create a named list of markers for each cell type
SG_cell_type_markers_list <- list(
          "Ductal/basal epithelial cell" = c("KRT15", "SOX2"),
          "Fibroblast" = c("DCN", "LUM"),
          "Duct cell" = c("100A2", "WFDC2"),
          "Lonocyte" = c("CFTR", "FOXI1"),
          "Mucous acini" = c("MUC5B", "BPIFB2"),
          "Pericyte" = "MYO1B",
          "Serous acini cell" = c("LPO", "ODAM"),
          "Artery" = "CLDN5",
          "Capillary" = "CA4",
          "Venule" = "AQP1",
          "Endothelial cells" = c("PECAM1", "ENG", "VWF"),
          "T/NK cell" = c("GZMA", "HCST"),
          "Smooth muscle cells" = c("ACTA2", "MYH11", "ACTG2", "CNN1", "CALD1", "MCAM", "TAGLN", "PDGFRB", "MYL9"),
          "Myeloid cell" = c("AIF1", "CD163", "LYZ"),
          "Muscle satellite cell" = c("PAX7", "CD82", "NCAM1", "MYF5"),
          "Mast cell" = c("MS4A2", "TPSB2", "GATA2"),
          "Skeletal muscle cell" = c("ACTA1", "NEB", "MYL2"),
          "Lymphatic endothelial cell" = c("PDPN", "PROX1", "LYVE1"),
          "Schwann cell" = c("NGFR", "SOX10", "GAP43", "CDH19"),
          "Plasma cells" = c("JCHAIN", "IGKC", "IGHA1", "IGHA2"),
          "T cell" = c("CD3G", "CD3D", "CD3E"),
          "B cell" = c("MS4A1", "BANK1", "CD37"),
          "Intercalated duct cell" = c("KIT", "AQP5","ANO1","SOX10"),
          "Intercalated duct cells" = c("ANO1", "SOX10"),
          "Myoepithelial cell" = c("ACTA2", "MYH11", "CNN1"),
          "Stem cells" = c("ALDH1A1", "CD44", "PROM1", "NANOG", "KIT", "NES", "KLF4", "CD55", "ALCAM", "NOTCH4", "WNT7A", "PDPN"),
          "Stem cells" = c("POU5F1","SOX2","MYC","KLF4","NANOG","ZFP42","PROM1","CTNNB1","ABCG2","CD34","CD44","ZSCAN4"),
          "Tc1" = c("D8A", "EOMES", "GMZA"),
          "Tc2" = c("CD8A", "SPRY1", "ITAG1"),
          "T helper cell" = c("CD40LG"),
          "Dendritic cell" = c("CCL13","CCL17","CCL22","CD209","HSD11B1","NPR1","PPFIBP2"),
          "Macrophage" = c("APOE","ATG7","BCAT1","CCL7","CD163","CD68","CD84","CHI3L1","CHIT1","CLEC5A","COL8A2","COLEC12","CTSK","CXCL5","CYBB","DNASE2B","EMP1","FDX1","FN1","GM2A","GPC4","KAL1","MARCO","ME1","MS4A4A","MSR1","PCOLCE2","PTGDS","RAI14","SCARB2","SCG5","SGMS1","SULT1C2")

); gs_list <- list(gs_positive = SG_cell_type_markers_list)

# load auto-detection function
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
# db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";

# guess a tissue type
# tissue_guess = auto_detect_tissue_type(path_to_db_file = db_,
#                                        seuratObject = scSGC_object,
#                                        scaled = TRUE,
#                                        assay = "RNA")  # if saled = TRUE, make sure the data is scaled.
# 
# get cell-type-specific gene sets from our in-built database (DB)
# gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Kidney")
# 
# # assign cell types
# # scRNAseqData = readRDS(gzcon(url('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/exampleData.RDS'))); #load example scRNA-seq matrix

es.max = sctype_score(scRNAseqData = t(SGC_scCounts_filter), scaled = TRUE, gs = gs_list$gs_positive)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(scSGC_object@meta.data$seurat_clusters), function(cl) {
          es.max.cl = sort(rowSums(es.max[ ,rownames(scSGC_object@meta.data[scSGC_object@meta.data$seurat_clusters == cl, ])]), decreasing = !0)
          head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scSGC_object@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
scSGC_object@meta.data$customclassif = ""

for(j in unique(sctype_scores$cluster)){
          cl_type = sctype_scores[sctype_scores$cluster==j,]; 
          scSGC_object@meta.data$customclassif[scSGC_object@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(scSGC_object, reduction = "tsne", label = TRUE, repel = TRUE, group.by = 'customclassif')