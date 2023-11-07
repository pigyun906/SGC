sc_Expr <- t(SGC_scCounts_filter[rownames(scSGC_object@meta.data),])
cell_type_labels <- scSGC_object@meta.data$customclassif
cell_state_labels <- cell_type_labels
################################################################

################################################################
# source("Pipe_Combine.bulk.data.R")         # TMP or normalized count ?
gene_cov <- get(load("/home/jyhong906/Public/Gene_cov/gene_cov.rda"))
combine_TPM_mat <- countToTpm(Keep_SGC_matrix,
                              keyType = "SYMBOL",
                              gene_cov = gene_cov)
bulk_Expr <- combine_TPM_mat
# bulk_Expr <- SGC_mat
################################################################

refPhi_obj <- refPrepare(sc_Expr, cell_type_labels, cell_state_labels, output_class = 'refPhi')
InstaPrism.res <- InstaPrism(input_type = 'refPhi',
                             bulk_Expr = bulk_Expr,
                             refPhi = refPhi_obj)

theta <- as.matrix(InstaPrism.res@Post.ini.ct@theta)
theta <- theta[which(!rownames(theta) %in% c("Intercalated duct cell","Plasma cells")),]

theta <- apply(theta, 2, function(x) x / sum(x))

theta_df <- melt(t(theta)); colnames(theta_df) <- c("sample","cell_type","fraction")
theta_df$duct <- SGC_duct_groups$Duct
theta_df$duct <- factor(theta_df$duct, levels = c("ID","ED"))
theta_df$cell_type <- factor(theta_df$cell_type, levels = c("Fibroblast",
                                                            "Ductal/basal epithelial cell",
                                                            "Dendritic cell",
                                                            # "Plasma cells",
                                                            "Myoepithelial cell",
                                                            "Stem cells"))
# save in /home/jyhong906/Project/SGC/Data/RNA/17.Cell type decomposition/Cell type decomposition.pdf (5 x 9)
ggplot(theta_df,
       aes(x = sample,
           y = fraction,
           fill = cell_type)) +

          # add an bar plot
          geom_bar(stat='identity') +

          # Remove gray background
          theme_bw() +

          # Add an color pallate
          scale_fill_manual(values = c('#D25959', # red
                                       "#FF884B", # orange
                                       "#FFD384", # yellow
                                       "#5F8D4E", # green
                                       # "#6B778D", # sky
                                       # "#293462", # blue
                                       "#905E96"  # purple
          )) + # "#FF9494","#4A3F35"

          # axis, main title
          labs(x     = "",
               y     = "Cell type fraction") +

          # eliminates background, gridlines, and chart border
          theme(axis.title       = element_text(size = 20,
                                                face = 'plain'),
                plot.background  = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border     = element_rect(colour = "black",
                                                fill   = NA,
                                                size   = 1,),
                legend.position = "bottom",
                legend.title = element_blank(),
                legend.text = element_text(size = 12),
                axis.line = element_line(color = 'black'),
                # axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.x = element_text(size = 5,
                                           colour = "black",
                                           hjust = 1,
                                           margin = ggplot2::margin(5,0,0,0), # 1번은 x, 2번은 y, 3번은 x margin, 4번은 y margin
                                           angle = 90),
                axis.text.y = element_text(size = 15,
                                           colour = "black",
                                           hjust = -0.1,
                                           margin = ggplot2::margin(0,5,0,0)),
                axis.title.x = element_text(size = 0,
                                            face = "bold",
                                            vjust = -1,
                                            margin = ggplot2::margin(0,0,0,0)),
                axis.title.y = element_text(size = 20,
                                            face = "bold",
                                            vjust= 4,
                                            margin = ggplot2::margin(0,10,0,30))) +

          facet_grid(~ duct, scales = "free") +

          # facet #
          theme(strip.text.x = element_text(size = 20,
                                            face = "bold"),
                strip.background=element_rect(fill=NA, color = NA))

# Create a list of cell types
cell_types <- unique(theta_df$cell_type)

# Loop through each cell type and perform the Wilcoxon rank-sum test
for(ct in cell_types){

          # Subset the data to get fractions for each duct type and cell type
          ID_fractions <- subset(theta_df, duct == "ID" & cell_type == ct)[,"fraction"]
          ED_fractions <- subset(theta_df, duct == "ED" & cell_type == ct)[,"fraction"]

          # print(shapiro.test(ID_fractions))
          # print(shapiro.test(ED_fractions)) # majority lower than 0.05

          # Perform the Wilcoxon rank-sum test
          test_result <- wilcox.test(ID_fractions, ED_fractions)
          boxplot(ID_fractions, ED_fractions)

          # Print the test result
          cat(paste("Cell type:", ct, "\n"))
          print(test_result)
          cat("\n")
}













































# 
# # load("/home/jyhong906/tutorial.gbm.rdata")
# 
# # bk.dat # The sample-by-gene raw count matrix of bulk RNA-seq expression
# # sc.dat # The cell-by-gene raw count matrix of bulk RNA-seq expression
# 
# # sort(table(cell.type.labels)) # cell.type.labels <- nrow(sc.dat) 
# 
# # SGC <- ComBat_seq(as.matrix(SGC), batch=SGC_groups$Batch, group=SGC_duct_groups$Duct); SGC
# 
# # sc.stat <- plot.scRNA.outlier(
# #           input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
# #           cell.type.labels=cell.type.labels,
# #           species="hs", #currently only human(hs) and mouse(mm) annotations are supported
# #           return.raw=TRUE #return the data used for plotting.
# #           #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
# # )
# 
# cell.type.labels <- scSGC_object@meta.data$customclassif
# sc.stat <- plot.scRNA.outlier(
#           input = SGC_scCounts_filter[rownames(scSGC_object@meta.data),], #make sure the colnames are gene symbol or ENSMEBL ID
#           cell.type.labels = cell.type.labels,
#           species = "hs", #currently only human(hs) and mouse(mm) annotations are supported
#           return.raw = TRUE #return the data used for plotting.
#           #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
# ); sc.stat
# 
# 
# bk.stat <- plot.bulk.outlier(
#           bulk.input = t(counts(SGC_duct_svadds)),# make sure the colnames are gene symbol or ENSMEBL ID 
#           sc.input = SGC_scCounts_filter[rownames(scSGC_object@meta.data),], #make sure the colnames are gene symbol or ENSMEBL ID 
#           cell.type.labels = scSGC_object@meta.data$customclassif,
#           species = "hs", #currently only human(hs) and mouse(mm) annotations are supported
#           return.raw = TRUE
# ); bk.stat
# 
# plot.bulk.vs.sc(
#           sc.input = SGC_scCounts_filter[rownames(scSGC_object@meta.data),],
#           bulk.input = t(counts(SGC_dds))
# )
# 
# 
# sc.dat.filtered.pc <-  select.gene.type (SGC_scCounts_filter[rownames(scSGC_object@meta.data),],
#                                          gene.type = "protein_coding")
# 
# cell.state.labels <- make.unique(scSGC_object@meta.data$customclassif, sep="_")
# 
# myPrism <- new.prism(
#           reference = sc.dat.filtered.pc, 
#           mixture = t(counts(SGC_dds)),
#           input.type = "count.matrix", 
#           cell.type.labels = cell.type.labels, 
#           cell.state.labels = cell.state.labels,
#           key = NULL
#           # outlier.cut = 0.01,
#           # outlier.fraction = 0.1,
# ); bp.res <- run.prism(prism = myPrism,
#                        n.cores = 72)
# 
# theta <- get.fraction (bp = bp.res,
#                        which.theta = "final",
#                        state.or.type = "type")
# 
# 
# for (idx in seq(nrow(as.data.frame(theta)))) {
#           print(colnames(theta)[apply(as.data.frame(theta)[idx,], 1, function(x) {which.max(x)})])
# }
# 
# theta_df <- melt(theta); colnames(theta_df) <- c("sample", "cell_type", "fraction")
# theta_df$duct <- factor(SGC_duct_groups$Duct, levels = c("ID", "ED"))
# theta_df$cell_type <- factor(theta_df$cell_type, levels = c("Fibroblast",
#                                                             "Ductal/basal epithelial cell",
#                                                             "Dendritic cell",
#                                                             "Macrophage",
#                                                             "Myoepithelial cell"))
# save(theta_df, file = "/home/jyhong906/Project/SGC/Script/RNA/17.Cell type decomposition/theta.RData")
# load("/home/jyhong906/Project/SGC/Script/RNA/17.Cell type decomposition/theta.RData")
# 
# # Save in "/home/jyhong906/Project/SGC/Data/RNA/17.Cell type decomposition/Fig. Cell type decomposition.pdf (5 x 10)"
# ggplot(theta_df,
#        aes(x = sample,
#            y = fraction,
#            fill = cell_type)) +
#           
#           # add an bar plot 
#           geom_bar(stat='identity') +
#           
#           # Remove gray background
#           theme_bw() + 
#           
#           # Add an color pallate
#           scale_fill_manual(values = c('#D25959', # red
#                                        "#FF884B", # orange 
#                                        "#FFD384", # yellow
#                                        "#5F8D4E", # green
#                                        # "#6B778D", # sky
#                                        # "#293462", # blue
#                                        "#905E96"  # purple
#                                        )) + # "#FF9494","#4A3F35"
#           
#           # axis, main title
#           labs(x     = "",
#                y     = "Cell type fraction") +
#           
#           # eliminates background, gridlines, and chart border
#           theme(axis.title       = element_text(size = 20,
#                                                 face = 'plain'),
#                 plot.background  = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border     = element_rect(colour = "black",
#                                                 fill   = NA,
#                                                 size   = 1,),
#                 legend.position = "bottom",
#                 legend.title = element_blank(),
#                 legend.text = element_text(size = 12),
#                 axis.line = element_line(color = 'black'),
#                 # axis.text.x=element_blank(),
#                 axis.ticks.x=element_blank(),
#                 axis.text.x = element_text(size = 10,
#                                            colour = "black",
#                                            hjust = 1,
#                                            margin = ggplot2::margin(5,0,0,0), # 1번은 x, 2번은 y, 3번은 x margin, 4번은 y margin
#                                            angle = 90),
#                 axis.text.y = element_text(size = 15,
#                                            colour = "black",
#                                            hjust = -0.1,
#                                            margin = ggplot2::margin(0,5,0,0)),
#                 axis.title.x = element_text(size = 0,
#                                             face = "bold",
#                                             vjust = -1,
#                                             margin = ggplot2::margin(0,0,0,0)),
#                 axis.title.y = element_text(size = 20,
#                                             face = "bold",
#                                             vjust= 4,
#                                             margin = ggplot2::margin(0,10,0,30))) +
#           
#           facet_grid(~ duct, scales = "free") +
#           
#           # facet #
#           theme(strip.text.x = element_text(size = 20,
#                                             face = "bold"),
#                 strip.background=element_rect(fill=NA, color = NA))
# 
# 
# # Create a list of cell types
# cell_types <- unique(theta_df$cell_type)
# 
# # Loop through each cell type and perform the Wilcoxon rank-sum test
# for(ct in cell_types){
#           
#           # Subset the data to get fractions for each duct type and cell type
#           ID_fractions <- subset(theta_df, duct == "ID" & cell_type == ct)[,"fraction"]
#           ED_fractions <- subset(theta_df, duct == "ED" & cell_type == ct)[,"fraction"]
#           
#           # print(shapiro.test(ID_fractions))
#           # print(shapiro.test(ED_fractions)) # majority lower than 0.05
#           
#           # Perform the Wilcoxon rank-sum test
#           test_result <- wilcox.test(ID_fractions, ED_fractions)
#           boxplot(ID_fractions, ED_fractions)
# 
#           # Print the test result
#           cat(paste("Cell type:", ct, "\n"))
#           print(test_result)
#           cat("\n")
# }
# 
# 
# pval <- c()
# for(ct in cell_types){
#           
#           # Subset the data to get fractions for each duct type and cell type
#           ID_fractions <- subset(MECA_het_df, duct == "ID_Like" & cell_type == ct)[,"fraction"]
#           ED_fractions <- subset(MECA_het_df, duct == "ED_Like" & cell_type == ct)[,"fraction"]
#           
#           # print(shapiro.test(ID_fractions))
#           # print(shapiro.test(ED_fractions)) # majority lower than 0.05
#           
#           # Perform the Wilcoxon rank-sum test
#           test_result <- wilcox.test(ID_fractions, ED_fractions)
#           pval <- c(pval, test_result$p.value)
#           boxplot(ID_fractions, ED_fractions)
#           
#           # Print the test result
#           cat(paste("Cell type:", ct, "\n"))
#           print(test_result)
#           cat("\n")
# } # p.adjust(pval, method = "fdr") # 0.7826915775 0.0748545369 0.8076910397 0.0010389879 0.0001159318
# 
# ggplot(MECA_het_df,
#        aes(x = sample,
#            y = fraction,
#            fill = cell_type)) +
#           
#           # add an bar plot 
#           geom_bar(stat='identity') +
#           
#           # Remove gray background
#           theme_bw() + 
#           
#           # Add an color pallate
#           scale_fill_manual(values = c('#D25959', # red
#                                        "#FF884B", # orange 
#                                        "#FFD384", # yellow
#                                        "#5F8D4E", # green
#                                        # "#6B778D", # sky
#                                        # "#293462", # blue
#                                        "#905E96"  # purple
#           )) + # "#FF9494","#4A3F35"
#           
#           # axis, main title
#           labs(x     = "",
#                y     = "Cell type fraction") +
#           
#           # eliminates background, gridlines, and chart border
#           theme(axis.title       = element_text(size = 20,
#                                                 face = 'plain'),
#                 plot.background  = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border     = element_rect(colour = "black",
#                                                 fill   = NA,
#                                                 size   = 1,),
#                 legend.position = "bottom",
#                 legend.title = element_blank(),
#                 legend.text = element_text(size = 12),
#                 axis.line = element_line(color = 'black'),
#                 # axis.text.x=element_blank(),
#                 axis.ticks.x=element_blank(),
#                 axis.text.x = element_text(size = 10,
#                                            colour = "black",
#                                            hjust = 1,
#                                            margin = ggplot2::margin(5,0,0,0), # 1번은 x, 2번은 y, 3번은 x margin, 4번은 y margin
#                                            angle = 90),
#                 axis.text.y = element_text(size = 15,
#                                            colour = "black",
#                                            hjust = -0.1,
#                                            margin = ggplot2::margin(0,5,0,0)),
#                 axis.title.x = element_text(size = 0,
#                                             face = "bold",
#                                             vjust = -1,
#                                             margin = ggplot2::margin(0,0,0,0)),
#                 axis.title.y = element_text(size = 20,
#                                             face = "bold",
#                                             vjust= 4,
#                                             margin = ggplot2::margin(0,10,0,30))) +
#           
#           facet_grid(~ duct, scales = "free") +
#           
#           # facet #
#           theme(strip.text.x = element_text(size = 20,
#                                             face = "bold"),
#                 strip.background=element_rect(fill=NA, color = NA))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # MECA heterogeneity #
# MECA_het_df <- theta_df[substr(theta_df$sample,1,4) == "MECA",]
# EDlike <- c(3,4,5,7,8,9,10,11,14,17,18,20,22,24,34,36,38); IDlike <- seq(40)[!seq(40) %in%EDlike]
# 
# MECA_het_df$duct <- "ED_Like"; MECA_het_df$duct[sapply(str_split(MECA_het_df$sample,"_T_"), function(x) {as.numeric(x[2])}) %in%  IDlike] <- "ID_Like"
# 
# pval <- c()
# for(ct in cell_types){
#           
#           # Subset the data to get fractions for each duct type and cell type
#           ID_fractions <- subset(MECA_het_df, duct == "ID_Like" & cell_type == ct)[,"fraction"]
#           ED_fractions <- subset(MECA_het_df, duct == "ED_Like" & cell_type == ct)[,"fraction"]
#           
#           # print(shapiro.test(ID_fractions))
#           # print(shapiro.test(ED_fractions)) # majority lower than 0.05
#           
#           # Perform the Wilcoxon rank-sum test
#           test_result <- wilcox.test(ID_fractions, ED_fractions)
#           pval <- c(pval, test_result$p.value)
#           boxplot(ID_fractions, ED_fractions)
#           
#           # Print the test result
#           cat(paste("Cell type:", ct, "\n"))
#           print(test_result)
#           cat("\n")
# } # p.adjust(pval, method = "fdr") # 0.7826915775 0.0748545369 0.8076910397 0.0010389879 0.0001159318
# 
# ggplot(MECA_het_df,
#        aes(x = sample,
#            y = fraction,
#            fill = cell_type)) +
#           
#           # add an bar plot 
#           geom_bar(stat='identity') +
#           
#           # Remove gray background
#           theme_bw() + 
#           
#           # Add an color pallate
#           scale_fill_manual(values = c('#D25959', # red
#                                        "#FF884B", # orange 
#                                        "#FFD384", # yellow
#                                        "#5F8D4E", # green
#                                        # "#6B778D", # sky
#                                        # "#293462", # blue
#                                        "#905E96"  # purple
#           )) + # "#FF9494","#4A3F35"
#           
#           # axis, main title
#           labs(x     = "",
#                y     = "Cell type fraction") +
#           
#           # eliminates background, gridlines, and chart border
#           theme(axis.title       = element_text(size = 20,
#                                                 face = 'plain'),
#                 plot.background  = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border     = element_rect(colour = "black",
#                                                 fill   = NA,
#                                                 size   = 1,),
#                 legend.position = "bottom",
#                 legend.title = element_blank(),
#                 legend.text = element_text(size = 12),
#                 axis.line = element_line(color = 'black'),
#                 # axis.text.x=element_blank(),
#                 axis.ticks.x=element_blank(),
#                 axis.text.x = element_text(size = 10,
#                                            colour = "black",
#                                            hjust = 1,
#                                            margin = ggplot2::margin(5,0,0,0), # 1번은 x, 2번은 y, 3번은 x margin, 4번은 y margin
#                                            angle = 90),
#                 axis.text.y = element_text(size = 15,
#                                            colour = "black",
#                                            hjust = -0.1,
#                                            margin = ggplot2::margin(0,5,0,0)),
#                 axis.title.x = element_text(size = 0,
#                                             face = "bold",
#                                             vjust = -1,
#                                             margin = ggplot2::margin(0,0,0,0)),
#                 axis.title.y = element_text(size = 20,
#                                             face = "bold",
#                                             vjust= 4,
#                                             margin = ggplot2::margin(0,10,0,30))) +
#           
#           facet_grid(~ duct, scales = "free") +
#           
#           # facet #
#           theme(strip.text.x = element_text(size = 20,
#                                             face = "bold"),
#                 strip.background=element_rect(fill=NA, color = NA))