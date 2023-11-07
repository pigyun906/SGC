source("/home/jyhong906/Project/SGC/Script/RNA/git/_13.signaling molecules/pipe_signaling molecules.R")

col_list <- c("B cell",'Dendritic cell', 'Macrophage', 'Total_T', "CYT", 'TIS', 'APM', 'Type')
Scatter_df <- as.data.frame(matrix(nrow=ncol(Immune_ssGSEA), ncol = length(col_list)))# New_df
colnames(Scatter_df) <- col_list; rownames(Scatter_df) <- colnames(Immune_ssGSEA)

# Professional Antigen presenting cells #
Scatter_df$`B cell` <- scale(as.numeric(Immune_ssGSEA[grep('B.cells',rownames(Immune_ssGSEA)),]))
Scatter_df$`Dendritic cell` <- scale(as.numeric(Immune_ssGSEA[grep('aDC',rownames(Immune_ssGSEA)),]))
Scatter_df$Macrophage <- scale(as.numeric(Immune_ssGSEA[grep('Macrophages',rownames(Immune_ssGSEA)),]))
Scatter_df$Total_T <- scale(as.numeric(Immune_ssGSEA[grep('T.cells',rownames(Immune_ssGSEA)),][2,]))

Scatter_df$Type <- c(rep('ID',58),rep('ED',36))

Scatter_df <- as.data.frame(matrix(nrow=ncol(Immune_ssGSEA), ncol = nrow(Immune_ssGSEA)))# New_df
colnames(Scatter_df) <- rownames(Immune_ssGSEA); rownames(Scatter_df) <- colnames(Immune_ssGSEA)

# Professional Antigen presenting cells #
for (idx in seq(nrow(Immune_ssGSEA))) {
          print(idx)
          Scatter_df[,idx] <- scale(Immune_ssGSEA[idx,])
}

unique_cell_list <- union(union(union(union(union(rownames(mcp_counter_df), rownames(epic_df)), rownames(quantiseq_df)), rownames(xcell_df)), rownames(cibersort_abs_df)), rownames(timer_df))

ex_cell <- intersect(unique_cell_list, rownames(new_immune_scaled_mat))

tool_df_list <- tool_df_list[-5]
new_immune_df <- as.data.frame(matrix(nrow = ncol(SGC_mat), ncol = length(ex_cell)))
colnames(new_immune_df) <- ex_cell
rownames(new_immune_df) <- colnames(SGC_mat)

for (idx in seq(ex_cell)) {
          print(idx)
          print(ex_cell[idx])
          
          method_df <- c()
          for (method_idx in seq(tool_df_list)) {
                    
                    if (ex_cell[idx] %in% rownames(get(tool_df_list[method_idx]))) {
                              method_df <- rbind(method_df, get(tool_df_list[method_idx])[ex_cell[idx],])
                    }
                    
          }

          if (length(method_df) >= 1) {
                    print(method_df)
                    new_immune_df[,ex_cell[idx]] <- scale(apply(method_df, 2, mean))         
          }
}
new_immune_df <- new_immune_df[,c(which(is.na(colSums(new_immune_df)) == F))]

# Polarization of macrophage #
Immune_df <- t(scale(t(cbind(new_immune_df[,-c(3,6,7,17,38,39)], Scatter_df[,c(7,24,25,27,34)]))))
# MQ_df <- as.data.frame(cbind(TIDE_ssGSEA[,rownames(SGC_groups)]["TAM.M2",],Immune_df[,"Macrophage M1"],Immune_df[,"Macrophage"]))
MQ_df <- as.data.frame(cbind(Immune_df[,"Macrophage M2"],Immune_df[,"Macrophage M1"],Immune_df[,"Macrophage"]))

colnames(MQ_df) <- c("M2","M1","Total_M")
# MQ_df <- t(scale(t(MQ_df)))

M_list <- c("M1", "M2")
col_list <- c('variable','subtype','value',"facet")
M_df <- as.data.frame(matrix(nrow = ncol(Immune_ssGSEA) * length(M_list), ncol = length(col_list)))
colnames(M_df) <- col_list
M_df$subtype <- rep(c(rep('ID',58),rep('ED',36)), length(M_list))
M_df$subtype <- factor(M_df$subtype, levels = c("ID", "ED"))
M_df$facet <- c(rep("M", ncol(Immune_ssGSEA) * 2))


M_df$value <-  c(apply(new_immune_scaled_mat[which(substr(rownames(new_immune_scaled_mat),1,10) == "Macrophage"),][c(3,6,9),], 2, mean),
                 apply(new_immune_scaled_mat[which(substr(rownames(new_immune_scaled_mat),1,10) == "Macrophage"),][c(4,7,10),], 2, mean))
M_df$variable <- c(rep('M1',ncol(Immune_ssGSEA)),rep('M2',ncol(Immune_ssGSEA)))

# save in /home/jyhong906/Project/SGC/Data/RNA/16.APC/Fig. Macrophage.pdf (5 x 5)
M_p <- ggplot(M_df,
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
          
          scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
          
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
          facet_wrap( ~ facet, scales="free"); M_p

# M2/M1 ratio

ID_idx <- M_df$subtype == "ID"
ED_idx <- M_df$subtype == "ED"
M1_idx <- M_df$variable == "M1"
M2_idx <- M_df$variable == "M2"
ID_ratio <- median(M_df[M_df$subtype == "ID" & M_df$variable == "M2",]$value) / median(M_df[M_df$subtype == "ID" & M_df$variable == "M1",]$value)
ED_ratio <- median(M_df[M_df$subtype == "ED" & M_df$variable == "M2",]$value) / median(M_df[M_df$subtype == "ED" & M_df$variable == "M1",]$value)

# Anti-inflammatory subtypes of DCs #
DC_list <- c("iDC", "TolDC")
col_list <- c('variable','subtype','value',"facet")
DC_df <- as.data.frame(matrix(nrow = ncol(Immune_ssGSEA) * length(DC_list), ncol = length(col_list)))
colnames(DC_df) <- col_list
DC_df$subtype <- rep(c(rep('ID',58),rep('ED',36)), length(DC_list))
DC_df$subtype <- factor(DC_df$subtype, levels = c("ID", "ED"))
DC_df$facet <- c(rep("DC", ncol(Immune_ssGSEA) * 2))
DC_df$value[1:94] <-  Immune_ssGSEA["iDC",] 
DC_df$value[95:188] <-  Immune_ssGSEA["tolDC",] 
DC_df$variable[1:188] <- c(rep('iDC',ncol(Immune_ssGSEA)),rep('tolDC',ncol(Immune_ssGSEA)))

# save in /home/jyhong906/Project/SGC/Data/RNA/16.APC/Fig. Dendritic cells.pdf (5 x 5)
DC_p <- ggplot(DC_df,
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
          facet_wrap( ~ facet, scales="free"); DC_p

# Immunosuppressive factors #
IS_list <- c("IL10", "TGFB", "IL6", "VEGF", "IDO", "PGE2")
col_list <- c('variable','subtype','value',"facet")
IS_df <- as.data.frame(matrix(nrow = ncol(Immune_ssGSEA) * length(IS_list), ncol = length(col_list)))
colnames(IS_df) <- col_list
IS_df$subtype <- rep(c(rep('ID',58),rep('ED',36)), length(IS_list))
IS_df$subtype <- factor(IS_df$subtype, levels = c("ID", "ED"))
IS_df$facet <- c(rep("IS", ncol(Immune_ssGSEA) * 2))
IS_df$value[1:94] <-  scale(SGC_mat[grep('IL10',rownames(SGC_mat)),][1,])
IS_df$value[95:188] <-  scale(SGC_mat[grep('TGFB1',rownames(SGC_mat)),][1,])
IS_df$value[189:282] <-  scale(SGC_mat[grep('IL6',rownames(SGC_mat)),][1,])
IS_df$value[283:376] <-  scale(SGC_mat[grep('VEGF',rownames(SGC_mat)),][1,])
IS_df$value[377:470] <-  scale(SGC_mat[grep('IDO',rownames(SGC_mat)),][2,])
IS_df$value[471:564] <-  scale(SGC_mat[grep('PTGES2',rownames(SGC_mat)),][1,])
IS_df$variable[1:564] <- c(rep('IL10',ncol(Immune_ssGSEA)),rep('TGFB1',ncol(Immune_ssGSEA)),rep('IL6',ncol(Immune_ssGSEA)),
                                   rep('VEGF',ncol(Immune_ssGSEA)),rep('IDO',ncol(Immune_ssGSEA)),rep('PTGES2',ncol(Immune_ssGSEA)))

# save in /home/jyhong906/Project/SGC/Data/RNA/16.APC/Fig. Immunosuppressive factors.pdf (5 x 10)
IS_p <- ggplot(IS_df,
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
          facet_wrap( ~ facet, scales="free"); IS_p