source("/home/jyhong906/Project/SGC/Script/RNA/git/_09.GSEA/")

# synLogin('ID', 'password')

# Maps ENSEMBL IDs to HUGO
# Use srcType = "ensembl_gene_id" for Ensembl IDs
# Use srcType = "entrezgene" for Entrez IDs
genes2hugo <- function( v, srcType = "ensembl_gene_id" )
{
  ## Retrieve the EMSEMBL -> HUGO mapping
  ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl" )
  ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )

  ## Make sure there was at least one mapping
  if( nrow(ID) < 1 ) top( "No IDs mapped successfully" )

  ## Drop empty duds
  j <- which( ID[,2] == "" )
  if( length(j) > 0 ) ID <- ID[-j,]
  stopifnot( all( ID[,1] %in% v ) )

  ID
}

# Load RNAseq data
synRNA <- synGet( "syn2701943", downloadLocation = "/home/jyhong906/Project/SGC/Data/RNA/12.Cancer stem cell/PCBC" )
X <- read.delim( synRNA$path ) %>%
  tibble::column_to_rownames( "tracking_id" ) %>% as.matrix

# Retrieve metadata
synMeta <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
Y <- as.data.frame(synMeta) %>%
  mutate( UID = gsub("-", ".", UID) ) %>%
  tibble::column_to_rownames( "UID" )
Y[1:4,]

Y <- as.data.frame(synMeta)
Y$UID <- gsub("-", ".", as.data.frame(synMeta)$UID)
Y <- Y %>% tibble::column_to_rownames( "UID" )

# Retrieve the labels from the metadata
y <- Y[colnames(X),]
rownames(y) <- colnames(X)

# Fix the missing labels by hand
y["SC11.014BEB.133.5.6.11",] <- "EB"
y["SC12.039ECTO.420.436.92.16",] <- "ECTO"

## Drop the splice form ID from the gene names
v <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()
rownames(X) <- v

xyy <- y$Diffname_short
names(xyy) <- rownames(y)
y <- xyy
head(y)

# Map Ensembl IDs to HUGO
V <- genes2hugo( rownames(X) )
X <- X[V[,1],]
rownames(X) <- V[,2]

# Reduce gene set to the provides list.
# fnGenes: 지우고 싶은 genes
# if(!is.null(fnGenes)){
#   vGenes <- read.delim( fnGenes, header=FALSE ) %>% as.matrix() %>% drop()
#   VE <- genes2hugo( vGenes, "entrezgene" )
#   X <- X[intersect( rownames(X), VE[,2] ),]
# }

m <- apply( X, 1, mean )
X <- X - m
j <- which( y == "SC" )
X.tr <- X[,j]
X.bk <- X[,-j]

mm <- gelnet( t(X.tr), NULL, 0, 1 )
fnSig <- "/home/jyhong906/Project/SGC/Data/RNA/12.Cancer stem cell/PCBC/PCBC-stemsig.tsv"
write.table(mm$w, file = fnSig, sep = "\t", quote = FALSE, col.names = FALSE)

## Perform leave-one-out cross-validation
# auc <- c()
# for(i in 1:ncol(X.tr)){
#   ## Train a model on non-left-out data
#   X1 <- X.tr[,-i]
#   m1 <- gelnet( t(X1), NULL, 0, 1 )
#
#   ## Score the left-out sample against the background
#   s.bk <- apply( X.bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
#   s1 <- cor( m1$w, X.tr[,i], method="sp" )
#
#   ## AUC = P( left-out sample is scored above the background )
#   auc[i] <- sum( s1 > s.bk ) / length(s.bk)
#   cat( "Current AUC: ", auc[i], "\n" )
#   cat( "Average AUC: ", mean(auc), "\n" )
# }

## fnOut - filename of the output signature
## fnGenes - [optional] filename of the list of entrez ID to consider
main.train <- function( fnOut = "pcbc-stemsig.tsv", fnGenes = NULL )
{
  ## Load RNAseq data
  synRNA <- synGet( "syn2701943", downloadLocation = "/data/PCBC" )
  X <- read.delim( synRNA@filePath ) %>%
    tibble::column_to_rownames( "tracking_id" ) %>%
    as.matrix()

  ## Retrieve metadata
  synMeta <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
  Y <- synMeta@values %>%
    mutate( UID = gsub("-", ".", UID) ) %>%
    tibble::column_to_rownames( "UID" )

  ## Retrieve the labels from the metadata
  y <- Y[colnames(X),]
  names(y) <- colnames(X)

  ## Fix the missing labels by hand
  y["SC11.014BEB.133.5.6.11"] <- "EB"
  y["SC12.039ECTO.420.436.92.16"] <- "ECTO"

  ## Drop the splice form ID from the gene names
  v <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()
  rownames(X) <- v

  ## Map Ensembl IDs to HUGO
  V <- genes2hugo( rownames(X) )
  X <- X[V[,1],]
  rownames(X) <- V[,2]

  ## Reduce the gene set to the provided list (if applicable)
  if( is.null( fnGenes ) == FALSE )
  {
    vGenes <- read.delim( fnGenes, header=FALSE ) %>% as.matrix() %>% drop()
    VE <- genes2hugo( vGenes, "entrezgene" )
    X <- X[intersect( rownames(X), VE[,2] ),]
  }

  ## Mean-center the data
  m <- apply( X, 1, mean )
  X <- X - m

  ## Identify stem cell samples
  j <- which( y == "SC" )
  X.tr <- X[,j]
  X.bk <- X[,-j]

  ## Train a one-class model
  mm <- gelnet( t(X.tr), NULL, 0, 1 )

  ## Store the signature to a file
  write.table(mm$w, file = fnOut, sep = "\t", quote = FALSE, col.names = FALSE)

  ## Perform leave-one-out cross-validation
  auc <- c()
  for( i in 1:ncol(X.tr) )
  {
    ## Train a model on non-left-out data
    X1 <- X.tr[,-i]
    m1 <- gelnet( t(X1), NULL, 0, 1 )

    ## Score the left-out sample against the background
    s.bk <- apply( X.bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
    s1 <- cor( m1$w, X.tr[,i], method="sp" )

    ## AUC = P( left-out sample is scored above the background )
    auc[i] <- sum( s1 > s.bk ) / length(s.bk)
    cat( "Current AUC: ", auc[i], "\n" )
    cat( "Average AUC: ", mean(auc), "\n" )
  }

  return(auc)
}

w <- read.delim(fnSig, header = FALSE, row.names = 1 ) %>% as.matrix() %>% drop()
w[1:10]

# s <- synGet( "syn4976369", downloadLocation = "C:/Users/jyhong906/Desktop/PanCanStem/data/" ) # Download

# Auxiliary function: Reduces HUGO|POSITION gene IDs to just HUGO
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )
# X <- read.delim( s$path, as.is=TRUE, check.names=FALSE ) %>%    ## Read the raw values
#   filter( !grepl( "\\?", gene_id ) ) %>%     ## Drop genes with no mapping to HUGO
#   mutate( gene_id = f( gene_id ) ) %>%       ## Clip gene ids to HUGO
#   filter( gene_id %in% names(w) )            ## Reduce to the signature's gene set
#
# j <- grep( "SLC35E2", X[,1] )
# if( length(j) > 1 ) X <- X[-j[-1],]
# rownames(X) <- NULL
# X <- X %>% tibble::column_to_rownames( "gene_id" ) %>% as.matrix()
# stopifnot( all( rownames(X) %in% names(w) ) )
# w <- w[ rownames(X) ]
# s <- apply( X, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
# s <- s - min(s)
# s <- s / max(s)
# fnOut <- "C:/Users/jyhong906/Desktop/PanCanStem/mRNA_StemScore.tsv"
# write.table( cbind(s), file=fnOut, sep="\t", quote=FALSE, col.names=FALSE )
X <- read.delim( '/home/jyhong906/Project/SGC/Data/RNA/02.DESeq/matrix.txt', as.is=TRUE, check.names=FALSE )
X$gene_id <- rownames(X)
X <- X %>%    ## Read the raw values
  filter( !grepl( "\\?", gene_id ) ) %>%     ## Drop genes with no mapping to HUGO
  mutate( gene_id = f( gene_id ) ) %>%       ## Clip gene ids to HUGO
  filter( gene_id %in% names(w) )            ## Reduce to the signature's gene set

j <- grep( "SLC35E2", X[,1] )
if( length(j) > 1 ) X <- X[-j[-1],]
rownames(X) <- NULL
X <- X %>% tibble::column_to_rownames( "gene_id" ) %>% as.matrix()
stopifnot( all( rownames(X) %in% names(w) ) )
w <- w[ rownames(X) ]
s <- apply( X, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
s <- s - min(s)
s <- s / max(s)
fnOut <- "/home/jyhong906/Project/SGC/Data/RNA/12.Cancer stem cell/SGC_mRNA_StemScore.tsv"
# write.table( cbind(s), file=fnOut, sep="\t", quote=FALSE, col.names=FALSE )


CS_df <- read.table('/home/jyhong906/Project/SGC/Data/RNA/12.Cancer stem cell/SGC_mRNA_StemScore.tsv', header = F)
colnames(CS_df) <- c('Samples', 'Stemness_Index')
CS_df$Duct <- as.character(SGC_duct_groups$Duct)
rownames(CS_df) <- CS_df$Samples; CS_df <- CS_df[,-1]

CS_df$Duct <- factor(CS_df$Duct, levels=c("ID","ED")); CS_df$facet <- "CS"

# Save in /home/jyhong906/Project/SGC/Data/RNA/12.Cancer stem cell/Fig. Cancer stemness.pdf (5 x 4.5)
CS_p <- ggplot(CS_df,
               aes(x = Duct,
                   y = Stemness_Index,
                   fill = Duct)) +
  
  # dot
  geom_boxplot(width = 0.5,
               alpha = 0.5) +
  
  geom_dotplot(method   = "dotdensity",
               alpha    = 1,
               binwidth = (max(CS_df$Stemness_Index) - min(CS_df$Stemness_Index))/50,
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

