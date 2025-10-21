## 20251007
## DGE (Psoriasis: Post-treatment vs. Pre-treatment)

# Loading packages
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("gplots")
library("ggplot2")
library("devtools")
library("AnnotationHub")
library("AnnotationDbi")
library("ensembldb")
library("dplyr")

directory <- "/Users/sara/Desktop/Autocore/Psoriasis/counts/"
setwd(directory)
outputPrefix <-"Psoriasis_"
meta_psoriasis<- read.csv("../Psoriasis_metadata.csv")
meta_psoriasis <- meta_psoriasis %>%
  mutate(
    condition = factor(isolate, levels =c("PreTreat", "PostTreat"))
  )
condition_by_run <- setNames(meta_psoriasis$condition, meta_psoriasis$Run)


### --- Re-run the analysis without the outlier 
sampleFiles <- list.files(pattern = "_counts.txt$")
sampleNames <- sub("_counts.txt$", "", sampleFiles)
sampleCondition <- condition_by_run[sampleNames]


table <- data.frame(sampleName = sampleNames,
                    fileName = sampleFiles,
                    condition = sampleCondition)
# Creating DESeq2 object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = table,
                                       design = ~condition)

ddsHTSeq$condition

treatments <- c("PreTreat","PostTreat")
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

ddsHTSeq$condition

# Run DESeq2 differential expression analysis
dds <- DESeq(ddsHTSeq)

View(counts(dds, normalized=T))

#### Quality Control ####
pdf("../plots/Psoriasis_Dispersion_Plot.pdf", width = 11.69, height = 8.27)  # inches
par(mar = c(10, 5, 2, 2))  # bottom, left, top, right
plotDispEsts(dds)
dev.off()

pdf("../plots/Psoriasis_Boxplot_normalized_counts.pdf", width = 11.69, height = 8.27) 
par(mar = c(10, 5, 2, 2)) 
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

# Transform counts for data visualization
vst <- varianceStabilizingTransformation(dds, blind=T)

# Extract the rlog matrix from the object and compute pairwise correlation values
vst_mat <- assay(vst)
vst_cor <- cor(vst_mat)

# Plot heatmap (SDM)
library(pheatmap)
pdf("../plots/Psoriasis_Heatmap1.pdf", width = 11.69, height = 8.27) 
pheatmap(vst_cor)
dev.off()


pdf("../plots/Psoriasis_PCA.pdf", width = 11.69, height = 8.27) 
plotPCA(vst, intgroup="condition")
dev.off()

## Convert
library(AnnotationDbi)
library(org.Hs.eg.db)

ensembl_gene_ids <- sub("\\.\\d+$", "", rownames(dds))  # remove Ensembl version suffixes
hgnc_symbols <- mapIds(org.Hs.eg.db,
                       keys      = ensembl_gene_ids,
                       column    = "SYMBOL",
                       keytype   = "ENSEMBL",
                       multiVals = "first")

mcols(dds)$hgnc_symbol <- hgnc_symbols

normalized_counts_matrix <- counts(dds, normalized = TRUE)
rownames(normalized_counts_matrix) <- ifelse(!is.na(hgnc_symbols) & hgnc_symbols != "",
                                             hgnc_symbols,
                                             rownames(dds))
View(normalized_counts_matrix)


res <- results(dds)
res$hgnc_symbol <- hgnc_symbols[rownames(res)]
summary(res)

# Log fold change shrinkage for visualization and ranking
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_PostTreat_vs_PreTreat", type="ashr")
resLFC

sig <- subset(res, padj < 0.05)
head(sig)
dim(sig)

# Order results by padj value (the most significant to the least)
resOrdered <- sig[order(sig$padj),]

resdata <- merge(as.data.frame(resOrdered), 
                 as.data.frame(counts(dds,normalized =TRUE)), 
                 by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'

write.csv(resdata, file = paste0(outputPrefix,"results-with-normalized.csv"))

pdf("../plots/MAplot.pdf", width = 11.69, height = 8.27)  # inches

par(mar = c(10, 5, 2, 2))  # bottom, left, top, right

plotMA(res, ylim=c(-15,15),main = "RNAseq experiment")

dev.off()
# MA plot of RNAseq data for shrunken Log fold change dataset

#plotMA(resLFC, ylim=c(-15,15),main = "RNAseq experiment")

##########
# Enhanced Volcano plot
library("ggrepel")
library("EnhancedVolcano")

pdf("../plots/Psoriasis_Volcano_Plot.pdf", width = 11.69, height = 8.27) 
EnhancedVolcano(res,
                lab = res$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01,
                FCcutoff = 2,
                xlim = c(-10, 10))
dev.off()
##########
# MA plot2
library("ggplot2")
library("ggrepel")

resdata=data.frame(res)
gene=res$hgnc_symbol
resdata <- cbind(gene, data.frame(res, row.names=NULL))
head(resdata)
resdata$'log2baseMean' <- log2(resdata[,2])

x = ggplot(data = resdata, aes(x = log2baseMean, y = log2FoldChange)) 
x
pdf("../plots/Psoriasis_MA_Plot2.pdf", width = 11.69, height = 8.27) 
x+geom_point(aes(colour = padj))+
  scale_color_gradient(low="red", high="green")+
  geom_text_repel(aes(label=ifelse(log2FoldChange> 5,as.character(gene),'')),hjust=0.5,vjust=0.5,size=3)+
  geom_text_repel(aes(label=ifelse(log2FoldChange<3*(-1),as.character(gene),'')),hjust=0,vjust=0,size=3)

dev.off()

# Sample Distance Matrix2
library("RColorBrewer")
library("gplots")

sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vst), vst$type, sep="")
colnames(sampleDistMatrix) <- paste(colnames(vst), vst$type, sep="")
colors <- colorRampPalette( rev(brewer.pal(8, "Blues")) )(255)
pdf("../plots/Psoriasis_Heatmap2.pdf", width = 11.69, height = 8.27) 
heatmap(sampleDistMatrix,col=colors,margin = c(8,8))
dev.off()

# Heatmap of the top 500 DEGs 
library("RColorBrewer")
library("gplots")

select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:500]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=500)
pdf("../plots/Psoriasis_Heatmap_Genes.pdf", width = 11.69, height = 8.27) 
heatmap.2(assay(vst)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="Heatmap of the top 500 DEGs in SLET vs. SLE")
dev.off()

# 1) choose genes (note: your current 'select' picks top-mean genes, not DEGs)
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:500]
mat <- assay(vst)[select, ]

# 2) remove zero-variance rows before scaling
keep <- apply(mat, 1, function(x) sd(x, na.rm = TRUE) > 0)
mat  <- mat[keep, , drop = FALSE]

# 3) compute row Z-scores explicitly
z <- t(scale(t(mat)))
z[!is.finite(z)] <- 0  # just in case

# 4) symmetric limits + matching breaks and palette
lim    <- max(abs(z), na.rm = TRUE)
breaks <- seq(-lim, lim, length.out = 501)               # = length(col)+1
cols   <- colorRampPalette(c("blue", "white", "red"))(500)

pdf("../plots/SLE_Heatmap_Genes2.pdf", width = 11.69, height = 8.27)
gplots::heatmap.2(z,
                  scale = "none",                 # already scaled
                  col   = cols, breaks = breaks,  # <-- ensures a proper, filled key
                  key   = TRUE, keysize = 1.2,    # enlarge if needed
                  symkey = FALSE,                 # symmetry handled by 'breaks'
                  key.title = "Row Z-Score", key.xlab = "",
                  density.info = "none", trace = "none",
                  cexCol = 0.6, labRow = FALSE,
                  main = "Heatmap of the top 500 genes (SLET vs SLE)")
dev.off()
keep <- apply(assay(vst)[select, ], 1, function(x) sd(x, na.rm=TRUE) > 0)
mat  <- assay(vst)[select, ][keep, ]

lim    <- 2.5  # or compute as above
breaks <- seq(-lim, lim, length.out = 501)
cols   <- colorRampPalette(c("blue","white","red"))(500)

pdf("../plots/alaki.pdf", width = 25.69, height = 21.27 )
heatmap.2(mat, scale="row", col=cols, breaks=breaks,
          key=TRUE, symkey=FALSE, density.info="none", trace="none")
dev.off()
##PCA with labels
BiocManager::install("genefilter")
BiocManager::install("grDevices",force =TRUE)
library("genefilter")
library("ggplot2")
library("grDevices")
library("ggrepel")

rld <- rlogTransformation(dds, blind=T)
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vst)[select,]))

# set condition
group <- dds$condition

scores <- data.frame(pc$x, group)
scores$group=group
scores$rn <- rownames(scores)
(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(group))))
  + geom_point(size = 4)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1")
  #+ geom_text(label=scores$rn,size = 3)
  + geom_text_repel(label=scores$rn, size = 3)
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(0.90,0.14),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))
