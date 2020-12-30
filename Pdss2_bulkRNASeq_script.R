
library("DESeq2",quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
library("biomaRt",quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
library("tximport",quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
library("dtplyr",quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
library("NOISeq",quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
library("edgeR",quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
library("pheatmap",quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
library("RColorBrewer",quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
library("gplots",quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
library("ggplot2",quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)
library("NMF",quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)

### Import and Format Data

#Import data from Dec 2017
countdata <- read.table(file="/Volumes/My Passport/Sequencing/RNASeq/Podocyte_Pdss2_shRNA/rsem.genes.counts.matrix.txt", sep="\t", header=TRUE)
colnames(countdata)[1] <- "geneID"
countdata[is.na(countdata)] <- 0 # change NA values to 0
countmatrix <- countdata[,-1]
rownames(countmatrix) <- countdata[,1]
# Convert countmatrix to matrix. this maintains colnames
countmatrix <- data.matrix(countmatrix)
#rearrange columns
countmatrix <- countmatrix[,c("Scr01_1","Scr02_1","Scr03_1","P101_1","P102_1","P103_1","P201_1","P202_1","P203_1")]
# Write count data matrix to a table
write.table(countmatrix, file="/Volumes/My Passport/Sequencing/RNASeq/Podocyte_Pdss2_shRNA/Analysis/countmatrix.csv")

# Additional biological annotations
data <- useMart(biomart = "ensembl", dataset="mmusculus_gene_ensembl")
bioanno_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","percentage_gene_gc_content", "chromosome_name","start_position","end_position", "gene_biotype"), filters = "external_gene_name", values = rownames(countmatrix), mart = data)
bioanno_transcript <- getBM(attributes = c("ensembl_gene_id", "cds_length"), filters = "external_gene_name", values = rownames(countmatrix), mart = data)
mylength <- bioanno_genes[,c("ensembl_gene_id", "start_position","end_position")]
mylength$genelength <- (mylength$end_position - mylength$start_position)  
mylength <- mylength[,c("ensembl_gene_id", "genelength")]
rownames(mylength) <- mylength$ensembl_gene_id
mygc <- bioanno_genes[,c("ensembl_gene_id", "percentage_gene_gc_content")]
rownames(mygc) <- mygc$ensembl_gene_id
mybiotypes <- bioanno_genes[,c("ensembl_gene_id", "gene_biotype")]
rownames(mybiotypes) <- mybiotypes$ensembl_gene_id
mychroms <- bioanno_genes[,c("ensembl_gene_id", "chromosome_name", "start_position","end_position")]
rownames(mychroms) <- mychroms$ensembl_gene_id
mychroms <- mychroms[, c("chromosome_name", "start_position","end_position")]

### Create DESeq2 data object
countmatrix <- read.table(file="/Volumes/My Passport/Sequencing/RNASeq/Podocyte_Pdss2_shRNA/Analysis/countmatrix.csv")
countmatrix <- data.matrix(countmatrix)
sample_info <- data.frame(condition = gsub("0[0-9]_1","",colnames(countmatrix)), row.names = colnames(countmatrix))

#Generate DESeqDataSet
dir <- "/Volumes/My Passport/Sequencing/RNASeq/Podocyte_Pdss2_shRNA/RSEM_mm10_ucsc_genomestudio_genes"
samples <- read.table("/Volumes/My Passport/Sequencing/RNASeq/Podocyte_Pdss2_shRNA/Analysis/samples.txt", header=TRUE)
files <- file.path(dir, paste0(samples$run, ".genes.results"))
txi.rsem <- tximport(files, type="rsem", txOut=FALSE)
txi.rsem$length[txi.rsem$length == 0] <- 1

DESeq.ds <- DESeqDataSetFromMatrix(countData = countdata,
                       colData = sample_info,
                       design = ~ condition)

DESeq.ds <- DESeqDataSetFromTximport (txi.rsem,
                                      colData = sample_info,
                                      design = ~ condition)
#Remove genes without any counts
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0,]
DESeq.rlog <- rlog(DESeq.ds, blind = TRUE)
rlog.norm.counts <- assay(DESeq.rlog)
write.table(rlog.norm.counts,file = "/Volumes/My Passport/Sequencing/RNASeq/Podocyte_Pdss2_shRNA/Analysis/rlog_norm_counts.csv")

#investigate different library sizes
library_sizes <- colSums(counts(DESeq.ds))
barplot(library_sizes, main = "library size comparison", ylab="sum(read counts)")
DESeq.ds <- estimateSizeFactors(DESeq.ds)
counts_normalized <- counts(DESeq.ds, normalized = TRUE)
library_size_norm <- colSums(counts_normalized)
barplot(library_size_norm, main="library size normalized", ylab="sum(read counts)")
write.table(counts_normalized,file = "/Volumes/My Passport/Sequencing/RNASeq/Podocyte_Pdss2_shRNA/Analysis/norm_counts.csv")

#PCA on normalized/transformed data
P <- plotPCA(DESeq.rlog)
P <- P + theme_bw() + ggtile("Rlog transformed counts")
print(P)

pca_data <- p
eig = pca_data$sdev^2
pca_data$explained = eig/sum(eig)
pca_xlab = paste("PC2,", pca_data$explained[2]*100, "%explained")

plot(pca_data$x[,2],pca_data$x[,3],
     col = colData(DESeq.ds)[,1],
     main = "PCA of seq. depth normalized \n and rlog - transformed read counts",
     xlab = paste("PC2,", round(pca_data$explained[2]*100, digits = 2), "%explained"),
     ylab = paste("PC3,", round(pca_data$explained[3]*100, digits = 2), "%explained"))

### DGE Analysis

DESeq.ds <- DESeq(DESeq.ds)

DGE.results <- results(DESeq.ds, contrast = c("condition","Scr","P1"),alpha = 0.05)
summary(DGE.results)
dge_p1results <- cbind(DGE.results$baseMean, DGE.results$log2FoldChange, DGE.results$lfcSE, DGE.results$stat, DGE.results$pvalue, DGE.results$padj)
colnames(dge_p1results) <- c("baseMean","log2FC","log2FCSE","Stat","Pvalue","AdjPvalue")
rownames(dge_p1results) <- rownames(DGE.results)

DGE.results <- results(DESeq.ds, contrast = c("condition","Scr","P2"),alpha = 0.05)
summary(DGE.results)
dge_p2results <- cbind(DGE.results$baseMean, DGE.results$log2FoldChange, DGE.results$lfcSE, DGE.results$stat, DGE.results$pvalue, DGE.results$padj)
colnames(dge_p2results) <- c("baseMean","log2FC","log2FCSE","Stat","Pvalue","AdjPvalue")
rownames(dge_p2results) <- rownames(DGE.results)

DGE.results <- results(DESeq.ds, contrast = c("condition","P1","P2"),alpha = 0.05)
summary(DGE.results)
dge_p1p2results <- cbind(DGE.results$baseMean, DGE.results$log2FoldChange, DGE.results$lfcSE, DGE.results$stat, DGE.results$pvalue, DGE.results$padj)
colnames(dge_p1p2results) <- c("baseMean","log2FC","log2FCSE","Stat","Pvalue","AdjPvalue")
rownames(dge_p1p2results) <- rownames(DGE.results)

write.table(dge_p1results,file = "/Volumes/My Passport/Sequencing/RNASeq/Podocyte_Pdss2_shRNA/Analysis/dge_p1results.csv")
write.table(dge_p2results,file = "/Volumes/My Passport/Sequencing/RNASeq/Podocyte_Pdss2_shRNA/Analysis/dge_p2results.csv")
write.table(dge_p1p2results,file = "/Volumes/My Passport/Sequencing/RNASeq/Podocyte_Pdss2_shRNA/Analysis/dge_p1p2results.csv")

#visualization
hist(log10(dge_p1results$pvalue), main = "frequencies of log10(p-values), Scr vs. P1")
hist(log10(dge_p2results$pvalue), main = "frequencies of log10(p-values), Scr vs. P2")

DGE.results.sorted <- DGE.results [order(DGE.results$padj),]
DGEgenes <- rownames(subset(DGE.results.sorted, padj < 0.05))
hm.mat_DGEgenes <- rlog.norm.counts[DGEgenes, ]
pheatmap(hm.mat_DGEgenes, 
         Rowv = TRUE, Colv = TRUE,
         distfun = "pearson", hclustfun = "average",
         scale = "row")

