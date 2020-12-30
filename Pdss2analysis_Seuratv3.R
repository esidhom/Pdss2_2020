install.packages('Seurat')
library(Seurat, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(dplyr, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(reshape2, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(ggplot2, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(pheatmap, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(biomaRt, quietly=TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(viridis, quietly=TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(ggrepel, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(colorspace, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(reticulate, quietly=TRUE, warn.conflicts = FALSE, verbose = FALSE)
reticulate::py_install(packages ='umap-learn')

setwd("/Users/esidhom/Dropbox (MIT)/Harvard_MD-PhD/GrekaLab/Experiments/Analysis_sc")
source("./AnalysisFunctions/CreateObject.R")

#### Load biomaRt for mouse <--> human gene transfer #####
human <- useMart(host="www.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

#### VlnPlot functions ####
plot.Vln <- function(plot.list){
  list.num <- length(plot.list)
  for (i in 1:(list.num-1)) {
    temp.title <- plot.list[[i]]$labels$title
    plot.list[[i]] <- plot.list[[i]] + NoLegend() + labs(y=temp.title) + theme(axis.title.x=element_blank(),
                                                                               axis.text.x=element_blank(),
                                                                               plot.title=element_blank(),
                                                                               plot.margin = unit(c(-1,0,-1,0), "cm"))
  }
  
  temp.title <- plot.list[[list.num]]$labels$title
  plot.list[[list.num]] <- plot.list[[list.num]] + NoLegend() + labs(y=temp.title) + theme(axis.title.x=element_blank(),
                                                                                           plot.title=element_blank(),
                                                                                           plot.margin = unit(c(-1,0,-1,0), "cm"))
  
  temp.plot <- cowplot::plot_grid(plotlist = plot.list,ncol = 1,align = "hv" ,axis = "lb") +
    theme(plot.margin = unit(c(1.5,0.5,1.5,0.5), "cm")) 
  print(temp.plot)
}

plot.Vln2 <- function(plot.list){
  list.num <- length(plot.list)
  for (i in 1:(list.num)) {
    temp.title <- plot.list[[i]]$labels$title
    plot.list[[i]] <- plot.list[[i]] + NoLegend() + labs(title = temp.title) + theme(axis.title.x=element_blank(),
                                                                                     axis.text.x = element_text(angle = 0, hjust = 0.5),
                                                                                     axis.title.y=element_blank(),
                                                                                     plot.margin = unit(c(-1,0,-1,0), "cm"))
  }
  
  temp.plot <- cowplot::plot_grid(plotlist = plot.list,nrow = 1,align = "hv" ,axis = "lb") +
    theme(plot.margin = unit(c(1.5,0.5,1.5,0.5), "cm")) 
  print(temp.plot)
}
#### Import Data and plot initial QC ####
#read in 10x data
print(paste("Current Time: ",format(Sys.time(), "%Y %b %d %H:%M:%S")," - Begin reading in 10x data"))
ctrl1.data <- Read10X(data.dir = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/1/outs/filtered_gene_bc_matrices/mm10-1.2.0_premrna")
ctrl2.data <- Read10X(data.dir = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/3/outs/filtered_gene_bc_matrices/mm10-1.2.0_premrna")
ctrl3.data <- Read10X(data.dir = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/5/outs/filtered_gene_bc_matrices/mm10-1.2.0_premrna")
kdkd1.data <- Read10X(data.dir = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/2/outs/filtered_gene_bc_matrices/mm10-1.2.0_premrna")
kdkd2.data <- Read10X(data.dir = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/4/outs/filtered_gene_bc_matrices/mm10-1.2.0_premrna")
kdkd3.data <- Read10X(data.dir = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/6/outs/filtered_gene_bc_matrices/mm10-1.2.0_premrna")

#convert to sparse matrix for efficiency
print(paste("Current Time: ",format(Sys.time(), "%Y %b %d %H:%M:%S")," - Begin conversion to sparse matrix"))
ctrl1.data <- as(as.matrix(ctrl1.data), "dgCMatrix")
ctrl2.data <- as(as.matrix(ctrl2.data), "dgCMatrix")
ctrl3.data <- as(as.matrix(ctrl3.data), "dgCMatrix")
kdkd1.data <- as(as.matrix(kdkd1.data), "dgCMatrix")
kdkd2.data <- as(as.matrix(kdkd2.data), "dgCMatrix")
kdkd3.data <- as(as.matrix(kdkd3.data), "dgCMatrix")

# Create and setup Seurat objects for each dataset
print(paste("Current Time: ",format(Sys.time(), "%Y %b %d %H:%M:%S")," - Begin creation of Seurat objects"))

ctrl1 <- CreateObject(ctrl1.data,"CTRL1","CTRL")
ctrl2 <- CreateObject(ctrl2.data,"CTRL2","CTRL")
ctrl3 <- CreateObject(ctrl3.data,"CTRL3","CTRL")
kdkd1 <- CreateObject(kdkd1.data,"KDKD1","KDKD")
kdkd2 <- CreateObject(kdkd2.data,"KDKD2","KDKD")
kdkd3 <- CreateObject(kdkd3.data,"KDKD3","KDKD")


#### Combine data and filter cells ####
merge.initial <- merge(ctrl1,c(ctrl2,ctrl3,kdkd1,kdkd2,kdkd3))

#Plots pre-filtering
ProjectName <- "COMBINED"
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_",ProjectName,"_","QCPlots.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
CreateQC(FILEPATH,merge.initial)

#Cell filtering
merge.filter <- subset(merge.initial, nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mito < 5)

#Plots post-filtering
ProjectName <- "COMBINED_postFILTER"
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_",ProjectName,"_","QCPlots.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
CreateQC(FILEPATH,merge.filter)


#Print QC stats to file
sink(file = "QCstats.txt")
table(Idents(object = merge.initial))
table(Idents(object = merge.filter))
QC.list <- merge.filter[[c("orig.ident","nCount_RNA","nFeature_RNA")]]
ident.list <- unique(QC.list$orig.ident)
ident.list <- sort(ident.list)
for (i in 1:length(ident.list)) {
  print(ident.list[i])
  print(colMeans(QC.list[which(QC.list$orig.ident == ident.list[i]),2:3]))
}
sink()

#define animal sex
orig.ident <- merge.filter[["orig.ident"]]
sex <- array(data = "", dim = dim(orig.ident))
sex[which(orig.ident == "CTRL1" | orig.ident == "CTRL2" | orig.ident == "KDKD1" | orig.ident == "KDKD2")] <- "MALE"
sex[which(orig.ident == "CTRL3" | orig.ident == "KDKD3")] <- "FEMALE"
merge.filter$sex <- sex

#### Perform initial data visualization with original Normalization, Scaling ####
merge.filter <- NormalizeData(merge.filter)

merge.filter <- FindVariableFeatures(merge.filter, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(merge.filter), 10)
print(top10)

## Test different resolutions @ dim30
merge.filter <- FindClusters(merge.filter, resolution = 0.2)
table(Idents(merge.filter),t(merge.filter[["orig.ident"]]))
VlnPlot(object = merge.filter, features = c("Synpo", "Nphs1","Dcdc2a","Dock10","nCount_RNA","percent.mito"), pt.size = 0)

merge.filter <- FindClusters(merge.filter, resolution = 0.15)
table(Idents(merge.filter),t(merge.filter[["orig.ident"]]))
VlnPlot(object = merge.filter, features = c("Synpo", "Nphs1","Dcdc2a","Dock10","nCount_RNA","percent.mito"), pt.size = 0)

## Test different resolutions @ dim20
merge.filter <- FindNeighbors(merge.filter, dims = 1:20)
merge.filter <- FindClusters(merge.filter, resolution = 0.2)
table(Idents(merge.filter),t(merge.filter[["orig.ident"]]))
VlnPlot(object = merge.filter, features = c("Synpo", "Nphs1","Dcdc2a","Dock10","nCount_RNA","percent.mito"), pt.size = 0)

merge.filter <- FindClusters(merge.filter, resolution = 0.3)
table(Idents(merge.filter),t(merge.filter[["orig.ident"]]))
VlnPlot(object = merge.filter, features = c("Synpo", "Nphs1","Dcdc2a","Dock10","nCount_RNA","percent.mito"), pt.size = 0)



#### Integrate datasets using new Seurat v3 - create CTRL-only/all samples ####
animal.list <- SplitObject(merge.filter, split.by = "orig.ident")
for (i in 1:length(animal.list)) {
  animal.list[[i]] <- SCTransform(animal.list[[i]], vars.to.regress = c("nFeature_RNA","percent.mito"))
}

##Save Seurat files to server to run on computing cluster
saveRDS(object = merge.initial,file = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/v3_mergeinitial.rds")
saveRDS(object = merge.filter,file = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/v3_mergefilter.rds")
saveRDS(object = animal.list,file = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/v3_objectlist.rds")

reference.list <- animal.list[c("CTRL1","CTRL2","CTRL3")]
anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
integrated.data.ctrl <- IntegrateData(anchorset=anchors,dims = 1:30)

##performed on computing cluster
anchors <- FindIntegrationAnchors(object.list = animal.list, dims = 1:30)
integrated.data.all <- IntegrateData(anchorset=anchors,dims = 1:30)

##load in data
integrated.data.all <- readRDS(file = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/v3_integrateddataall_nFeaturemito.rds")

#### Perform clustering on all integrated data ####
DefaultAssay(integrated.data.all) <- "integrated"
integrated.data.all <- ScaleData(object=integrated.data.all,vars.to.regress = c("nFeature_RNA","percent.mito"))
integrated.data.all <- RunPCA(integrated.data.all)

ProjectName <- "INTEGRATED_ALL_regressnFeature_percentmito"
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_",ProjectName,"_","PCA_plot.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file = FILEPATH)
DimPlot(integrated.data.all, group.by = "orig.ident", cols = c('#1D3842','#3B6A76','#749EAA','#5E1742','#962E40','#C9463D'))
ElbowPlot(object=integrated.data.all, ndims = 50)
dev.off()

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_",ProjectName,"_","PCA2_plot.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file = FILEPATH,height = 12)
DimHeatmap(integrated.data.all, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(integrated.data.all, dims = 16:30, cells = 500, balanced = TRUE)
#Podocyte markers at PC#14
dev.off()

##Test 10, 15, 20, 30 for PCA visualization with 0.8 resolution to check for podocyte cluster
test.genelist <- c("Vcam1","Dock10","Synpo","Cubn","Slc12a1","Slc12a3","Slc8a1","Aqp2","Slc26a4","Kit","Emcn","Lama2","Dcdc2a")
test.pcs <- c(10,15,20,30)
res <- 0.8
for (i in 1:length(test.pcs)) {
  print(paste("Numnber of PCS:",test.pcs[i]))
  
  integrated.data.all <- FindNeighbors(integrated.data.all, dims = 1:test.pcs[i])
  integrated.data.all <- FindClusters(integrated.data.all, resolution = res)
  print(table(Idents(integrated.data.all),t(integrated.data.all[["orig.ident"]])))
  plot <- DotPlot(object = integrated.data.all, features = test.genelist, assay = "RNA", cols = c("lightgrey", "navy")) + RotatedAxis()
  print(plot)
}

## Test different resolutions @ dim15
pcs <- 15
test.res <- c(0.2, 0.3, 0.4, 0.5)
integrated.data.all <- FindNeighbors(integrated.data.all, dims = 1:pcs)
for (i in 1:length(test.res)) {
  print(paste("Resolution:",test.res[i]))
  integrated.data.all <- FindClusters(integrated.data.all, resolution = test.res[i])
  print(table(Idents(integrated.data.all),t(integrated.data.all[["orig.ident"]])))
  plot <- DotPlot(object = integrated.data.all, features = test.genelist, assay = "RNA", cols = c("lightgrey", "navy")) + RotatedAxis()
  print(plot)
}

##Final parameters: 15 PCs, Resolution 0.2
clustering.names <- grep("integrated_snn", names(integrated.data.all@meta.data),value = T)
for (i in 1:length(clustering.names)) {
  integrated.data.all@meta.data[[clustering.names[i]]] <- NULL
}

pcs <- 15
res <- 0.2
integrated.data.all <- FindNeighbors(integrated.data.all, dims = 1:pcs)
integrated.data.all <- FindClusters(integrated.data.all, resolution = res)
sink('Clustering_PCs15_res0_2.txt')
table(Idents(integrated.data.all),t(integrated.data.all[["orig.ident"]]))
sink()

#### Perform visualization + cluster identifications ####

##Run UMAP visualization and create feature plots
integrated.data.all <- RunUMAP(integrated.data.all, dims = 1:pcs)
Idents(integrated.data.all) <- "integrated_snn_res.0.2"
ProjectName <- "INTEGRATED_ALL_regressnFeaure_PCs15_res0_2"
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_",ProjectName,"_","UMAP_plot.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file = FILEPATH,useDingbats=FALSE)
DimPlot(integrated.data.all, group.by = "orig.ident", cols = c('#1D3842','#3B6A76','#749EAA','#5E1742','#962E40','#C9463D'))
DimPlot(integrated.data.all, group.by = "stim", cols = c("#749EAA","#C9463D"))
DimPlot(integrated.data.all, label = TRUE) + NoLegend()
dev.off()


FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_",ProjectName,"_","UMAP_features_QC.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file = FILEPATH,width = 15,useDingbats=FALSE)
VlnPlot(object = integrated.data.all, features = 
          c("percent.mito","nFeature_RNA","nCount_RNA"), pt.size = 0)
FeaturePlot(object = integrated.data.all, features = c("percent.mito","nFeature_RNA","nCount_RNA"), 
            cols = c("light grey", "navy"),ncol=3)
dev.off()

genelist1 <- c("Dock10","Synpo","Cubn","Slc12a1")
genelist2 <- c("Slc12a3","Slc8a1","Aqp2","Slc26a4")
genelist3 <-  c("Kit","Emcn","Lama2","Dcdc2a")
genelist <- c(genelist1, genelist2, genelist3)

DefaultAssay(integrated.data.all) <- "RNA"

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_",ProjectName,"_","UMAP_features_markergenes.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file = FILEPATH,useDingbats=FALSE,width=15)
FeaturePlot(object = integrated.data.all, features = genelist1, cols = c("light grey", "navy"),min.cutoff=c('q10'),max.cutoff=c('q90'),ncol = 2)
FeaturePlot(object = integrated.data.all, features = genelist2, cols = c("light grey", "navy"),min.cutoff=c('q10'),max.cutoff=c('q90'),ncol = 2)
FeaturePlot(object = integrated.data.all, features = genelist3, cols = c("light grey", "navy"),min.cutoff=c('q10'),max.cutoff=c('q90'),ncol = 2)
DotPlot(object = integrated.data.all, features = c(genelist,"Dock10","Pax8"), assay = "RNA", cols = c("lightgrey", "navy")) + RotatedAxis()
VlnPlot(object = integrated.data.all, features = genelist, pt.size = 0, assay = "RNA", ncol=4)
dev.off()

##PT analysis
pt.all <- subset(x=integrated.data.all, idents=c("0","1","5"))
pt.all.markers <- FindAllMarkers(object = pt.all)
write.table(pt.all.markers,file = "./ForMoran/markergenelist_pt.csv",sep = ',', row.names = FALSE)

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_",ProjectName,"_","PTanalysis.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file = FILEPATH,useDingbats=FALSE,width=10)
VlnPlot(object = pt.all, features = c("Cubn","Slc5a2","Slc5a12","Slc7a7","Steap2","Ttc36","Slc22a12","Atp11a","Slc5a10"), pt.size = 0, assay = "RNA", ncol=3)
VlnPlot(object = pt.all, features = c("Cubn","Slc5a10","Atp11a","Cyp7b1","Gramd1b","Fgf1","Napsa","Aadat","Rnf24","Ghr","Mat2a"),pt.size = 0,assay = "RNA")
dev.off()

#Create proportions plot
#make proportions plot by sample
p <- DimPlot(integrated.data.all, label = TRUE)
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]]
pdata <- pdata[order(pdata$group),]
ucols <- rev(unique(pdata$colour))

Cluster.Proportions = prop.table(x = table(Idents(integrated.data.all),t(integrated.data.all[["orig.ident"]])), margin = 2)
Cluster.Proportions = Cluster.Proportions[order(nrow(Cluster.Proportions):1),]
print(Cluster.Proportions)
Cluster.Proportions.melt <- melt(Cluster.Proportions, id.vars = "row")
ggplot(Cluster.Proportions.melt, aes(x = Var2, y = value, fill = as.factor(Var1))) + 
  geom_bar(stat = "identity", colour="black", size=0.2) + labs(y = "Cluster Proportion", x = "Sample", fill = "Cluster Identity") + 
  scale_fill_manual(values=ucols) + theme_classic()

##Save Seurat object to load on cluster for marker gene analysis
saveRDS(object = integrated.data.all, file = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/v3_integrateddataall_nFeaturemito.rds")

#Perform FindMarkers - performed on computing cluster
integrated.data.all.markers <- FindAllMarkers(object = integrated.data.all)
write.table(integrated.data.all.markers,file = "markergenelist_integratedata_all_pca_dim15_res0_2.csv",sep = ',', row.names = FALSE)

integrated.data.all.markers <- read.table("./GeneLists/markergenelist_integratedata_all_pca_dim15_res0_2.csv", header=TRUE, sep=",")
top5 <- integrated.data.all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
print(top5)
DotPlot(object=integrated.data.all, assay = "RNA",features = unique(top5$gene),cols = c("light grey","navy")) + RotatedAxis()

#### Use DoubletFinder to remove doublets - did not use after regressing out percent.mito
library(devtools, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')

#### Identify doublet clusters ####
#best practices: apply DoubletFinder on individual samples 
#pK identification per sample
DefaultAssay(integrated.data.all) <- "integrated"
animal.list <- SplitObject(integrated.data.all, split.by = "orig.ident")
pk.list <- array(0,dim=c(length(animal.list),1))
pk.plot <- array(0,dim=c(33,length(animal.list)+1))
pk.plot[,1] <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
colnames(pk.plot) <- c("pK",unique(integrated.data.all$orig.ident))

for (i in 1:length(animal.list)) {
  temp.object <- animal.list[[i]]
  #pK Identification
  sweep.res.list_kidney <- DoubletFinder::paramSweep_v3(temp.object, PCs = 1:20, sct = T)
  sweep.stats_kidney <- DoubletFinder::summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- DoubletFinder::find.pK(sweep.stats_kidney)
  pk.list[i] <- bcmvn_kidney$pK[which(bcmvn_kidney$BCmetric == max(bcmvn_kidney$BCmetric))] %>% paste %>% as.numeric
  pk.plot[,i+1] <- bcmvn_kidney$BCmetric
}
pK.max <- apply(pk.plot,2,max)
pk.plot.norm <- t(apply(pk.plot,1, function(x) {x/pK.max}))
pk.plot.norm[,1] <- pk.plot[,1]
df <- as.data.frame(pk.plot.norm) %>% melt(id.vars = "pK")
ggplot(df, aes(x=pK, y=value, color=variable)) + geom_line() + labs(y='Sample Norm. BCMVN') + theme_minimal()


for (i in 1:length(animal.list)) {
  temp.object <- animal.list[[i]]
  
  #Homotypic Doublet Proportion Estimate 
  annotations <- temp.object$integrated_snn_res.0.2
  homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
  nExp_poi <- round(0.075*length(annotations))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  #Run doubletFinder with different stringencies 
  temp.object <- DoubletFinder::doubletFinder_v3(temp.object, 
                                                 PCs = 1:20, pN = 0.25, pK = pk.list[i], nExp = nExp_poi, 
                                                 reuse.pANN = FALSE,sct=T)
  temp.pANN <- grep("pANN",names(temp.object@meta.data),value = T)
  temp.object <- DoubletFinder::doubletFinder_v3(temp.object, 
                                                 PCs = 1:20, pN = 0.25, pK = pk.list[i], nExp = nExp_poi.adj, 
                                                 reuse.pANN = temp.pANN, sct=T) 
  #Save Results 
  animal.list[[i]] <- temp.object
}

#visualize doublets
merge.doubletfinder <- merge(x=animal.list[[1]], y=c(animal.list[[2]],animal.list[[3]],animal.list[[4]],animal.list[[5]],animal.list[[6]]))
doublet.metadata.names <- grep("DF.classifications",names(merge.doubletfinder@meta.data),value = T)[c(1,3,5,7,9,11)]
doublet.metadata <- FetchData(object=merge.doubletfinder, vars = doublet.metadata.names)
doublet.metadata$cells <- rownames(doublet.metadata)
doublet.metadata.melt <- melt(doublet.metadata, id.vars = "cells")
doublet.metadata.clean <- na.omit(doublet.metadata.melt)

ProjectName <- "INTEGRATED_ALL_regressnFeaure_PCs20_res0_2"
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_",ProjectName,"_","DoubletFinder_summary.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file = FILEPATH)
merge.doubletfinder <- AddMetaData(object = merge.doubletfinder, metadata = doublet.metadata.clean$value,col.name = "DF.classifications.merge")
doublet.prop <- as.data.frame(prop.table(x=table(merge.doubletfinder$integrated_snn_res.0.2,merge.doubletfinder$DF.classifications.merge), margin = 1))
ggplot(doublet.prop, aes(x=factor(x = Var1, levels = c(0:15)), y=Freq, fill=Var2)) + geom_bar(stat="identity") + 
  labs(y = "Cluster Proportion", x = "Cluster", fill = "DF Classification") + theme_minimal()
doublet.prop.bysample <- as.data.frame(prop.table(x=table(merge.doubletfinder$integrated_snn_res.0.2,merge.doubletfinder$orig.ident,merge.doubletfinder$DF.classifications.merge), margin=1))
ggplot(doublet.prop.bysample, aes(x=factor(x = Var1, levels = c(0:15)), y=Freq, fill=Var2)) + geom_bar(stat="identity") + 
  labs(y = "Cluster Proportion", x = "Cluster", fill = "DF Classification") + facet_grid(. ~ Var3) + theme_minimal()
dev.off()

#### Look at cluster correlations ####
#calc avg expr matrix of highly variable genes
avgexpr.integrated <- AverageExpression(object = integrated.data.all, assays = "RNA", features = VariableFeatures(object = integrated.data.all))
#reorder columns in numerical order
avgexpr.integrated$RNA <- avgexpr.integrated$RNA[, order(as.integer(colnames(avgexpr.integrated$RNA)))]

#calcualte spearman correlation matrix
cormatrix <- cor(avgexpr.integrated$RNA, method="spearman")
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormatrix)
cormat.melted <- melt(upper_tri, na.rm = TRUE)

#plot data
ggplot(data = cormat.melted, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "navy", high = "firebrick3", mid = "white",
                       midpoint = mean(c(min(cormat.melted$value),1)), limit = c(min(cormat.melted$value),1), space = "Lab",
                       name="Spearman\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(0, 13, by = 1)) +
  scale_y_continuous(breaks = seq(0, 13, by = 1))


##### Re-naming of clusters #####
Idents(integrated.data.all) <- integrated.data.all$integrated_snn_res.0.2
integrated.data.all <- RenameIdents(object = integrated.data.all, 
                                    "8"="Interstitial",
                                    "7"="Endothelial",
                                    "9"="CD-betaIC",
                                    "11"="CD-alphaIC",
                                    "10"="CD-PCs",
                                    "3"="CNT",
                                    "4"="DCT",
                                    "12"="TAL/DCT",
                                    "2"="TAL",
                                    "5"="PT-S3",
                                    "0"="PT-S2",
                                    "1"="PT-S1",
                                    "13"="Podocyte",
                                    "6"="Dock10+")
integrated.data.all[["celltype"]] <- Idents(object = integrated.data.all)

##### Make figures with new annotation #####
Idents(integrated.data.all) <- integrated.data.all$celltype
DefaultAssay(integrated.data.all) <- "SCT"

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_FeaturePlots_renamed.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,width = 9)

#remake t-SNE plot with ordered and labeled clusters
p <- DimPlot(object = integrated.data.all, cols = rev(c("#bfd3e6",viridis(20)[18],plasma(18)[seq(from=4, to=8, by=2)],viridis(20)[2],
                                                        "#4eb3d3",viridis(20)[seq(from=4, to=16, by=2)])), pt.size = 0.5, label = T)


print(p)
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]]
pdata <- pdata[order(pdata$group),]
ucols <- rev(unique(pdata$colour))

#make dotplots with selected features
markers.to.plot <- c("Dock10","Nphs1","Synpo","Cubn","Slc5a2","Slc5a12","Atp11a","Slc5a10",
                     "Slc12a1","Umod","Slc12a3","Slc8a1","Scnn1b","Aqp2","Slc4a9","Atp6v0d2","Slc4a1","Slc26a4",
                     "Flt1","Emcn","Lama2","Meis1","Ptprc","Cd74")
DotPlot(object = integrated.data.all, assay = "SCT",features = rev(markers.to.plot), cols = c("light grey","navy")) + RotatedAxis()
DotPlot(object = integrated.data.all, assay = "SCT",features = rev(markers.to.plot), cols = c("#749EAA","#C9463D"), split.by = "stim") + RotatedAxis()
dev.off()

Idents(integrated.data.all) <- factor(x = Idents(integrated.data.all), levels = rev(levels(Idents(integrated.data.all))))

#make combined VlnPlots with selected features
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_FeaturePlots_renamed_Vlns.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,height = 10)
VlnPlot(object = integrated.data.all, assay = "SCT",features = markers.to.plot[1:6], 
        cols = ucols, pt.size = 0, combine = F) %>% plot.Vln
VlnPlot(object = integrated.data.all, assay = "SCT",features = markers.to.plot[7:12], 
        cols = ucols, pt.size = 0, combine = F) %>% plot.Vln
VlnPlot(object = integrated.data.all, assay = "SCT",features = markers.to.plot[13:18], 
        cols = ucols, pt.size = 0, combine = F) %>% plot.Vln
VlnPlot(object = integrated.data.all, assay = "SCT",features = markers.to.plot[19:24], 
        cols = ucols, pt.size = 0, combine = F) %>% plot.Vln
dev.off()

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_FeaturePlots_renamed_Vlns_combined.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,width = 10, height = 25)
VlnPlot(object = integrated.data.all, assay = "SCT",features = markers.to.plot, 
                     cols = ucols, pt.size = 0, combine = F) %>% plot.Vln
dev.off()

#make proportions plot by sample
Idents(integrated.data.all) <- integrated.data.all$celltype
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_renamed_Proportions.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
Cluster.Proportions = prop.table(x = table(Idents(integrated.data.all), integrated.data.all$orig.ident), margin = 2)
Cluster.Proportions = Cluster.Proportions[order(nrow(Cluster.Proportions):1),]
names(ucols) <- row.names(Cluster.Proportions)
print(Cluster.Proportions)
Cluster.Proportions.melt <- melt(Cluster.Proportions, id.vars = "row")
ggplot(Cluster.Proportions.melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") + labs(y = "Cluster Proportion", x = "Sample", fill = "Cluster Identity") + 
  scale_fill_manual(values=ucols) + theme_classic()

#make proportions plot by stimulation
Cluster.Proportions = prop.table(x = table(Idents(integrated.data.all), integrated.data.all$stim), margin = 2)
Cluster.Proportions = Cluster.Proportions[order(nrow(Cluster.Proportions):1),]
names(ucols) <- row.names(Cluster.Proportions)
print(Cluster.Proportions)
Cluster.Proportions.melt <- melt(Cluster.Proportions, id.vars = "row")
ggplot(Cluster.Proportions.melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") + labs(y = "Cluster Proportion", x = "stim", fill = "Cluster Identity") + 
  scale_fill_manual(values=ucols) + theme_classic()

dev.off()

#Podo markers
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_renamed_Podomarkers.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE, width = 10, height = 5)
Idents(integrated.data.all) <- "celltype"
VlnPlot(object = integrated.data.all, features = c("Synpo","Wt1","Plce1"), cols = c("#749EAA","#C9463D"), 
        pt.size = 0, idents = "Podocyte", split.by = "stim", assay = "SCT", combine = F) %>% plot.Vln2()
dev.off()

##### DE Expression tests #####
##look at DE genes in podocyte cluster & PT clusters
#create new meta-data structure to hold both cell type and stimulation 
integrated.data.all$celltype.stim <- paste(integrated.data.all$celltype, "_", 
                                               integrated.data.all$stim, sep = "") %>% as.factor()
integrated.data.all$celltype.stim <- factor(integrated.data.all$celltype.stim, levels = c("Interstitial_KDKD",
                                                                                          "Interstitial_CTRL",
                                                                                          "Endothelial_KDKD",
                                                                                          "Endothelial_CTRL",
                                                                                          "CD-betaIC_KDKD",
                                                                                          "CD-betaIC_CTRL",
                                                                                          "CD-alphaIC_KDKD",
                                                                                          "CD-alphaIC_CTRL",
                                                                                          "CD-PCs_KDKD",
                                                                                          "CD-PCs_CTRL",
                                                                                          "CNT_KDKD",
                                                                                          "CNT_CTRL",
                                                                                          "DCT_KDKD",
                                                                                          "DCT_CTRL",
                                                                                          "TAL/DCT_KDKD",
                                                                                          "TAL/DCT_CTRL",
                                                                                          "TAL_KDKD",
                                                                                          "TAL_CTRL",
                                                                                          "PT-S3_KDKD",
                                                                                          "PT-S3_CTRL",
                                                                                          "PT-S2_KDKD",
                                                                                          "PT-S2_CTRL",
                                                                                          "PT-S1_KDKD",
                                                                                          "PT-S1_CTRL",
                                                                                          "Podocyte_KDKD",
                                                                                          "Podocyte_CTRL",
                                                                                          "Dock10+_KDKD",
                                                                                          "Dock10+_CTRL"))
Idents(integrated.data.all) <- integrated.data.all$celltype.stim
DefaultAssay(integrated.data.all) <- "SCT"

#perform CTRL vs KDKD comparisons
kdkd.response.podo <- FindMarkers(integrated.data.all, ident.1 = "Podocyte_CTRL", ident.2 = "Podocyte_KDKD",print.bar = TRUE,logfc.threshold = 0)
head(kdkd.response.podo, 15)
write.table(kdkd.response.podo,file = paste("./GeneLists/kdkdresponse_",ProjectName,"_podocyte_logfc_0.csv",sep = ""),sep = ',')

kdkd.response.pt <- FindMarkers(integrated.data.all, ident.1 = "PT-S1_CTRL", ident.2 = "PT-S1_KDKD",print.bar = TRUE)
head(kdkd.response.pt, 15)
write.table(kdkd.response.pt,file = paste("./GeneLists/kdkdresponse_",ProjectName,"_pts1_logfc_0_25.csv",sep = ""),sep = ',')

kdkd.response.pt <- FindMarkers(integrated.data.all, ident.1 = "PT-S2_CTRL", ident.2 = "PT-S2_KDKD",print.bar = TRUE)
head(kdkd.response.pt, 15)
write.table(kdkd.response.pt,file = paste("./GeneLists/kdkdresponse_",ProjectName,"_pts2_logfc_0_25.csv",sep = ""),sep = ',')

kdkd.response.pt <- FindMarkers(integrated.data.all, ident.1 = "PT-S3_CTRL", ident.2 = "PT-S3_KDKD",print.bar = TRUE)
head(kdkd.response.pt, 15)
write.table(kdkd.response.pt,file = paste("./GeneLists/kdkdresponse_",ProjectName,"_pts3_logfc_0_25.csv",sep = ""),sep = ',')

kdkd.response.dock10 <- FindMarkers(integrated.data.all, ident.1 = "Dock10+_CTRL", ident.2 = "Dock10+_KDKD",logfc.threshold=0,print.bar = TRUE)
head(kdkd.response.dock10, 15)
write.table(kdkd.response.dock10,file = paste("./GeneLists/kdkdresponse_",ProjectName,"_dock10_logfc_0.csv",sep = ""),sep = ',')

markergenes_UnknownDock10.25 <- FindMarkers(integrated.data.all, ident.1 = "Dock10+_KDKD",print.bar = TRUE)
head(markergenes_UnknownDock10.25, 15)
write.table(markergenes_UnknownDock10,file = paste("./GeneLists/markergenes_",ProjectName,"_dock10kdkd_logfc_0.csv",sep = ""),sep = ',')

#All markergenes of Dock10+ for GSEA analysis
Idents(integrated.data.all) <- integrated.data.all$celltype
markergenes_UnknownDock10 <- FindMarkers(integrated.data.all, ident.1 = "Dock10+",print.bar = TRUE,
                                         logfc.threshold = 0)
head(markergenes_UnknownDock10, 15)
write.table(markergenes_UnknownDock10,file = paste("./GeneLists/markergenes_",ProjectName,"_dock10_logfc_0.csv",sep = ""),sep = ',')

##### Dock10 Analysis ####
Idents(integrated.data.all) <- "celltype"
DefaultAssay(integrated.data.all) <- "SCT"

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_Dock10Analysis_VlnPlots.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,height=3)
VlnPlot(object = integrated.data.all, features = c("Dock10","Vcam1","Pax8"),
        cols = c("#749EAA","#C9463D"),split.by = "stim",pt.size = 0,assay = "SCT",
        idents = "Dock10+",combine = F) %>% plot.Vln2()
dev.off()

## Dot Plot of Immune/PEC Marker Genes
Idents(integrated.data.all) <- integrated.data.all$celltype

AdaptiveList <- c("Cd3g","Cd79b","Cd79a")
InnateList <- c("Ly6c1","Cxcr4","Cd86","Cd80","C1qb","C1qa","Itgax",
                "Fcgr3","Ptprc","Csf1r","Ncr1")
PECs <- c("Cldn1","Pax8","Cd44")

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_Dock10Analysis_DotPlot.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
DotPlot(object = integrated.data.all, features = rev(c(AdaptiveList,InnateList,PECs)),
        cols = c("lightgrey", "navy"), assay = "SCT") + RotatedAxis()
dev.off()

## Correlation plot
avgexpr.integrated <- AverageExpression(object = integrated.data.all, assays = "SCT", features = VariableFeatures(object = integrated.data.all, assay = "integrated"))

#calcualte spearman correlation matrix
cormatrix <- cor(avgexpr.integrated$RNA, method="spearman")
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormatrix)
cormat.melted <- melt(upper_tri, na.rm = TRUE)

#plot data
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_Dock10Analysis_CorrPlot.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE, height=4, width=5.5)
ggplot(data = cormat.melted, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "navy", high = "firebrick3", mid = "white",
                       midpoint = mean(c(min(cormat.melted$value),1)), limit = c(min(cormat.melted$value),1), space = "Lab",
                       name="Spearman\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_fixed() 
dev.off()

## Cycling cells analysis
G1S.genes <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL",
               "PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP",
               "CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1",
               "BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8")
G2M.genes <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B",
               "MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB",
               "BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","HJURP","CDCA3","HN1",
               "CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2",
               "KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3",
               "GAS2L3","CBX5","CENPA")

# Get mouse gene names based on human genes from HMDB
DefaultAssay(integrated.data.all) <- "SCT"
annotation <- getLDS(attributes = c("hgnc_symbol"),
                     filters = "hgnc_symbol",
                     values = G1S.genes,
                     mart = human,
                     attributesL = c("mgi_symbol"), martL = mouse) 

G1S.genes <- annotation$MGI.symbol

annotation <- getLDS(attributes = c("hgnc_symbol"),
                     filters = "hgnc_symbol",
                     values = G2M.genes,
                     mart = human,
                     attributesL = c("mgi_symbol"), martL = mouse) 
G2M.genes <- annotation$MGI.symbol

integrated.data.all <- CellCycleScoring(object = integrated.data.all, s.features = G1S.genes, 
                                        g2m.features = G2M.genes)

CellCycle.Data <- FetchData(object = integrated.data.all, vars = c("S.Score","G2M.Score","celltype","stim"))
CellCycle.Filtered <- CellCycle.Data[which(CellCycle.Data$S.Score>0.15 | CellCycle.Data$G2M.Score>0.15),]

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_Dock10Analysis_CellCycle.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE, width=14)

temp.plot <- ggplot(CellCycle.Data,aes(x=S.Score,y=G2M.Score, color=CellCycle.Data$celltype=="Dock10+")) + 
  geom_point(aes(alpha=0.5)) +  labs(title ="all cells - Dock10+ highlighted") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("grey","red")) + guides(color = "none",alpha="none") +
  geom_abline(intercept = 0.15, slope = 0) + geom_vline(xintercept=0.15) 
cowplot::plot_grid(temp.plot, ncol = 2)

CellCycle.Filtered.Data <- FetchData(object=integrated.data.all,cells = rownames(CellCycle.Filtered),vars = c(G1S.genes,G2M.genes), slot = "data")
p <- pheatmap(mat = t(CellCycle.Filtered.Data),scale = "column",cluster_rows = F,show_colnames = F,
             colorRampPalette(c("navy", "white", "firebrick3"))(100), border_color=NA, fontsize_row = 4, 
             treeheight_col=10, silent=T,gaps_row = length(G1S.genes))
plot.new()
print(p)

CellCycle.Filtered.Prop <- prop.table(x = table(CellCycle.Filtered$celltype, CellCycle.Filtered$stim), margin = 2)
CellCycle.Filtered.Prop.Melt <- melt(CellCycle.Filtered.Prop)
CellCycle.Filtered.Prop.Melt$Var1 <- factor(CellCycle.Filtered.Prop.Melt$Var1, levels=rev(levels(CellCycle.Filtered.Prop.Melt$Var1)))
temp.1 <- ggplot(CellCycle.Filtered.Prop.Melt, aes(x = Var2, y = value, fill=Var1)) + geom_bar(stat = "identity", color="black") +
  labs(title = "Cycling Cells",x="",y="Prop. of Cells", fill = "Cluster Identity") + 
  scale_fill_manual(values=ucols) + theme_classic()

CellCycle.Filtered.Tbl <- table(CellCycle.Filtered$celltype, CellCycle.Filtered$stim)
CellCycle.Filtered.Tbl.Melt <- melt(CellCycle.Filtered.Tbl) 
CellCycle.Filtered.Tbl.Melt$Var1 <- factor(CellCycle.Filtered.Tbl.Melt$Var1, levels=rev(levels(CellCycle.Filtered.Tbl.Melt$Var1)))
temp.2 <- ggplot(CellCycle.Filtered.Tbl.Melt, aes(x = Var2, y = value, fill=Var1)) + 
  geom_bar(stat = "identity", color="black") + labs(title = "Cycling Cells",x="",y="#cells", fill = "Cluster Identity") + 
  scale_fill_manual(values=ucols) + theme_classic()

cowplot::plot_grid(temp.1, temp.2, ncol=2)
dev.off()

## Volcano plot of KDKD versus CTRL ##
#Data to import gene list
FILENAME = "GeneLists/markergenes_INTEGRATED_ALL_regressnFeaure_PCs15_res0_2_dock10_logfc_0.csv"
FILEPATH = paste("./",FILENAME, sep = "")

#Data to export final files 
LISTNAME = "markergenes_INTEGRATED_ALL_regressnFeaure_PCs15_res0_2_dock10_logfc_0" #filename to append to saved files
NAME = "KEGG"
VolcanoList <- c("Dock10","Vcam1","Cxcl1","Cd74","Rhoa","Rock2","Actg1",
                 "Actb","Acta2","Itgav","Itgb6","Itgb8","Itga1","Itgb1",
                 "Kras","Nras","Ifkbk","Nfkb1","Nfkbia")

geneList <-read.table(FILEPATH,header = TRUE, sep = ",")
geneList$SIGN <- apply(as.array(geneList$p_val_adj), MARGIN=1, function(x) ifelse(x < 0.05, 1, 0))
geneList$SIGN[is.element(rownames(geneList),VolcanoList)] <- 2
geneList$logpadj <- -log10(geneList$p_val_adj)

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_Dock10Analysis_Volcano.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
ggplot(data=geneList,aes(x=avg_logFC, y=logpadj))+geom_point(aes(color=as.factor(SIGN))) + 
  labs(title=LISTNAME,x="LogFC (Enriched in cluster)",y="-log10(p-adj)") + 
  theme(legend.position="top", legend.justification = "center") + 
  scale_color_manual(values = c("#CBD2D6","#749EAA","#C9463D"),name="", 
                     labels=c("non-sign.","sign. (p-adj<0.05)","selected genes")) + 
  geom_text_repel(aes(label=ifelse(is.element(rownames(geneList), VolcanoList), 
                                   as.character(rownames(geneList)),'')), 
                  force=1,point.padding=0.1) + theme_classic()
dev.off()

##### PUFA Analysis #####
DefaultAssay(integrated.data.all) <- "SCT"
Idents(integrated.data.all) <- "celltype"
PodoPT <- subset(x = integrated.data.all, idents = c("Podocyte","PT-S1","PT-S2","PT-S3"))


genelist1 <- c("Pla2g4a","Pla2g15","Pla2g6","Pla2g7","Pla2g16","Pla2g12a","Dgat2")

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_PUFAAnalysis_Pla2enzymes.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)

Idents(PodoPT) <- "celltype"
DotPlot(object = PodoPT, assay = "SCT", features = rev(genelist1),split.by = "stim",
        cols = c("#749EAA","#C9463D")) + RotatedAxis()
Idents(PodoPT) <- "celltype.stim"
DotPlot(object = PodoPT, assay = "SCT", features = rev(genelist1), 
        cols = c("light grey","navy")) + RotatedAxis()


Podocyte.mergePT <- PodoPT
celltype <- Podocyte.mergePT$celltype
celltype <- gsub("-S[0-9]","",celltype)
Podocyte.mergePT$celltype <- celltype
Idents(Podocyte.mergePT) <- "celltype"
Podocyte.mergePT$celltype.stim <- paste(Idents(Podocyte.mergePT), "_", 
                                                  Podocyte.mergePT$stim, sep = "")

DotPlot(object = Podocyte.mergePT, assay = "SCT", features = rev(genelist1),split.by = "stim",
        cols = c("#749EAA","#C9463D")) + RotatedAxis()
Idents(Podocyte.mergePT) <- "celltype.stim"
DotPlot(object = Podocyte.mergePT, assay = "SCT", features = rev(genelist1), 
        cols = c("light grey","navy")) + RotatedAxis()
dev.off()

##PUFA signature
PUFA.up <- c("GSTP1","GCLC","AOX1","DLL1","PSMD12","GCLM","EPHX1","PSMD4","DTX2","EIF4G1","PSMB3",
             "PSMD1","PSMA7","PSMC5","NUMBL","PSMD7","PSMB2","PSMD2","PSMB4","PSMC6","PSMD14","CYP7B1",
             "PSMD5","PSMD3","PSME3","BPNT1")
PUFA.down <- c("PAPSS1","PLK1","ATM","CDC14A","SEC61G","CASP3","FMO1","PSMB9","CDC6","MAT1A","ANAPC5","FCER1G","CALR","SEC61B")
annotation <- getLDS(attributes = c("hgnc_symbol"),
                     filters = "hgnc_symbol",
                     values = PUFA.up,
                     mart = human,
                     attributesL = c("mgi_symbol"), martL = mouse) 

PUFA.up <- annotation$MGI.symbol
PUFA.up <- intersect(PUFA.up,rownames(integrated.data.all))

annotation <- getLDS(attributes = c("hgnc_symbol"),
                     filters = "hgnc_symbol",
                     values = PUFA.down,
                     mart = human,
                     attributesL = c("mgi_symbol"), martL = mouse) 
PUFA.down <- annotation$MGI.symbol
PUFA.down <- intersect(PUFA.down,rownames(integrated.data.all))

integrated.data.all <- AddModuleScore(object = integrated.data.all, features = list(PUFA.up), assay = "SCT", name = "PUFA.up")
integrated.data.all <- AddModuleScore(object = integrated.data.all, features = list(PUFA.down), assay = "SCT", name = "PUFA.down")

PodoPT <- subset(x = integrated.data.all, idents = c("Podocyte","PT-S1","PT-S2","PT-S3"))
Podocyte.mergePT <- PodoPT
celltype <- Podocyte.mergePT$celltype
celltype <- gsub("-S[0-9]","",celltype)
Podocyte.mergePT$celltype <- celltype
Idents(Podocyte.mergePT) <- "celltype"
Podocyte.mergePT$celltype.stim <- paste(Idents(Podocyte.mergePT), "_", 
                                        Podocyte.mergePT$stim, sep = "")


FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_PUFAAnalysis_PUFAsign.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)

#make VlnPlots
Idents(Podocyte.mergePT) <- "celltype"
VlnPlot(object = Podocyte.mergePT, features =  c("PUFA.up1","PUFA.down1"), pt.size = 0, split.by = "stim",cols = c("#749EAA","#C9463D"))

Idents(Podocyte.mergePT) <- "celltype.stim"
VlnPlot(object = Podocyte.mergePT, features =  c("PUFA.up1","PUFA.down1"), pt.size = 0, cols = rep(c("#749EAA","#C9463D"),2))

#make Heatmaps
PUFA.geneexpression <- AverageExpression(object = Podocyte.mergePT,assays = "SCT", features = PUFA.up,use.scale = F)
p <- pheatmap(mat = PUFA.geneexpression$SCT ,scale = "row", cluster_cols = F, main="PUFA.up",silent = T,
              colorRampPalette(c("navy", "white", "firebrick3"))(100))
plot.new()
print(p)

PUFA.geneexpression <- AverageExpression(object = Podocyte.mergePT,assays = "SCT", features = PUFA.down,use.scale = F)
p <- pheatmap(mat = PUFA.geneexpression$SCT,scale = "row",cluster_cols = F, main="PUFA.down",silent = T,
              colorRampPalette(c("navy", "white", "firebrick3"))(100))
plot.new()
print(p)

##PUFA Wilcoxon rank sum test
Idents(integrated.data) <- "celltype"
celltypes <- rev(levels(integrated.data.all$celltype))
PUFA.pvalue <- array(0, dim = c(length(celltypes),2))
rownames(PUFA.pvalue) <- celltypes
colnames(PUFA.pvalue) <- c("PUFAup","PUFAdown")

PUFA.meanshift <- array(0, dim = c(length(celltypes),2))
rownames(PUFA.meanshift) <- celltypes
colnames(PUFA.meanshift) <- c("PUFAup","PUFAdown")

sink("PUFA_Sign_Stats_regressnGene.txt")
for (i in 1:length(celltypes)) {
  print(celltypes[i])
  temp_data <- subset(x = integrated.data.all,idents = celltypes[i])
  Idents(temp_data) <- "stim"
  Ctrl_cells <- WhichCells(object = temp_data, idents = "CTRL")
  Kdkd_cells <- WhichCells(object = temp_data, idents = "KDKD")
  
  PUFA_ctrl <- FetchData(object = temp_data, vars = c("PUFA.up1","PUFA.down1"), cells = Ctrl_cells)
  PUFA_kdkd <- FetchData(object = temp_data, vars = c("PUFA.up1","PUFA.down1"), cells = Kdkd_cells)
  
  PUFA.pvalue[i,1] <- wilcox.test(PUFA_ctrl$PUFA.up1,PUFA_kdkd$PUFA.up1,alternative = "less")$p.value
  PUFA.meanshift[i,1] <- mean(PUFA_kdkd$PUFA.up1,na.rm = T) - mean(PUFA_ctrl$PUFA.up,na.rm = T)
  
  PUFA.pvalue[i,2] <- wilcox.test(PUFA_ctrl$PUFA.down1,PUFA_kdkd$PUFA.down1,alternative = "greater")$p.value
  PUFA.meanshift[i,2] <- mean(PUFA_kdkd$PUFA.down1,na.rm = T) - mean(PUFA_ctrl$PUFA.down1,na.rm = T)
}
print("Wilcox Rank Sum p-values per cluster")
print(PUFA.pvalue)
print("-------------------------------------")
print("Mean-shift per cluster")
print(PUFA.meanshift)
sink()

#Plot mean-shifts for all clusters
PUFA.meanshift <- data.frame(rownames(PUFA.meanshift),PUFA.meanshift)
colnames(PUFA.meanshift)[1] <- "celltype"
PUFA.meanshift$celltype <- factor(PUFA.meanshift$celltype, levels =rev(levels(integrated.data.all$celltype)))
PUFA.meanshift.melt <- melt(PUFA.meanshift,id.vars = "celltype")


ggplot(PUFA.meanshift.melt, aes(x=celltype,y=value,fill=variable)) +
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+
  scale_fill_manual(values = c("#C9463D","#749EAA")) +
  theme_minimal() + theme(legend.title=element_blank(), 
                          axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y="mean-shift [KDKD - CTRL]",x="Cluster Identity")
dev.off()

#Plot cohens-d for all clusters
cohens_d <- function(group_A_values, group_B_values) {
  n_A <- length(group_A_values)
  n_B <- length(group_B_values)
  mu_A <- mean(group_A_values)
  mu_B <- mean(group_B_values)
  var_A <- var(group_A_values)
  var_B <- var(group_B_values)
  pooled_sd <- sqrt(((n_A - 1) * var_A + (n_B - 1) * var_B) / (n_A + n_B - 2))
  cohen_d <- (mu_A - mu_B) / pooled_sd
  cohen_d
}

celltypes <- rev(levels(integrated.data.all$celltype))
Idents(integrated.data.all) <- integrated.data.all$celltype.stim
PUFA.cohensd <- lapply(celltypes,function(x){
  Ctrl_cells <- WhichCells(object = integrated.data.all, idents = paste0(x,"_CTRL"))
  Kdkd_cells <- WhichCells(object = integrated.data.all, idents = paste0(x,"_KDKD"))
  
  PUFA_ctrl <- FetchData(object = integrated.data.all, vars = c("PUFA.up1","PUFA.down1"), cells = Ctrl_cells)
  PUFA_kdkd <- FetchData(object = integrated.data.all, vars = c("PUFA.up1","PUFA.down1"), cells = Kdkd_cells)
  
  cohens_d.up <- cohens_d(PUFA_kdkd$PUFA.up1, PUFA_ctrl$PUFA.up1)
  cohens_d.down <- cohens_d(PUFA_kdkd$PUFA.down1, PUFA_ctrl$PUFA.down1)
  cd <- cbind(x,cohens_d.up,cohens_d.down) %>% as.data.frame()
  return(cd)
}) %>% do.call("rbind", args = .)

PUFA.cohensd.df <- melt(PUFA.cohensd, id.vars = "x")

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_PUFAAnalysis_cohensd.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
ggplot(PUFA.cohensd.df, aes(x=x,y=as.numeric(value),fill=variable)) +
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+
  scale_fill_manual(values = c("#C9463D","#749EAA")) +
  theme_minimal() + theme(legend.title=element_blank(), 
                          axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y="cohen's d",x="Cluster Identity")
dev.off()


##PUFA Wilcoxon rank sum test using downsampling + Monte-Carlo simulation
#Look at effect of MCS repetitions on p-value for 102 cells #
#interate through cell clusters with significant p-value for PUFAup
Idents(integrated.data.all) <- "celltype"
celltypes <- rev(levels(integrated.data.all$celltype))
celltypes <- celltypes[which(celltypes != "Podocyte")]
ncells=length(WhichCells(object = integrated.data.all, idents = "Podocyte")) #select number of repitiions for MC simulation
nrep = 1000
rm(pval.summary2)
rm(mcs_raw.data)

for (j in 1:length(celltypes)) {
  print(celltypes[j])
  temp_data <- subset(x = integrated.data.all,idents = celltypes[j])
  Idents(temp_data) <- "stim"
  temp.ctrl <- WhichCells(object = temp_data, idents = "CTRL")
  temp.kdkd <- WhichCells(object = temp_data, idents = "KDKD")
  temp.ratio <- length(temp.ctrl)/(length(temp.kdkd)+length(temp.ctrl))
  ctrl.no <- round(temp.ratio*ncells)
  kdkd.no <- ncells - ctrl.no
  
  temp.array <- array(0, dim=c(nrep,2)) #array to output p-values for Wilcox rank sum test
  colnames(temp.array) <- c("PUFA.up","PUFA.down")
  
  #repeat nrep repititions
  for (i in 1:nrep) {
    Ctrl.subset <- sample(temp.ctrl,ctrl.no)
    Kdkd.subset <- sample(temp.kdkd,kdkd.no)
    
    
    PUFA_ctrl <- FetchData(temp_data, vars = c("PUFA.up1","PUFA.down1"), cells = Ctrl.subset)
    PUFA_kdkd <- FetchData(temp_data, vars = c("PUFA.up1","PUFA.down1"), cells = Kdkd.subset)
    
    #calculate p-value from Wilcox rank sum test with alternative hypothesis that CTRL is left-shifted to KDKD
    temp.array[i,1] <- wilcox.test(PUFA_ctrl$PUFA.up1,PUFA_kdkd$PUFA.up1,alternative = "less")$p.value
    temp.array[i,2] <- wilcox.test(PUFA_ctrl$PUFA.down1,PUFA_kdkd$PUFA.down1,alternative = "greater")$p.value
  }
  
  #save rawdata
  if (!is.element("mcs_raw.data",ls())) {
    mcs_raw.data <-  data.frame(temp.array, rep(celltypes[j],1000))
    colnames(mcs_raw.data) <- c("PUFA.up","PUFA.down","celltype")
  }
  else {
    temp.df <- data.frame(temp.array, rep(celltypes[j],1000))
    colnames(temp.df) <- c("PUFA.up","PUFA.down","celltype")
    
    mcs_raw.data <- rbind(mcs_raw.data, temp.df)
  }
  
  #calculate average & sd of p-values  from repetitions & save to data-frame
  if (!is.element("pval.summary2",ls())) {
    
    pval.summary2 <- data.frame(as.array(c(as.array(1:9),as.array(1:100)*10)), 
                                apply(as.array(c(as.array(1:9),as.array(1:100)*10)),MARGIN = 1, function (x) mean(temp.array[1:x,1])),
                                apply(as.array(c(as.array(1:9),as.array(1:100)*10)),MARGIN = 1, function (x) sd(temp.array[1:x,1])),
                                rep(celltypes[j],109),
                                rep("PUFA.up",109))
    colnames(pval.summary2) <- c("rep","pval.avg","pval.sd","celltype","signature")
    
    temp.df <- data.frame(as.array(c(as.array(1:9),as.array(1:100)*10)), 
               apply(as.array(c(as.array(1:9),as.array(1:100)*10)),MARGIN = 1, function (x) mean(temp.array[1:x,2])),
               apply(as.array(c(as.array(1:9),as.array(1:100)*10)),MARGIN = 1, function (x) sd(temp.array[1:x,2])),
               rep(celltypes[j],109),
               rep("PUFA.down",109))
    colnames(temp.df) <- c("rep","pval.avg","pval.sd","celltype","signature")
    pval.summary2 <- rbind(pval.summary2, temp.df)
    
  } else {
    temp.df1 <- data.frame(as.array(c(as.array(1:9),as.array(1:100)*10)), 
                          apply(as.array(c(as.array(1:9),as.array(1:100)*10)),MARGIN = 1, function (x) mean(temp.array[1:x,1])),
                          apply(as.array(c(as.array(1:9),as.array(1:100)*10)),MARGIN = 1, function (x) sd(temp.array[1:x,1])),
                          rep(celltypes[j],109),
                          rep("PUFA.up",109))
    colnames(temp.df1) <- c("rep","pval.avg","pval.sd","celltype","signature")
    
    temp.df2 <- data.frame(as.array(c(as.array(1:9),as.array(1:100)*10)), 
                          apply(as.array(c(as.array(1:9),as.array(1:100)*10)),MARGIN = 1, function (x) mean(temp.array[1:x,2])),
                          apply(as.array(c(as.array(1:9),as.array(1:100)*10)),MARGIN = 1, function (x) sd(temp.array[1:x,2])),
                          rep(celltypes[j],109),
                          rep("PUFA.down",109))
    colnames(temp.df2) <- c("rep","pval.avg","pval.sd","celltype","signature")
    
    pval.summary2 <- rbind(pval.summary2, temp.df1, temp.df2)
  }
  
  colnames(pval.summary2) <- c("rep","pval.avg","pval.sd","celltype","signature")
}

#log transform data
pval.summary2$logpval <- -log10(pval.summary2$pval.avg)
pval.summary2$logpsd <- abs(pval.summary2$pval.sd/(pval.summary2$pval.avg*log(10)))
pval.summary2$logrep <- log10(pval.summary2$rep)

#setup dataframes for ggplot2
temp.data <- subset(x = integrated.data.all, idents = "Podocyte") #subset data for cluster
Idents(temp.data) <- "stim"
Ctrl.subset <- WhichCells(object = temp.data, idents = "CTRL")
Kdkd.subset <- WhichCells(object = temp.data, idents = "KDKD")

PUFA_ctrl <- FetchData(temp.data,vars = c("PUFA.up1","PUFA.down1"),cells = Ctrl.subset)
PUFA_kdkd <- FetchData(temp.data,vars = c("PUFA.up1","PUFA.down1"),cells = Kdkd.subset)
Podocyte.log10pvalue <- c(-log10(wilcox.test(PUFA_ctrl$PUFA.up1,PUFA_kdkd$PUFA.up1,alternative = "less")$p.value),
                     -log10(wilcox.test(PUFA_ctrl$PUFA.down1,PUFA_kdkd$PUFA.down1,alternative = "greater")$p.value))
names(Podocyte.log10pvalue) <- c("PUFA.up", "PUFA.down") 
Podocyte.log10pvalue <- as.data.frame(Podocyte.log10pvalue)
Podocyte.log10pvalue$signature <- rownames(Podocyte.log10pvalue)

mcs_raw.data.melt <- melt(mcs_raw.data, id.vars = "celltype")
mcs_raw.data.melt$logpval <- -log10(mcs_raw.data.melt$value)
names(mcs_raw.data.melt) <- c("celltype","signature","value","logpval")
                     
#plot results
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_PUFAAnalysis_PUFAsign_montecarlo.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE, width=17)

ggplot(pval.summary2, aes(x=logrep,color=celltype)) + 
  geom_line(aes(y=logpval), alpha=0.5) + 
  geom_errorbar(aes(ymax=(logpval+logpsd),ymin=(logpval-logpsd)), alpha=0.5) +
  labs(title="Monte Carlo: random samples of 102 cells", x="log10(rep#)", y="-log10(nominal p-value)") +
  facet_grid(. ~ signature) +
  geom_hline(data=Podocyte.log10pvalue, aes(yintercept = Podocyte.log10pvalue), alpha=0.8) +
  scale_color_manual(values = ucols[c(1,3:14)]) +
  theme_minimal()

pval.summary.filter <- pval.summary2[which(pval.summary2$rep==1000),]
ggplot(pval.summary.filter, aes(x=celltype,fill=celltype)) + 
  geom_bar(aes(y=logpval), stat="identity", alpha=0.8,width = 0.9) + 
  geom_errorbar(aes(ymax=(logpval+logpsd),ymin=(logpval-logpsd)), alpha=0.8,width = 0.6) +
  labs(title="Monte Carlo: random samples of 102 cells; 1000 replicates", x="Celltype", y="-log10(nominal p-value)") +
  facet_grid(. ~ signature) +
  geom_hline(data=Podocyte.log10pvalue, aes(yintercept = Podocyte.log10pvalue), alpha=0.8) +
  geom_hline(aes(yintercept = -log10(0.05)), alpha=0.8, color="grey") +
  scale_fill_manual(values = ucols[c(1,3:14)]) + guides(fill = "none") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggplot(mcs_raw.data.melt, aes(x=celltype,fill=celltype)) + 
  geom_boxplot(aes(y=logpval)) +
  labs(title="Monte Carlo: random samples of 102 cells; 1000 replicates", x="Celltype", y="-log10(nominal p-value)") +
  facet_grid(. ~ signature) +
  geom_hline(data=Podocyte.log10pvalue, aes(yintercept = Podocyte.log10pvalue), alpha=0.8) +
  geom_hline(aes(yintercept = -log10(0.05)), alpha=0.8, color="grey") +
  scale_fill_manual(values = ucols[c(1,3:14)]) +  guides(fill = "none") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

dev.off()

##### Podocyte-PT comparison #####
etc.composite.list <- c(grep("Nduf",rownames(integrated.data.all),value = T),
                        grep("Cox",rownames(integrated.data.all),value = T),
                        grep("Uqcr",rownames(integrated.data.all),value = T),
                        grep("Sdh",rownames(integrated.data.all),value = T))


Idents(Podocyte.mergePT) <- "celltype.stim"
DefaultAssay(Podocyte.mergePT) <- "SCT"
ETC.geneexpression <- AverageExpression(object = Podocyte.mergePT, assays = "SCT", features = etc.composite.list)

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_PodoPT_ETC.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,height=14)
p <- pheatmap(mat = ETC.geneexpression$SCT,scale = "row",cluster_cols = F, main="ETC genes",silent = T,
              colorRampPalette(c("navy", "white", "firebrick3"))(100))
plot.new()
print(p)
dev.off()

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_PodoPT_ETCDotPlot.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,width=18, height = 3)
DotPlot(object = Podocyte.mergePT, assay = "SCT", features = p[["tree_row"]]$labels[p[["tree_row"]]$order], 
        cols = c("light grey", "navy")) + RotatedAxis()
dev.off()

#create ETC signature score
integrated.data.all <- AddModuleScore(object = integrated.data.all, features = list(etc.composite.list), assay = "SCT", name = "ETC.geneexpression")
PodoPT <- subset(x = integrated.data.all, idents = c("Podocyte","PT-S1","PT-S2","PT-S3"))
Podocyte.mergePT <- PodoPT
celltype <- Podocyte.mergePT$celltype
celltype <- gsub("-S[0-9]","",celltype)
Podocyte.mergePT$celltype <- celltype
Idents(Podocyte.mergePT) <- "celltype"
Podocyte.mergePT$celltype.stim <- paste(Idents(Podocyte.mergePT), "_",
                                        Podocyte.mergePT$stim, sep = "")
Idents(Podocyte.mergePT) <- "celltype"
DefaultAssay(Podocyte.mergePT) <- "SCT"
VlnPlot(object = Podocyte.mergePT, features =  "ETC.geneexpression1", pt.size = 0, split.by = "stim",cols = c("#749EAA","#C9463D"))

#Braf/Mapk expression plots
Idents(Podocyte.mergePT) <- "celltype.stim"
DefaultAssay(Podocyte.mergePT) <- "SCT"
Braf.geneexpression <- AverageExpression(object = Podocyte.mergePT,assays = "SCT", features = c("Braf","Raf1","Mapk1"))
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_PodoPT_Braf.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,height=4)
p <- pheatmap(mat = Braf.geneexpression$SCT,scale = "row",cluster_cols = F, main="Podo-sp",silent = T, 
              colorRampPalette(c("navy", "white", "firebrick3"))(100))
plot.new()
print(p)

DotPlot(object = Podocyte.mergePT, assay = "SCT", features = c("Braf","Raf1","Mapk1"), cols = c("light grey", "navy"))
dev.off()

Idents(Podocyte.mergePT) <- "celltype"
VlnPlot(object = Podocyte.mergePT, features = c("Braf","Raf1","Nras"),
        cols = c("#749EAA","#C9463D"),split.by = "stim",pt.size = 0.5,assay = "RNA",
        combine = F) %>% plot.Vln()


#Gpx4 plots
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_PodoPT_Gpx4.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,height=4, width = 4.5)
Idents(Podocyte.mergePT) <- "celltype.stim"
DotPlot(object = Podocyte.mergePT, assay = "SCT", features = c("Gpx4"), cols = c("light grey", "navy"))

Idents(Podocyte.mergePT) <- "celltype"
VlnPlot(object = Podocyte.mergePT, features = c("Gpx4"),
        cols = c("#749EAA","#C9463D"),split.by = "stim",pt.size = 0,assay = "RNA",
        combine = F) 
dev.off()

#### TAL/DCT cluster ####
DefaultAssay(integrated.data.all) <- "SCT"
Idents(integrated.data.all) <- "celltype"

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),ProjectName,"_TALDCT_Dcdc2aFeaturePlots.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,width=16)
FeaturePlot(object = integrated.data.all, features = c("Slc12a1","Slc12a3","Dcdc2a"), order=T, 
            cols = c("light grey","navy"), min.cutoff = "q10", max.cutoff = "q90", ncol = 3)
FeaturePlot(object = integrated.data.all, features = c("Slc12a1","Slc12a3", "Dcdc2a"), order=T,
            cols = c("light grey","navy"), min.cutoff = "q10", max.cutoff = "q90", ncol = 3,  
            cells = WhichCells(object = integrated.data.all, idents = c("TAL","TAL/DCT","DCT")))
FeaturePlot(object = integrated.data.all, features = c("Slc12a1","Slc12a3", "Dcdc2a"), order=T, 
            cols = c("light grey","navy"), min.cutoff = "q10", max.cutoff = "q90", ncol = 3, 
            cells = WhichCells(object = integrated.data.all, idents = c("TAL/DCT")))
VlnPlot(object = integrated.data.all, features = c("Slc12a1","Dcdc2a","Slc12a3"), cols = ucols, 
        pt.size = 0, assay = "SCT", combine = F) %>% plot.Vln()
dev.off()

##### Save data ####


saveRDS(object = integrated.data.all, file = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/v3_integrateddataall_nFeaturemito_final.rds")
integrated.data.all <- readRDS(file = "/Volumes/broad_grekalab/eriene/Pdss2_nucSeq/count_mm10premRNA/v3_integrateddataall_nFeaturemito_final.rds")
