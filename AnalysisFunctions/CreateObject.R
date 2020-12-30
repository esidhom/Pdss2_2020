CreateObject <- function(DATA.10X,ProjectName,CONDITION) {
  #create Seurate object
  SeuratObject <- CreateSeuratObject(counts = DATA.10X, project = ProjectName, min.cells = 5)
  SeuratObject$stim <- CONDITION
  
  #generate QC plots and save them as PDFs
  SeuratObject[["percent.mito"]] <- PercentageFeatureSet(SeuratObject, pattern = "^mt-")
  
  FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_",ProjectName,"_","QCPlots.pdf", sep = "")
  FILEPATH <- paste("Plots/",FILENAME,sep="")
  CreateQC(FILEPATH,SeuratObject)
  
  return(SeuratObject)
}

CreateQC <- function(FILEPATH,SeuratObject){
  pdf(file=FILEPATH,useDingbats=FALSE)
  plot1 <- VlnPlot(object = SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))
  plot2 <- VlnPlot(object = SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), pt.size=0)
  print(plot1)
  print(plot2)
  
  plot1 <- FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "percent.mito")
  plot2 <- FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(CombinePlots(plots = list(plot1, plot2),legend="right"))
  
  dev.off()
  
}