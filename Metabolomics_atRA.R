library(dplyr, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(reshape2, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(ggplot2, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(pheatmap, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)

setwd("/Users/esidhom/Dropbox (MIT)/Harvard_MD-PhD/GrekaLab/Experiments/Analysis_metabolomics_final")

##Load all Analysis Functions
source("./Analysis_Functions/MissingValues.R")
source("./Analysis_Functions/NormalizationData.R")
source("./Analysis_Functions/DataQualityPlots.R")
source("./Analysis_Functions/PCAanalysis.R")
source("./Analysis_Functions/DEanalysis_atRA.R")
source("./Analysis_Functions/Lipidomicsanalysis.R")
source("./Analysis_Functions/PlotSelectedMetabolites.R")

#Create classes for data
Metabolomics.Final <- setClass(Class = "Metabolomics.Final", slots=list(Raw="data.frame",Imputed="data.frame",Normalized="data.frame",
                                                                        PCA="ANY",DE="ANY",Lipidomics="ANY")) 
de.results <- setClass(Class = "de.results", slots=list(pvalue="matrix",fdr="matrix",log2fc="matrix"))
lipidomics.results <- setClass(Class = "lipidomics.results", slots=list(Sum="matrix",Log2FC="matrix",FC="matrix",
                                                                        TGSUM="ANY",TGLog2="ANY",PLSUM="ANY",PLLog2="ANY",
                                                                        TGSATSUM="ANY",TGSATLog2="ANY",PLSATSUM="ANY",PLSATLog2="ANY"))


#Read in Data
Metabolomics.Cell <- read.csv(file = "./2018Oct_Metabolomics_atRA_Cells.csv",stringsAsFactors = F)
Metabolomics.Media <- read.csv(file = "./2018Oct_Metabolomics_atRA_Media.csv",stringsAsFactors = F)

###### Cell Extract Analysis #####
Sample.Groups <- rbind(c(1:3),c(4:6))
Group.Names <- c("DMSO","atRA")
Methods <- unique(Metabolomics.Cell$Method)
Methods <- Methods[which(Methods != '')]
Samples <- names(Metabolomics.Cell)[7:12]
Sample.No <- length(Samples)
MethodList <- Metabolomics.Cell$Method

for (i in 1:Sample.No) {
  temp.sample <- Samples[i]
  if (is.integer(Metabolomics.Cell[[temp.sample]])) {
    Metabolomics.Cell[[temp.sample]] <- as.numeric(Metabolomics.Cell[[temp.sample]])
  }
}

#impute missing values
Imputed.Cell <- MissingValues(Metabolomics.Cell,Sample.No)
#remove "internal standards", log2 norm, mean-center
Normalized.Cell <- NormalizationData(Imputed.Cell,Methods,Sample.No)
#make and save data quality plots
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_atRACell_DataQuality.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
DataQualityPlots(Imputed.Cell,Normalized.Cell,Methods,MethodList,Samples,Sample.No,Sample.Groups,Group.Names)
dev.off()

#pca analysis
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_atRACell_PCA.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
PCA.cell <- PCAanalysis(Normalized.Cell,3,4,Samples,Sample.No,Sample.Groups,Group.Names)
dev.off()

#DE analysis
##filter data for named metabolites only
Filtered.Cell <- Normalized.Cell[which(Normalized.Cell$Metabolite != ''),]

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_atRACell_DEanalysis.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
DE.Cell <- DEanalysis_atRA(Filtered.Cell,Sample.No,Sample.Groups,Group.Names,0.25)
dev.off()

#Lipidomics Analysis
Lipidomics.Cell <- Filtered.Cell[which(Filtered.Cell$Method=="C8-pos"),]
Lipid.classes <- c("CE","LPC","LPE","SM","Ceramide","MAG","DAG","TAG","PC","PE","PS","PI","plasmalogen",
                   "Q9","Q10")

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_atRACell_Lipidomicsanalysis.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
Lipidomics.CellResults <- Lipidomicsanalysis(Lipidomics.Cell,Lipid.classes,Samples,Sample.No,Sample.Groups,Group.Names)
dev.off()

#Create final class with all data
Metabolomics.Final <- setClass(Class = "Metabolomics.Final", slots=list(Raw="data.frame",Imputed="data.frame",Normalized="data.frame",
                                                                        PCA="ANY",DE="ANY",Lipidomics="ANY")) 
Cell.Final <- Metabolomics.Final(Raw=Metabolomics.Cell,Imputed=Imputed.Cell,Normalized=Normalized.Cell,
                                 PCA=PCA.cell,DE=DE.Cell,Lipidomics=Lipidomics.CellResults)
saveRDS(object=Cell.Final,file = "./atRA_Cell_Final.rds")

###### Media Analysis #####
Sample.Groups <- rbind(c(1:3),c(4:6))
Group.Names <- c("DMSO","atRA")
Methods <- unique(Metabolomics.Media$Method)
Methods <- Methods[which(Methods != '')]
Samples <- names(Metabolomics.Media)[7:12]
Sample.No <- length(Samples)
MethodList <- Metabolomics.Media$Method

for (i in 1:Sample.No) {
  temp.sample <- Samples[i]
  if (is.integer(Metabolomics.Media[[temp.sample]])) {
    Metabolomics.Media[[temp.sample]] <- as.numeric(Metabolomics.Media[[temp.sample]])
  }
}

#impute missing values
Imputed.Media <- MissingValues(Metabolomics.Media,Sample.No)
#remove "internal standards", log2 norm, mean-center
Normalized.Media <- NormalizationData(Imputed.Media,Methods,Sample.No)
#make and save data quality plots
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_atRAMedia_DataQuality.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
DataQualityPlots(Imputed.Media,Normalized.Media,Methods,MethodList,Samples,Sample.No,Sample.Groups,Group.Names)
dev.off()

#pca analysis
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_atRAMedia_PCA.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
PCA.media <- PCAanalysis(Normalized.Media,2,6,Samples,Sample.No,Sample.Groups,Group.Names)
dev.off()

#DE analysis
##filter data for named metabolites only
Filtered.Media <- Normalized.Media[which(Normalized.Media$Metabolite != ''),]

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_atRAMedia_DEanalysis.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,height=4)
DE.Media <- DEanalysis_atRA(Filtered.Media,Sample.No,Sample.Groups,Group.Names,0.1)
dev.off()

write.csv(x = DE.Media@fdr, file = "./GeneLists/atRA_Media_DE_fdr.csv")
write.csv(x = DE.Media@log2fc, file = "./GeneLists/atRA_Media_DE_log2fc.csv")

#Selected metabolites Analysis
Metabolites.Selected <- c("arachidonate","adrenate","docosapentaenoate")

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_atRAMedia_SelectedPUFAs.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,height=4)
Selected.Media <- PlotSelectedMetabolites(Filtered.Media,Metabolites.Selected,noCOL=3,Samples,Sample.No,Sample.Groups,Group.Names)
dev.off()

#save final data structure
Media.Final <- Metabolomics.Final(Raw=Metabolomics.Media,Imputed=Imputed.Media,Normalized=Normalized.Media,
                                 PCA=PCA.media,DE=DE.Media,Lipidomics = Selected.Media)
saveRDS(object=Media.Final,file = "./atRA_Media_Final.rds")
