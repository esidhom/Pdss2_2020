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
source("./Analysis_Functions/DEanalysis.R")
source("./Analysis_Functions/Lipidomicsanalysis.R")

#Read in Data
Metabolomics.Cell <- read.csv(file = "./2017Nov_Metabolomics_shRNA_Cells.csv",stringsAsFactors = F)
Metabolomics.Cell[,19:length(names(Metabolomics.Cell))] <- NULL
Metabolomics.Media <- read.csv(file = "./2017Nov_Metabolomics_shRNA_Media.csv",stringsAsFactors = F)

#Create classes for data
Metabolomics.Final <- setClass(Class = "Metabolomics.Final", slots=list(Raw="data.frame",Imputed="data.frame",Normalized="data.frame",
                                                                        PCA="ANY",DE="ANY",Lipidomics="ANY")) 
de.results <- setClass(Class = "de.results", slots=list(pvalue="matrix",fdr="matrix",log2fc="matrix"))
lipidomics.results <- setClass(Class = "lipidomics.results", slots=list(Sum="matrix",Log2FC="matrix",FC="matrix",
                                                                        TGSUM="ANY",TGLog2="ANY",PLSUM="ANY",PLLog2="ANY",
                                                                        TGSATSUM="ANY",TGSATLog2="ANY",PLSATSUM="ANY",PLSATLog2="ANY"))

###### Cell Extract Analysis #####
Sample.Groups <- rbind(c(1:4),c(5:8),c(9:12))
Group.Names <- c("Scr","P1","P2")
Methods <- unique(Metabolomics.Cell$Method)
Methods <- Methods[which(Methods != '')]
Samples <- names(Metabolomics.Cell)[7:18]
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
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_Pdss2shRNACell_DataQuality.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
DataQualityPlots(Imputed.Cell,Normalized.Cell,Methods,MethodList,Samples,Sample.No,Sample.Groups,Group.Names)
dev.off()

#pca analysis
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_Pdss2shRNACell_PCA.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
PCA.cell <- PCAanalysis(Normalized.Cell,3,4,Samples,Sample.No,Sample.Groups,Group.Names)
dev.off()

#DE analysis
##filter data for named metabolites only
Filtered.Cell <- Normalized.Cell[which(Normalized.Cell$Metabolite != ''),]

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_Pdss2shRNACell_DEanalysis.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
DE.Cell <- DEanalysis(Filtered.Cell,Sample.No,Sample.Groups,Group.Names,0.1)
dev.off()

write.csv(x = DE.Cell@fdr, file = "./GeneLists/shRNA_Cell_DE_fdr.csv")
write.csv(x = DE.Cell@log2fc, file = "./GeneLists/shRNA_Cell_DE_log2fc.csv")

#Lipidomics Analysis
Lipidomics.Cell <- Filtered.Cell[which(Filtered.Cell$Method=="CP"),]
Lipid.classes <- c("CE","LPC","LPE","SM","Ceramide","MAG","DAG","TAG","PC","PE","PS","PI","plasmalogen")

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_Pdss2shRNACell_Lipidomicsanalysis.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
Lipidomics.CellResults <- Lipidomicsanalysis(Lipidomics.Cell,Lipid.classes,Samples,Sample.No,Sample.Groups,Group.Names)
dev.off()

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_Pdss2shRNACell_Lipidomicsanalysis_heatmap.pdf", 
                  sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE, height = 3, width=10)
p <- pheatmap(mat =Lipidomics.CellResults@Log2FC[sort(rownames(Lipidomics.CellResults@Log2FC)),], 
              scale = "row", fontsize = 6, cluster_rows = F,
         colorRampPalette(c("navy", "white", "firebrick3"))(100), border_color=NA, 
         treeheight_col=5, silent = T)
plot.new()
print(p)
p <- pheatmap(mat =Lipidomics.CellResults@TGSATLog2, scale = "row", fontsize = 6, cluster_rows = F, cluster_cols = F,
         colorRampPalette(c("navy", "white", "firebrick3"))(100), border_color=NA, 
         treeheight_col=5, silent = T)
plot.new()
print(p)
dev.off()


#print statistics for L2FC of lipidomics analysis (one-way anova + Tukey multiple pairwise-comparisons)
sink(file="Lipidomics_Statistics.txt")
temp.data <- Lipidomics.CellResults@Log2FC
melt.temp.data <- melt(temp.data)
levels(melt.temp.data$Var2) <- list(Scr = "Scr.1", Scr="Scr.2",Scr="Scr.3",Scr="Scr.4",
                                    P1 = "P1.1", P1 = "P1.2", P1 = "P1.3", P1 = "P1.4",
                                    P2 = "P2.1", P2 = "P2.2", P2 = "P2.3", P2 = "P2.4")
lapply(levels(melt.temp.data$Var1), function(X) {
  print(X)
  temp.data <- subset(melt.temp.data, melt.temp.data$Var1 == X)
  temp.aov <- aov(value ~ Var2, data = temp.data)
  print(TukeyHSD(temp.aov))
  print("----------------")
})
sink()

#separately add one additional plot of PU TGs & PU PLs together
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_Pdss2shRNACell_PULipids.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,width=5)
PU.only <- rbind(Lipidomics.CellResults@TGSATLog2[3,],Lipidomics.CellResults@PLSATLog2[3,])
df <- as.data.frame(PU.only)
df$lipid <- factor(c("PU.TG","PU.PL"), levels = c("PU.TG","PU.PL"))
df.melt <- melt(df, id.vars = "lipid")
df.melt$variable <- gsub("[[:punct:]][0-9]","",df.melt$variable)

colors.boxplot <- c("#000000","#7894A3","#507382")
ggplot(df.melt,aes(x=factor(variable, levels=Group.Names),y=value, fill=factor(variable, levels=Group.Names))) + 
  geom_boxplot()  + labs(title = "", x="", y="Log2FC") + theme_classic () +
  facet_grid(~lipid) + scale_fill_manual(values = colors.boxplot) + guides(fill="none")
dev.off()

#Save data
Cell.Final <- Metabolomics.Final(Raw=Metabolomics.Cell,Imputed=Imputed.Cell,Normalized=Normalized.Cell,
                                 PCA=PCA.cell,DE=DE.Cell,Lipidomics=Lipidomics.CellResults)
saveRDS(object=Cell.Final,file = "./shRNA_Cell_Final.rds")

###### Media Analysis #####
Sample.Groups <- rbind(c(1:4),c(5:8),c(9:12))
Group.Names <- c("Scr","P1","P2")
Methods <- unique(Metabolomics.Media$Method)
Methods <- Methods[which(Methods != '')]
Samples <- names(Metabolomics.Media)[7:18]
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
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_Pdss2shRNAMedia_DataQuality.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
DataQualityPlots(Imputed.Media,Normalized.Media,Methods,MethodList,Samples,Sample.No,Sample.Groups,Group.Names)
dev.off()

#pca analysis
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_Pdss2shRNAMedia_PCApdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE)
PCA.Media <- PCAanalysis(Normalized.Media,1,2,Samples,Sample.No,Sample.Groups,Group.Names)
dev.off()

#DE analysis
##filter data for named metabolites only
Filtered.Media <- Normalized.Media[which(Normalized.Media$Metabolite != ''),]

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"_MetabolomicsAnalysis_Pdss2shRNAMedia_DEanalysis.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,height=2.5)
DE.Media <- DEanalysis(Filtered.Media,Sample.No,Sample.Groups,Group.Names,0.1)
dev.off()

write.csv(x = DE.Media@fdr, file = "./GeneLists/shRNA_Media_DE_fdr.csv")
write.csv(x = DE.Media@log2fc, file = "./GeneLists/shRNA_Media_DE_log2fc.csv")

#Save data
Media.Final <- Metabolomics.Final(Raw=Metabolomics.Media,Imputed=Imputed.Media,Normalized=Normalized.Media,
                                  PCA=PCA.Media,DE=DE.Media)
saveRDS(object=Media.Final,file = "./shRNA_Media_Final.rds")












