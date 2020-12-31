library(dplyr, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(reshape2, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(ggplot2, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(viridis, quietly=TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(stringr, quietly=TRUE, warn.conflicts = FALSE, verbose = FALSE)


setwd("/Users/esidhom/Dropbox (MIT)/Harvard_MD-PhD/GrekaLab/Experiments/Analysis_sc/NephroSeqData")

#### Import Nephroseq data for statistics on Braf, Raf1 and Nras expression ####
FileList <- intersect (list.files(pattern="Analysis_"), list.files(pattern=".csv"))
imported.data <- lapply(FileList, function(X) read.table(file = X, header = T, sep = ",", skip=3))
names(imported.data) <- FileList

#### Organize statistics data into single dataframe ####
imported.data <- lapply(FileList, function(X) {
  temp.data <- imported.data[[X]]
  Summary.name <- paste(temp.data$Dataset, temp.data$Analysis, sep = "//")
  Gene.name <- (strsplit(X,split="_")[[1]][2] %>% strsplit(., split=".c"))[[1]][1]
  return(cbind(temp.data, Summary.name, Gene.name))
})
names(imported.data) <- FileList
all.data <- do.call("rbind",imported.data)
temp.order <- order(as.character(all.data$Summary.name), as.character(all.data$Analysis.Synopsis), 
                    as.character(all.data$Gene.name), all.data$p.Value, all.data$Fold.Change, all.data$r.Value)
all.data <- all.data[temp.order,]
all.data$Summary.name <- factor(all.data$Summary.name, 
                                rev(levels(all.data$Summary.name)[order(as.character(levels(all.data$Summary.name)))]))

#### Visualize statistics data ####
#Disease v Control Plots
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"NephroSeqv5_finalforpaper.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,width = 11.5,height=7)

df <- all.data
df$Fold.Change[which(is.na(df$Fold.Change))] <- 0
df$r.Value[which(is.na(df$r.Value))] <- 0
df$p.Value[which(df$p.Value >= 0.05)] <- 1

df <- subset (df, df$Analysis.Synopsis == "Disease vs. Control Analysis" | 
                df$Analysis.Synopsis == "Proteinuria Analysis")
df <- subset(df, df$p.Value < 0.05)
df <- subset(df, abs(df$r.Value) > 0.5 | abs(df$Fold.Change) > 1.5)

gradient.values <- c(1.1*min(df$Fold.Change),-1.5,1.5,1.1*max(df$Fold.Change))
gradient.values <- gradient.values - gradient.values[1]
gradient.values <- gradient.values/gradient.values[4]

area.breaks <- c(0, -log10(0.05), max(-log10(df$p.Value)))
area.labels <- c("p < 0.05", as.character(round(-log10(0.05), digits = 3)), 
                 as.character(round(max(-log10(df$p.Value)), digits = 3)))

ggplot(data = df) + geom_point(aes(x=as.factor(Gene.name), y=as.factor(Summary.name), color = Fold.Change, size=-log10(p.Value))) +
  scale_color_gradientn(colours=c("navy","white","white","firebrick3"), values=gradient.values, breaks = c(-1.5, 0, 1.5)) + 
  scale_size_area() + 
  facet_grid(.~Analysis.Synopsis, labeller = label_wrap_gen(width = 15, multi_line = TRUE)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1 ))

gradient.values <- c(1.1*min(df$r.Value),-0.5,0.5,1.1*max(df$r.Value))
gradient.values <- gradient.values - gradient.values[1]
gradient.values <- gradient.values/gradient.values[4]

ggplot(data = df) + geom_point(aes(x=as.factor(Gene.name), 
                                   y=as.factor(Summary.name), 
                                   color = r.Value, size=-log10(p.Value))) +
  scale_color_gradientn(colours=c("navy","white","white","firebrick3"), values=gradient.values) +
  scale_size_area() + 
  facet_grid(.~Analysis.Synopsis, labeller = label_wrap_gen(width = 15, multi_line = TRUE)) +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1 ))
dev.off()

#### Import Raf1/proteinuria data ####
FileList_RAF1 <- list.files(path = "./RAF1exp", pattern=".csv")
imported.raf1 <- lapply(FileList_RAF1, function(X) read.table(file = paste0("./RAF1exp/",X), header = T, sep = ",", skip=3))
names(imported.raf1) <- FileList_RAF1

#### Organize data into single data frame ####
imported.raf1 <- lapply(FileList_RAF1, function(X) {
  temp.data <- imported.raf1[[X]]
  Expt.name <-  strsplit(strsplit(X,split="RAF1Exp_")[[1]][2], split=".csv")[[1]][1]
  Expt.name <- array(Expt.name, dim = c(dim(temp.data)[1],1))
  temp.data <- cbind(temp.data, Expt.name)
  colnames(temp.data) <- c("Sample.Name","Proteinuria","RAF1.Expression","Expt.name")
  return(temp.data)
})
raf1.data <- do.call("rbind",imported.raf1)

#### Make Raf1 scatter plots ####
FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"NephroSeqv5_proteinuriaSCATTER_forPaper.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,width = 8,height=4)
ggplot(raf1.data, aes(x=Proteinuria, y=RAF1.Expression)) + geom_point() + geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Expt.name, nrow=1, scales = "free") + theme_bw()
dev.off()

#### Make Hodgin scatter plots ####
FileList_Hodgin <- list.files(path = "./Hodginexp", pattern=".csv")
imported.hodgin <- lapply(FileList_Hodgin, function(X) read.table(file = paste0("./Hodginexp/",X), header = T, sep = ",", skip=3))
names(imported.hodgin) <- FileList_Hodgin

hodgin.data.col <- do.call("cbind",imported.hodgin)
imported.hodgin <- lapply(FileList_Hodgin, function(X) {
  temp.data <- imported.hodgin[[X]]
  Gene.name <-  strsplit(X,split="Exp_")[[1]][1]
  Gene.name <- array(Gene.name, dim = c(dim(temp.data)[1],1))
  temp.data <- cbind(temp.data, Gene.name)
  colnames(temp.data) <- c("Sample.Name","Proteinuria","Gene.Expression","Gene.Name")
  return(temp.data)
})
hodgin.data.row <- do.call("rbind",imported.hodgin)

FILENAME <- paste(format(Sys.time(), "%Y%b%d_%H%M"),"NephroSeqv5_proteinuriaSCATTER_GPX4_forPaper.pdf", sep = "")
FILEPATH <- paste("Plots/",FILENAME,sep="")
pdf(file=FILEPATH,useDingbats=FALSE,width = 8,height=4)
p <- ggplot(subset(hodgin.data.row, Gene.Name == "GPX4"), aes(x=Proteinuria, y=Gene.Expression)) + 
  geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_bw() + labs(y="GPX4.expression")

q <- ggplot(hodgin.data.col, aes(x=GPX4Exp_Hodgin.csv.Median.centered.Log2.GPX4..reporter.g4504106_3p_at..Expression.Value, 
                            y=RAF1Exp_Hodgin.csv.Median.centered.Log2.RAF1..reporter.g4506400_3p_a_at..Expression.Value)) + 
  geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
  labs(x="GPX4.expression",y="RAF1.expression")

cowplot::plot_grid(p,q,nrow = 1)
dev.off()

#print statistics
summary(lm(RAF1Exp_Hodgin.csv.Median.centered.Log2.RAF1..reporter.g4506400_3p_a_at..Expression.Value ~
             + GPX4Exp_Hodgin.csv.Median.centered.Log2.GPX4..reporter.g4504106_3p_at..Expression.Value, data=hodgin.data.col))

summary(lm(Proteinuria ~ Gene.Expression, subset(hodgin.data.row, Gene.Name == "GPX4")))


