PlotSelectedMetabolites <- function(Data,Metabolites.Selected,noCOL,Samples,Sample.No,Sample.Groups,Group.Names) {
  
  if (dim(Sample.Groups)[1] == 2) {colors.boxplot <- c("#90908F","#6D2712")} 
  else if (dim(Sample.Groups)[1] == 3) {colors.boxplot <- c("#000000","#7894A3","#507382")}
  
  temp.data <- as.matrix((Data[,7:(7+Sample.No-1)]))
  temp.names <- Data$Metabolite
  Metabolite.indexes <- lapply(Metabolites.Selected, function(x) grep(x,temp.names))
  Metabolite.indexes <- Metabolite.indexes[which(lapply(Metabolite.indexes, length) > 0)]
  Metabolite.indexes <- lapply(Metabolite.indexes, function(x) if (length(x) > 1) {
    temp.indexes <- as.integer(unlist(x))
    selected.index <- grep("TF",Data$Compound[temp.indexes],value=F)
    return(temp.indexes[selected.index])} else {return(x)})
  names(Metabolite.indexes) <- Metabolites.Selected
  
  Met.data <- array(0,dim=c(length(Metabolite.indexes),Sample.No))
  for (i in 1:length(Metabolite.indexes)) {
    temp.index <- as.integer(unlist(Metabolite.indexes[i]))
    Met.data[i,] <- temp.data[temp.index,]
  }
  rownames(Met.data) <- Metabolites.Selected
  colnames(Met.data) <- Samples
  df <- as.data.frame(Met.data)
  df$met <- rownames(df)
  df.melt <- melt(df,id.vars = "met")
  df.melt$variable <- gsub("[[:punct:]][0-9]","",df.melt$variable)
  
  p <- ggplot(df.melt,aes(x=factor(variable, levels=Group.Names),y=value, fill=factor(variable, levels=Group.Names))) + 
    geom_boxplot()  + labs(title = "", x="", y="") + theme_minimal () +
    facet_wrap(~met, scales="free",ncol = noCOL) +
    scale_fill_manual(values = colors.boxplot) + guides(fill="none")
  print(p)
  
  return(Met.data)
}