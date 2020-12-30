Lipidomicsanalysis <- function(Data,Lipid.classes,Samples,Sample.No,Sample.Groups,Group.Names) {
  #find classes of lipids
  temp.data <- as.matrix((Data[,7:(7+Sample.No-1)]))
  temp.names <- Data$Metabolite
  Lipid.indexes <- lapply(Lipid.classes, function(x) grep(x,temp.names))
  names(Lipid.indexes) <- Lipid.classes
  Lipid.indexes <- Lipid.indexes[which(lapply(Lipid.indexes, length) > 0)]
  Lipid.indexes$PL <- as.integer(unlist(Lipid.indexes[c("PC","PE","PS")]))
  Lipid.classes <- names(Lipid.indexes)
  
  if (dim(Sample.Groups)[1] == 2) {colors.boxplot <- c("#90908F","#6D2712")} 
  else if (dim(Sample.Groups)[1] == 3) {colors.boxplot <- c("#000000","#7894A3","#507382")}
  
  Sum.byclass <- array(0,dim=c(length(Lipid.classes),Sample.No))
  for (i in 1:length(Lipid.classes)) {
    temp.array <- temp.data[Lipid.indexes[[i]],]
    if (length(Lipid.indexes[[i]]) > 1) {
      Sum.byclass[i,] <- log2(colSums(2^(temp.array)))
    } else {
      Sum.byclass[i,] <- temp.array
    }
  }
  rownames(Sum.byclass) <- Lipid.classes
  colnames(Sum.byclass) <- Samples
  df <- as.data.frame(Sum.byclass)
  df$class <- rownames(df)
  df.melt <- melt(df)
  df.melt$variable <- gsub("[[:punct:]][0-9]","",df.melt$variable)
  
  p <- ggplot(df.melt,aes(x=factor(variable, levels=Group.Names),y=value, fill=factor(variable, levels=Group.Names))) + 
    geom_boxplot() + facet_wrap(~class, scales = "free") + labs(title = "log2norm(sum)", x="", y="") +
    scale_fill_manual(values = colors.boxplot) + guides(fill="none") + theme_minimal()
  
  print(p)
  
  Log2fc.byclass <- array(0,dim=c(length(Lipid.classes),Sample.No))
  for (i in 1:length(Lipid.classes)) {
    temp.array <- Sum.byclass[i,]
    temp.mean <- mean(temp.array[1:dim(Sample.Groups)[2]])
    Log2fc.byclass[i,] <- temp.array - temp.mean
  }
  
  rownames(Log2fc.byclass) <- Lipid.classes
  colnames(Log2fc.byclass) <- Samples
  df <- as.data.frame(Log2fc.byclass)
  df$class <- rownames(df)
  df.melt <- melt(df)
  df.melt$variable <- gsub("[[:punct:]][0-9]","",df.melt$variable)
  
  p <- ggplot(df.melt,aes(x=factor(variable, levels=Group.Names),y=value, fill=factor(variable, levels=Group.Names))) + 
    geom_boxplot() + facet_wrap(~class, scales = "free") + labs(title = "log2FC", x="", y="") + 
    scale_fill_manual(values = colors.boxplot) + guides(fill="none") + theme_minimal()
  print(p)
  
  FC.byclass <- 2^Log2fc.byclass
  rownames(FC.byclass) <- Lipid.classes
  colnames(FC.byclass) <- Samples
  df <- as.data.frame(FC.byclass)
  df$class <- rownames(df)
  df.melt <- melt(df)
  df.melt$variable <- gsub("[[:punct:]][0-9]","",df.melt$variable)
  
  p <- ggplot(df.melt,aes(x=factor(variable, levels=Group.Names),y=value, fill=factor(variable, levels=Group.Names))) + 
    geom_boxplot() + facet_wrap(~class, scales = "free") + labs(title = "FC", x="", y="") +
    scale_fill_manual(values = colors.boxplot) + guides(fill="none") + theme_minimal()
  print(p)
  
  ##### Analyze TGs and PLs by #double bonds ####
  
  ### Analyze TGs
  tg.metabolites <- Data[Lipid.indexes$TAG,]$Metabolite
  tg.data <- temp.data[Lipid.indexes$TAG,]
  tg.db.count <- array(0,dim=c(length(tg.metabolites),1))
  for (i in 1:length(tg.db.count)) {
    temp.metab <- tg.metabolites[i]
    temp.metab <- gsub(" TAG","",temp.metab)
    temp.metab <- unlist(strsplit(temp.metab,""))
    colon.location <- grep(":",temp.metab)
    tg.db.count[i] <- as.integer(paste(temp.metab[(colon.location+1):length(temp.metab)],collapse=""))
  }
  
  unique.tg.db <- unique(tg.db.count)
  Sum.tg.bydb <- array(0, dim=c(length(unique.tg.db),Sample.No))
  for (i in 1:length(unique.tg.db)) {
    temp.array <- tg.data[tg.db.count == unique.tg.db[i],]
    no.rows <- length(temp.array)/Sample.No
    if (no.rows > 1) {
      Sum.tg.bydb[i,] <- log2(colSums(2^(temp.array)))
    } else {
      Sum.tg.bydb[i,] <- temp.array
    }
  }
  colnames(Sum.tg.bydb) <- Samples
  rownames(Sum.tg.bydb) <- paste("dbs_",unique.tg.db,sep="")
  
  Log2fc.tg.bydb <- array(0,dim=c(length(unique.tg.db),Sample.No))
  for (i in 1:length(unique.tg.db)) {
    temp.array <- Sum.tg.bydb[i,]
    temp.mean <- mean(temp.array[1:dim(Sample.Groups)[2]])
    Log2fc.tg.bydb[i,] <- temp.array - temp.mean
  }
  colnames(Log2fc.tg.bydb) <- Samples
  rownames(Log2fc.tg.bydb) <- paste("dbs_",unique.tg.db,sep="")
  
  df <- data.frame(Log2fc.tg.bydb)
  unique.tg.db <- sort(unique.tg.db)
  df$dbs <- factor(rownames(Log2fc.tg.bydb),levels=paste("dbs_",unique.tg.db,sep=""))
  df.melt <- melt(df,id.vars = "dbs")
  df.melt$variable <- gsub("[[:punct:]][0-9]","",df.melt$variable)
  
  p1 <- ggplot(df.melt,aes(x=factor(variable, levels=Group.Names),y=value, fill=factor(variable, levels=Group.Names))) + 
    geom_boxplot() + labs(title = "TG by DB: Sum", x="", y="") +  facet_wrap(~dbs, ncol=length(unique.tg.db)) + theme_minimal () +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_fill_manual(values = colors.boxplot) + guides(fill="none")
  
  #group by saturation
  Sum.tg.bysat <- array(0, dim=c(3,Sample.No))
  unique.tg.db <- rownames(Sum.tg.bydb)
  unique.tg.db <- as.integer(gsub("dbs_","",unique.tg.db))
  if (is.element(0,unique.tg.db)) {Sum.tg.bysat[1,] = Sum.tg.bydb[which(unique.tg.db == 0),]}
  if (is.element(1,unique.tg.db)) {Sum.tg.bysat[2,] = Sum.tg.bydb[which(unique.tg.db == 1),]}
  if (length(which(unique.tg.db == 4 | unique.tg.db == 5 | unique.tg.db == 6 | unique.tg.db == 7 | unique.tg.db == 8)) > 0) {
    Sum.tg.bysat[3,] =log2(colSums(2^(Sum.tg.bydb[which(unique.tg.db == 4 | unique.tg.db == 5 | unique.tg.db == 6 | unique.tg.db == 7 | unique.tg.db == 8),])))}
  rownames(Sum.tg.bysat) <- c("Saturated","Monounsaturated","Polyunsaturated")
  colnames(Sum.tg.bysat) <- Samples
  
  df <- data.frame(Sum.tg.bysat)
  df$sat <- factor(rownames(Sum.tg.bysat),levels=c("Saturated","Monounsaturated","Polyunsaturated"))
  df.melt <- melt(df,id.vars = "sat")
  df.melt$variable <- gsub("[[:punct:]][0-9]","",df.melt$variable)
  
  p3 <- ggplot(df.melt,aes(x=factor(variable, levels=Group.Names),y=value, fill=factor(variable, levels=Group.Names))) + 
    geom_boxplot() + labs(title = "TG by sat: Sum", x="", y="") + theme_minimal () +
    facet_wrap(~sat, scales="free",ncol = 3) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = colors.boxplot) + guides(fill="none")
  
  Log2fc.tg.bysat <- array(0,dim=c(3,Sample.No))
  for (i in 1:3) {
    temp.array <- Sum.tg.bysat[i,]
    temp.mean <- mean(temp.array[1:dim(Sample.Groups)[2]])
    Log2fc.tg.bysat[i,] <- temp.array - temp.mean
  }
  colnames(Log2fc.tg.bysat) <- Samples
  rownames(Log2fc.tg.bysat) <- c("Saturated","Monounsaturated","Polyunsaturated")
  
  df <- data.frame(Log2fc.tg.bysat)
  df$sat <- factor(rownames(Sum.tg.bysat),levels=c("Saturated","Monounsaturated","Polyunsaturated"))
  df.melt <- melt(df,id.vars = "sat")
  df.melt$variable <- gsub("[[:punct:]][0-9]","",df.melt$variable)
  
  p4 <- ggplot(df.melt,aes(x=factor(variable, levels=Group.Names),y=value, fill=factor(variable, levels=Group.Names))) + 
    geom_boxplot() + labs(title = "TG by sat: Log2FC", x="", y="") + theme_minimal () +
    facet_wrap(~sat, scales="free",ncol = 3) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = colors.boxplot) + guides(fill="none")
  
  print(cowplot::plot_grid(p3,p4,nrow=2))
  
  
  ### Analyze PLs
  pl.metabolites <- Data[Lipid.indexes$PL,]$Metabolite
  pl.data <- temp.data[Lipid.indexes$PL,]
  pl.db.count <- array(0,dim=c(length(pl.metabolites),1))
  for (i in 1:length(pl.db.count)) {
    temp.metab <- unlist(strsplit(pl.metabolites[i],""))
    space.location <- c(grep(" ",temp.metab),grep("-",temp.metab))
    temp.metab <- temp.metab[1:(min(space.location)-1)]
    colon.location <- grep(":",temp.metab)
    pl.db.count[i] <- as.integer(paste(temp.metab[(colon.location+1):length(temp.metab)],collapse=""))
  }
  
  unique.pl.db <- unique(pl.db.count)
  Sum.pl.bydb <- array(0, dim=c(length(unique.pl.db),Sample.No))
  for (i in 1:length(unique.pl.db)) {
    temp.array <- pl.data[pl.db.count == unique.pl.db[i],]
    no.rows <- length(temp.array)/Sample.No
    if (no.rows > 1) {
      Sum.pl.bydb[i,] <- log2(colSums(2^(temp.array)))
    } else {
      Sum.pl.bydb[i,] <- temp.array
    }
  }
  colnames(Sum.pl.bydb) <- Samples
  rownames(Sum.pl.bydb) <- paste("dbs_",unique.pl.db,sep="")
  
  Log2fc.pl.bydb <- array(0,dim=c(length(unique.pl.db),Sample.No))
  for (i in 1:length(unique.pl.db)) {
    temp.array <- Sum.pl.bydb[i,]
    temp.mean <- mean(temp.array[1:dim(Sample.Groups)[2]])
    Log2fc.pl.bydb[i,] <- temp.array - temp.mean
  }
  colnames(Log2fc.pl.bydb) <- Samples
  rownames(Log2fc.pl.bydb) <- paste("dbs_",unique.pl.db,sep="")
  
  df <- data.frame(Log2fc.pl.bydb)
  unique.pl.db <- sort(unique.pl.db)
  df$dbs <- factor(rownames(Log2fc.pl.bydb),levels=paste("dbs_",unique.pl.db,sep=""))
  df.melt <- melt(df,id.vars = "dbs")
  df.melt$variable <- gsub("[[:punct:]][0-9]","",df.melt$variable)
  
  p2 <- ggplot(df.melt,aes(x=factor(variable, levels=Group.Names),y=value, fill=factor(variable, levels=Group.Names))) + 
    geom_boxplot() + facet_wrap(~dbs, ncol=length(unique.pl.db)) + labs(title = "PL by DB: Sum", x="", y="") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_fill_manual(values = colors.boxplot) + guides(fill="none")
  
  
  #group by saturation
  Sum.pl.bysat <- array(0, dim=c(3,Sample.No))
  unique.pl.db <- rownames(Sum.pl.bydb)
  unique.pl.db <- as.integer(gsub("dbs_","",unique.pl.db))
  if (is.element(0,unique.pl.db)) {Sum.pl.bysat[1,] = Sum.pl.bydb[which(unique.pl.db == 0),]}
  if (is.element(1,unique.pl.db)) {Sum.pl.bysat[2,] = Sum.pl.bydb[which(unique.pl.db == 1),]}
  if (length(which(unique.tg.db == 4 | unique.tg.db == 5 | unique.tg.db == 6 | unique.tg.db == 7 | unique.tg.db == 8)) > 0) {
    Sum.pl.bysat[3,] =log2(colSums(2^(Sum.pl.bydb[which(unique.tg.db == 4 | unique.tg.db == 5 | unique.tg.db == 6 | unique.tg.db == 7 | unique.tg.db == 8),])))}
  rownames(Sum.pl.bysat) <- c("Saturated","Monounsaturated","Polyunsaturated")
  colnames(Sum.pl.bysat) <- Samples
  
  df <- data.frame(Sum.pl.bysat)
  df$sat <- factor(rownames(Sum.pl.bysat),levels=c("Saturated","Monounsaturated","Polyunsaturated"))
  df.melt <- melt(df,id.vars = "sat")
  df.melt$variable <- gsub("[[:punct:]][0-9]","",df.melt$variable)
  
  p3 <- ggplot(df.melt,aes(x=factor(variable, levels=Group.Names),y=value, fill=factor(variable, levels=Group.Names))) + 
    geom_boxplot() + labs(title = "PL by sat: Sum", x="", y="") + theme_minimal () +
    facet_wrap(~sat, scales="free",ncol = 3) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = colors.boxplot) + guides(fill="none")
  
  Log2fc.pl.bysat <- array(0,dim=c(3,Sample.No))
  for (i in 1:3) {
    temp.array <- Sum.pl.bysat[i,]
    temp.mean <- mean(temp.array[1:dim(Sample.Groups)[2]])
    Log2fc.pl.bysat[i,] <- temp.array - temp.mean
  }
  colnames(Log2fc.pl.bysat) <- Samples
  rownames(Log2fc.pl.bysat) <- c("Saturated","Monounsaturated","Polyunsaturated")
  
  df <- data.frame(Log2fc.pl.bysat)
  df$sat <- factor(rownames(Sum.pl.bysat),levels=c("Saturated","Monounsaturated","Polyunsaturated"))
  df.melt <- melt(df,id.vars = "sat")
  df.melt$variable <- gsub("[[:punct:]][0-9]","",df.melt$variable)
  
  p4 <- ggplot(df.melt,aes(x=factor(variable, levels=Group.Names),y=value, fill=factor(variable, levels=Group.Names))) + 
    geom_boxplot()  + labs(title = "PL by sat: Log2FC", x="", y="") + theme_minimal () +
    facet_wrap(~sat, scales="free",ncol = 3) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = colors.boxplot) + guides(fill="none")
  
  print(cowplot::plot_grid(p3,p4,nrow=2))
  
  
  #print the plots bydb
  print(cowplot::plot_grid(p1,p2,nrow = 2))
  
  
  
  lipidomics.results <- setClass(Class = "lipidomics.results", slots=list(Sum="matrix",Log2FC="matrix",FC="matrix",
                                                                          TGSUM="ANY",TGLog2="ANY",PLSUM="ANY",PLLog2="ANY",
                                                                          TGSATSUM="ANY",TGSATLog2="ANY",PLSATSUM="ANY",PLSATLog2="ANY"))
  final.results <- lipidomics.results(Sum = Sum.byclass, Log2FC=Log2fc.byclass,FC=FC.byclass,
                                      TGSUM = Sum.tg.bydb, TGLog2 = Log2fc.tg.bydb, 
                                      PLSUM = Sum.pl.bydb, PLLog2 = Log2fc.pl.bydb,
                                      TGSATSUM = Sum.tg.bysat, TGSATLog2 = Log2fc.tg.bysat, 
                                      PLSATSUM = Sum.pl.bysat, PLSATLog2 = Log2fc.pl.bysat)
}