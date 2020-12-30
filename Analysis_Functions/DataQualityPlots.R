DataQualityPlots <- function(Data.Imputed,Data.Norm,Methods,MethodList,Samples,Sample.No,Sample.Groups,Group.Names) {
  data.imputed <- as.matrix((Data.Imputed[,7:(7+Sample.No-1)]))
  data.norm <- as.matrix((Data.Norm[,7:(7+Sample.No-1)]))
  
  
  ########## Make data quality plots for raw imputed data #######
  ###plot sum/mean of all metabolites per method per sample
  sum.array <- array(0, dim = c(length(Methods),length(Samples)))
  colnames(sum.array) <- Samples
  mean.array <- array(0, dim = c(length(Methods),length(Samples)))
  colnames(mean.array) <- Samples
  for (i in 1:length(Methods)){
    temp.method <- Methods[i]
    temp.indexes <- which(Data.Imputed$Method == temp.method)
    temp.array <- data.imputed[temp.indexes,]
    sum.array[i,] <- colSums(temp.array)
    mean.array[i,] <- colMeans(temp.array)
  }
  sum.df <- data.frame(Methods,sum.array)
  mean.df <- data.frame(Methods,mean.array)
  
  sum.df <- melt(sum.df,id.vars = 'Methods')
  mean.df <- melt(mean.df,id.vars = 'Methods')
  
  p <- ggplot(sum.df) + geom_bar(stat="identity",aes(x=variable,y=value)) + 
    facet_wrap(~Methods, scales="free") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="Sum", x="", y="")
  print(p)
  
  p <- ggplot(mean.df) + geom_bar(stat="identity",aes(x=variable,y=value)) + 
    facet_wrap(~Methods, scales="free") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="Mean", x="", y="") 
  print(p)
  
  ###Plot histogram of samples by methods
  hist.df <- melt(data.frame(Data.Imputed$Method,log10(data.imputed)))
  p <- ggplot(hist.df) + facet_wrap(~Data.Imputed.Method) + 
    geom_histogram(aes(x=value,fill=variable), 
                   alpha=0.2, 
                   position="identity") +
    theme_minimal() + guides(alpha="none") +
    labs(x="log10(Abundance)")
  print(p)
  
  ###plot std vs raw mean peak area on log-scale
  temp.std <- log10(apply(data.imputed, 1, sd))
  temp.mean <- log10(apply(data.imputed,1,mean))
  df <- data.frame(temp.std,temp.mean)
  p <- ggplot(df, aes(x=temp.mean,y=temp.std)) + geom_point(aes(alpha=0.1)) +
    guides(alpha="none") + theme_minimal() + labs(title="", x="log10(Mean)", y="log10(SD)")
  print(p)
  
  
  ########## Make data quality plots for normalized imputed data #######
  ###plot sum/mean of all metabolites per method per sample
  sum.array <- array(0, dim = c(length(Methods),length(Samples)))
  colnames(sum.array) <- Samples
  mean.array <- array(0, dim = c(length(Methods),length(Samples)))
  colnames(mean.array) <- Samples
  for (i in 1:length(Methods)){
    temp.method <- Methods[i]
    temp.indexes <- which(Data.Norm$Method == temp.method)
    temp.array <- data.norm[temp.indexes,]
    sum.array[i,] <- colSums(temp.array)
    mean.array[i,] <- colMeans(temp.array)
  }
  sum.df <- data.frame(Methods,sum.array)
  mean.df <- data.frame(Methods,mean.array)
  
  sum.df <- melt(sum.df,id.vars = 'Methods')
  mean.df <- melt(mean.df,id.vars = 'Methods')
  
  p <- ggplot(sum.df) + geom_bar(stat="identity",aes(x=variable,y=value)) + 
    facet_wrap(~Methods, scales="free") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="Sum", x="", y="") 
  print(p)
  
  p <- ggplot(mean.df) + geom_bar(stat="identity",aes(x=variable,y=value)) + 
    facet_wrap(~Methods, scales="free") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="Mean", x="", y="") 
  print(p)
  
  ###Plot histogram of samples by methods
  hist.df <- melt(data.frame(Data.Norm$Method,data.norm))
  p <- ggplot(hist.df) + facet_wrap(~Data.Norm.Method) + 
    geom_histogram(aes(x=value,fill=variable), 
                   alpha=0.2, 
                   position="identity") +
    theme_minimal() + guides(alpha="none") +
    labs(x="log2(Abundance), mean-centered")
  print(p)
  
  ###plot std vs raw mean peak area on log-scale
  temp.std <- apply(data.norm, 1, sd)
  temp.mean <- apply(data.norm,1,mean)
  df <- data.frame(temp.std,temp.mean)
  p <- ggplot(df, aes(x=temp.mean,y=temp.std)) + geom_point(aes(alpha=0.1)) +
    guides(alpha="none") + theme_minimal() + labs(title="", x="norm. mean", y="norm. SD") 
  print(p)
}