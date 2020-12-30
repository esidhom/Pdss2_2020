PCAanalysis <- function(Data,pc1,pc2,Samples,Sample.No,Sample.Groups,Group.Names) {
  temp.data <- as.matrix((Data[,7:(7+Sample.No-1)]))
  pca.data <- prcomp(temp.data,center=F,scale=T)
  
  variance.explained <- ((pca.data$sdev)^2/sum((pca.data$sdev)^2))*100
  variance.explained <- round(variance.explained,3)
  print('%Variance Explained')
  print(variance.explained)
  
  pca.loadings <- as.data.frame(pca.data$rotation)
  row.names = vector('character')
  for (i in 1:length(Group.Names)) {row.names <- c(row.names,rep(Group.Names[i],dim(Sample.Groups)[2]))}
  pca.loadings$sample.groups <- row.names
  p<- ggplot(pca.loadings,aes(x=eval(as.symbol(paste('PC',pc1,sep=""))),
                              y=eval(as.symbol(paste('PC',pc2,sep=""))),
                              color=sample.groups))+
    geom_point() + labs(x=paste('PC',pc1,": ", variance.explained[pc1],"% variance explained",sep=""),
                        y=paste('PC',pc2,": ", variance.explained[pc2],"% variance explained",sep=""))
  print(p)
  return(pca.data)
  
}