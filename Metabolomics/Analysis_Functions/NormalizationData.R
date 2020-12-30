NormalizationData <- function(Data,Methods,Sample.No) {
  Data.Norm <- Data
  Std.Indexes <- grep("Standard",Data$HMDB.ID...representativeID.)
  
  #remove internal standards
  if (length(Std.Indexes) > 0) {Data.Norm <- Data.Norm[-Std.Indexes,]}
  
  
  temp.data <- as.matrix((Data.Norm[,7:(7+Sample.No-1)]))
  
  for (i in 1:length(Methods)) {
    temp.method <- Methods[i]
    temp.indexes <- which(Data.Norm$Method == temp.method)
    print(temp.method)
    print(paste('Number of metabolites:', length(temp.indexes)))
    
    temp.array <- temp.data[temp.indexes,]
    
    #log2transformation
    temp.array <- log2(temp.array)
    
    #mean-centering
    temp.mean <- colMeans(temp.array)
    temp.array <- temp.array - temp.mean
    
    temp.data[temp.indexes,] <- temp.array
    Data.Norm[,7:(7+Sample.No-1)] <- temp.data
  }
  return(Data.Norm)
}