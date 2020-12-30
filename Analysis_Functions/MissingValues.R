MissingValues <- function(Data, Sample.No) {
  Data.Imputed <- Data
  temp.data <- as.matrix((Data[,7:(7+Sample.No-1)]))
  
  #determine fraction of samples missing metabilite. 
  #If >30% --> delete; if <=30% --> impute missing value
  #Imputation: NA <- 1/2*min detection peak
  MissingIndexes <- vector(mode = "numeric");
  MetaboliteNo <- dim(Data)[1]
  na.Values <- c(which(is.na(temp.data)),which(temp.data == 0))
  
  
  na.Columns <- floor(na.Values/MetaboliteNo)+1;
  na.Indexes <- na.Values%%MetaboliteNo
  
  for (i in 1:length(na.Indexes)) {
    if (na.Indexes[i] == 0) {
      na.Indexes[i] <- MetaboliteNo
    }
  }
  
  UniqueIndexes <- unique(na.Indexes)
  
  for (i in 1:length(UniqueIndexes)) {
    temp.index <- UniqueIndexes[i]
    temp.count <- length(which(na.Indexes == temp.index))
    temp.fraction <- temp.count/Sample.No
    if (temp.fraction > 0.3) {
      MissingIndexes <- c(MissingIndexes,temp.index)
    } else {
      temp.array <- Data[temp.index,7:(7+Sample.No-1)]
      temp.array[is.na(temp.array)] <- 0
      temp.min <- min(temp.array)
      
      if (temp.min == 0) {
        temp.zero <- which(temp.array != 0)
        temp.min <- min(temp.array[temp.zero])
      }
      
      temp.columns <- na.Columns[which(na.Indexes == temp.index)]
      temp.array[,temp.columns] <- 1/2*temp.min
      Data.Imputed[temp.index,7:(7+Sample.No-1)] <- temp.array
    }
  }
  print('Missing Metabolites:')
  Missing.Metabolites <- Data$Metabolite[MissingIndexes]
  Missing.Metabolites <- Missing.Metabolites[which(Missing.Metabolites != '')]
  print(unique(Missing.Metabolites))
  Data.Imputed <- Data.Imputed[-MissingIndexes,]
  
  #check for successful elimination of missing values
  temp.data <- as.matrix(Data.Imputed[,7:(7+Sample.No-1)])
  na.Values <- which(is.na(temp.data))
  if (length(na.Values) == 0) {print('Successful removal of all missing values.')}
  
  return(Data.Imputed)
}