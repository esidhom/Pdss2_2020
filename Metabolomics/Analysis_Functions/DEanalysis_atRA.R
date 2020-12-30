DEanalysis_atRA <- function(Data,Sample.No,Sample.Groups,Group.Names,FDR) {
  temp.data <- as.matrix((Data[,7:(7+Sample.No-1)]))
  row.names = vector('character')
  for (i in 1:length(Group.Names)) {row.names <- c(row.names,rep(Group.Names[i],dim(Sample.Groups)[2]))}
  
  pval.results <- array(0,dim= c(dim(temp.data)[1],1))
  log2fc.results <- array(0,dim= c(dim(temp.data)[1],1))
  for (i in 1:dim(temp.data)[1]) {
    temp.test <- temp.data[i,]
    names(temp.test) <- row.names
    
    #calc log2fcs
    log2fc.results[i,1] <- mean(temp.test[which(names(temp.test)==Group.Names[2])]) - mean(temp.test[which(names(temp.test)==Group.Names[1])])
    
    #test demso v atRA
    p.var <- var.test(temp.test[which(names(temp.test)==Group.Names[1])],
                      temp.test[which(names(temp.test)==Group.Names[2])])$p.value
    
    if (p.var > 0.05) {
      pval.results[i,1] <- t.test(temp.test[which(names(temp.test)==Group.Names[1])],
                                  temp.test[which(names(temp.test)==Group.Names[2])],var.equal = T)$p.value
    } else {
      pval.results[i,1] <- t.test(temp.test[which(names(temp.test)==Group.Names[1])],
                                  temp.test[which(names(temp.test)==Group.Names[2])], var.equal = F)$p.value
    }
  }
  
  rownames(log2fc.results) <- Data$Metabolite
  rownames(pval.results) <- Data$Metabolite
  fdr.results <- array(0,dim= c(dim(temp.data)[1],1))
  fdr.results[,1] <- p.adjust(pval.results,method="fdr")
  rownames(fdr.results) <- Data$Metabolite
  atRA.sign <- which(fdr.results<FDR)
  
  if (length(atRA.sign) > 0) {
  rownames(temp.data) <- Data$Metabolite
  pheatmap(mat =temp.data[atRA.sign,], scale = "row", fontsize = 6, 
           colorRampPalette(c("#0b3542", "white", "#6b1212"))(100), border_color=NA, 
           treeheight_col=5, main=paste("DE, B-H < ",FDR*100,"%"))
  } else {print("No significant metabolites.")}
  
  de.results <- setClass(Class = "de.results", slots=list(pvalue="ANY",fdr="ANY",log2fc="ANY"))
  final.results <- de.results(pvalue = pval.results, fdr=fdr.results,log2fc=log2fc.results)
}