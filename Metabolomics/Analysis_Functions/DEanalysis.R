DEanalysis <- function(Data,Sample.No,Sample.Groups,Group.Names,FDR) {
  temp.data <- as.matrix((Data[,7:(7+Sample.No-1)]))
  row.names = vector('character')
  for (i in 1:length(Group.Names)) {row.names <- c(row.names,rep(Group.Names[i],dim(Sample.Groups)[2]))}
  
  pval.results <- array(0,dim= c(dim(temp.data)[1],3))
  log2fc.results <- array(0,dim= c(dim(temp.data)[1],2))
  for (i in 1:dim(temp.data)[1]) {
    temp.test <- temp.data[i,]
    names(temp.test) <- row.names
    
    #calc log2fcs
    log2fc.results[i,1] <- mean(temp.test[which(names(temp.test)==Group.Names[2])]) - mean(temp.test[which(names(temp.test)==Group.Names[1])])
    log2fc.results[i,2] <- mean(temp.test[which(names(temp.test)==Group.Names[3])]) - mean(temp.test[which(names(temp.test)==Group.Names[1])])
    
    #test scr v p1
    p.var <- var.test(temp.test[which(names(temp.test)==Group.Names[1])],
                      temp.test[which(names(temp.test)==Group.Names[2])])$p.value
    
    if (p.var > 0.05) {
      pval.results[i,1] <- t.test(temp.test[which(names(temp.test)==Group.Names[1])],
                                  temp.test[which(names(temp.test)==Group.Names[2])],var.equal = T)$p.value
    } else {
      pval.results[i,1] <- t.test(temp.test[which(names(temp.test)==Group.Names[1])],
                                  temp.test[which(names(temp.test)==Group.Names[2])], var.equal = F)$p.value
    }
    
    #test scr v p2
    p.var <- var.test(temp.test[which(names(temp.test)==Group.Names[1])],
                      temp.test[which(names(temp.test)==Group.Names[3])])$p.value
    
    if (p.var > 0.05) {
      pval.results[i,2] <- t.test(temp.test[which(names(temp.test)==Group.Names[1])],
                                  temp.test[which(names(temp.test)==Group.Names[3])],var.equal = T)$p.value
    } else {
      pval.results[i,2] <- t.test(temp.test[which(names(temp.test)==Group.Names[1])],
                                  temp.test[which(names(temp.test)==Group.Names[3])], var.equal = F)$p.value
    }
    
    #test p1 v p2
    p.var <- var.test(temp.test[which(names(temp.test)==Group.Names[2])],
                      temp.test[which(names(temp.test)==Group.Names[3])])$p.value
    
    if (p.var > 0.05) {
      pval.results[i,3] <- t.test(temp.test[which(names(temp.test)==Group.Names[2])],
                                  temp.test[which(names(temp.test)==Group.Names[3])],var.equal = T)$p.value
    } else {
      pval.results[i,3] <- t.test(temp.test[which(names(temp.test)==Group.Names[2])],
                                  temp.test[which(names(temp.test)==Group.Names[3])], var.equal = F)$p.value
    }
  }
  
  colnames(log2fc.results) <- c("Scr.v.P1","Scr.v.P2")
  rownames(log2fc.results) <- Data$Metabolite
  log2fc.combined <- sign(log2fc.results[,1])*sign(log2fc.results[,2])
  
  colnames(pval.results) <- c("Scr.v.P1","Scr.v.P2","P1.v.P2")
  rownames(pval.results) <- Data$Metabolite
  fdr.results <- apply(pval.results,2,function(x) p.adjust(x,method="fdr"))
  colnames(fdr.results) <- c("Scr.v.P1","Scr.v.P2","P1.v.P2")
  rownames(fdr.results) <- Data$Metabolite
  scrp1.sign <- which(fdr.results[,1]<FDR)
  scrp2.sign <- which(fdr.results[,2]<FDR)
  p1p2.nonsign <- which(fdr.results[,3]>FDR)
  p1p2.samedir <- which(log2fc.combined == 1)
  
  p1p2.sign <- intersect(scrp1.sign,scrp2.sign)
  p1p2.sign <- intersect(intersect(scrp1.sign,scrp2.sign),p1p2.samedir)
  
  rownames(temp.data) <- Data$Metabolite
  pheatmap(mat =temp.data[p1p2.sign,], scale = "row", fontsize = 6, 
           colorRampPalette(c("navy", "white", "firebrick3"))(100), border_color=NA, 
           treeheight_col=5, main=paste("DE, B-H < ",FDR*100,"%"))
  de.results <- setClass(Class = "de.results", slots=list(pvalue="matrix",fdr="matrix",log2fc="matrix"))
  
  #browser()
  final.results <- de.results(pvalue = pval.results, fdr=fdr.results,log2fc=log2fc.results)
}