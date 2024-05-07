#function for calculate Sample GSZ
Sample.GSZ <- function(gene.set,ex.data){
  GSZ <- data.frame(Samples=colnames(ex.data))
  for(i in 1:length(gene.set)){
    if(length(gene.set[[i]]) > 1){
      tmp <- sqrt(length(gene.set[[i]])) * ( colMeans( na.omit(ex.data[rownames(ex.data) %in% gene.set[[i]],] )) - colMeans( ex.data ) ) / matrixStats::colSds(as.matrix(ex.data))
    }else{
      tmp <- sqrt(length(gene.set[[i]])) * ( ex.data[rownames(ex.data) %in% gene.set[[i]],] - colMeans( ex.data ) ) / matrixStats::colSds(as.matrix(ex.data))
    }
    names(tmp) <- colnames(ex.data)
    tmp <- as.data.frame(tmp)
    colnames(tmp) <- names(gene.set)[i]
    GSZ <- cbind(GSZ,tmp)
  }
  return(GSZ)
  
}
